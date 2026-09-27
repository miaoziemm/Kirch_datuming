#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <chrono>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

static int str_ieq(const char *a, const char *b)
{
    if (a == NULL || b == NULL) return 0;
    while (*a != '\0' && *b != '\0') {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return 0;
        a++;
        b++;
    }
    return (*a == '\0' && *b == '\0');
}

static inline int idx2d(int iz, int ix, int nz)
{
    return ix * nz + iz;
}

static inline float clampf_local(float v, float lo, float hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

static int value_to_nearest_bin(float value, float g0, float dg, int ng)
{
    if (ng <= 0 || dg <= 0.0f || !isfinite(value)) return -1;
    int ig = (int)floorf((value - g0) / dg + 0.5f);
    if (ig < 0 || ig >= ng) return -1;
    return ig;
}

/*
 * Distribute one continuous angle value to the two neighboring angle bins.
 *
 * Example:
 *   angle = 25.4 deg, g0 = 0 deg, dg = 1 deg
 * gives
 *   ig0 = 25, w0 = 0.6
 *   ig1 = 26, w1 = 0.4
 *
 * At the first/last output sample the whole contribution is kept in the
 * boundary bin, so no energy is lost only because of floating-point roundoff.
 */
static int value_to_linear_bins(float value, float g0, float dg, int ng,
                                int *ig0, int *ig1, float *w0, float *w1)
{
    if (ng <= 0 || dg <= 0.0f || !isfinite(value)) return 0;

    const float u = (value - g0) / dg;
    if (u < -1.0e-5f || u > (float)(ng - 1) + 1.0e-5f) return 0;

    if (u <= 0.0f) {
        *ig0 = 0;
        *ig1 = 0;
        *w0 = 1.0f;
        *w1 = 0.0f;
        return 1;
    }

    if (u >= (float)(ng - 1)) {
        *ig0 = ng - 1;
        *ig1 = ng - 1;
        *w0 = 1.0f;
        *w1 = 0.0f;
        return 1;
    }

    *ig0 = (int)floorf(u);
    *ig1 = *ig0 + 1;
    *w1 = u - (float)(*ig0);
    *w0 = 1.0f - *w1;
    return 1;
}

/*
 * Optional cosine taper for large reflection angles.
 *
 * start < 0 or end <= start : disabled
 * angle <= start            : weight = 1
 * start < angle < end       : cosine taper
 * angle >= end              : weight = 0
 */
static float angle_taper_weight(float angle, float start, float end)
{
    if (!isfinite(angle)) return 0.0f;
    if (start < 0.0f || end <= start) return 1.0f;
    if (angle <= start) return 1.0f;
    if (angle >= end) return 0.0f;

    const float x = (angle - start) / (end - start);
    return 0.5f * (1.0f + cosf((float)M_PI * x));
}

/*
 * Mute the early-time fan that opens to both sides of a shot.
 *
 * The mute boundary for a receiver at xr is
 *
 *     normal:   t_boundary = intercept + slope * |xr - xs|
 *     inverted: t_boundary = intercept - slope * |xr - xs|.
 *
 * Samples before the boundary are removed.  A raised-cosine ramp after the
 * boundary avoids the high-frequency noise caused by a hard time cut.  The
 * routine operates on the input trace before the integrations used by the
 * Kirchhoff kernel, so muted energy cannot leak back into the angle gather.
 */
static void apply_direct_wave_mute(float *trace, int nt, float t0, float dt,
                                   float source_x, float receiver_x,
                                   float slope, float intercept,
                                   float taper_length, int inverted)
{
    const float offset_time = slope * fabsf(receiver_x - source_x);
    const float boundary = inverted ? intercept - offset_time
                                    : intercept + offset_time;

    for (int it = 0; it < nt; it++) {
        const float time = t0 + (float)it * dt;
        float weight = 0.0f;

        if (time >= boundary + taper_length) {
            weight = 1.0f;
        } else if (time > boundary) {
            const float u = (time - boundary) / taper_length;
            weight = 0.5f * (1.0f - cosf((float)M_PI * u));
        }

        trace[it] *= weight;
    }
}

static float finite_diff_x_radius(const float *tab, int nz, int nx,
                                  int iz, int ix, float dx, int radius)
{
    if (nx <= 1 || dx == 0.0f) return 0.0f;
    if (radius < 1) radius = 1;

    const int il = MAX(0, ix - radius);
    const int ir = MIN(nx - 1, ix + radius);

    if (ir == il) return 0.0f;

    return (tab[idx2d(iz, ir, nz)] - tab[idx2d(iz, il, nz)]) /
           ((float)(ir - il) * dx);
}

static float finite_diff_z_radius(const float *tab, int nz, int nx,
                                  int iz, int ix, float dz, int radius)
{
    (void)nx;
    if (nz <= 1 || dz == 0.0f) return 0.0f;
    if (radius < 1) radius = 1;

    const int it = MAX(0, iz - radius);
    const int ib = MIN(nz - 1, iz + radius);

    if (ib == it) return 0.0f;

    return (tab[idx2d(ib, ix, nz)] - tab[idx2d(it, ix, nz)]) /
           ((float)(ib - it) * dz);
}

/*
 * A radius larger than one suppresses grid-scale fluctuations in the
 * traveltime gradient.  This is especially useful for Marmousi, where
 * first-arrival traveltime tables can contain locally rapid changes in
 * gradient direction.  The traveltime itself is not smoothed; only the
 * gradient used for angle classification is stabilized.
 */
static float half_opening_angle_deg_from_grad(float sx, float sz, float rx, float rz)
{
    float ns = sqrtf(sx * sx + sz * sz);
    float nr = sqrtf(rx * rx + rz * rz);
    if (ns <= 1.0e-12f || nr <= 1.0e-12f) return NAN;

    /*
       Here sx,sz and rx,rz are the source- and receiver-traveltime
       gradients at the image point. The half-opening angle is

           alpha = 0.5 * angle(grad(Ts), grad(Tr)).

       atan2(|cross|, dot) is mathematically equivalent to acos(dot/(|a||b|))
       but is numerically better behaved when the opening angle is small.
    */
    const float dot = sx * rx + sz * rz;
    const float cross = sx * rz - sz * rx;
    const float opening = atan2f(fabsf(cross), dot);

    return 0.5f * opening * 180.0f / (float)M_PI;
}

static float half_opening_angle_deg_from_tables(const float *stable,
                                                const float *rtable,
                                                int nz, int nx,
                                                int iz, int ix,
                                                float dx, float dz,
                                                int radius)
{
    float sx = finite_diff_x_radius(stable, nz, nx, iz, ix, dx, radius);
    float sz = finite_diff_z_radius(stable, nz, nx, iz, ix, dz, radius);
    float rx = finite_diff_x_radius(rtable, nz, nx, iz, ix, dx, radius);
    float rz = finite_diff_z_radius(rtable, nz, nx, iz, ix, dz, radius);
    return half_opening_angle_deg_from_grad(sx, sz, rx, rz);
}

static void copy_n(float *dst, const float *src, off_t n)
{
    for (off_t i = 0; i < n; i++) dst[i] = src[i];
}

static void load_table_slice(sep_t *fp, int index, float *buf, off_t nxy,
                             int *cached_index, const char *name)
{
    if (cached_index != NULL && *cached_index == index) return;
    const off_t byte_off = (off_t)index * nxy * (off_t)sizeof(float);
    if (se_fsio_seek(fp->data->io, byte_off) != CODE_SUCCESS)
        ERROR(("Failed to seek %s slice %d", name, index));
    if (se_fsio_read_float(fp->data->io, buf, (size_t)nxy) != CODE_SUCCESS)
        ERROR(("Failed to read %s slice %d", name, index));
    if (cached_index != NULL) *cached_index = index;
}

/* Preserve the original volume interpolation exactly while streaming slices. */
static void interp_volume_linear_stream(float *out, sep_t *fp, int npos,
                                        off_t nxy, float coord, float o,
                                        float d, float *slice0, float *slice1,
                                        const char *name)
{
    int i0 = 0;

    if (npos > 1 && d != 0.0f) {
        const float f = (coord - o) / d;
        i0 = (int)floorf(f);

        if (i0 > 0 && i0 < npos - 1) {
            const float w = f - (float)i0;
            load_table_slice(fp, i0, slice0, nxy, NULL, name);
            load_table_slice(fp, i0 + 1, slice1, nxy, NULL, name);
            for (off_t i = 0; i < nxy; i++)
                out[i] = (1.0f - w) * slice0[i] + w * slice1[i];
            return;
        }
    }

    const int index = i0 >= npos - 1 ? npos - 1 : 0;
    load_table_slice(fp, index, slice0, nxy, NULL, name);
    copy_n(out, slice0, nxy);
}

static int parse_output_axis(int is_angle, int nh,
                             float h0, float dh,
                             float *g0, float *dg, int *ng)
{
    float hend = h0 + (float)(nh - 1) * dh;
    float vmin = is_angle ? 0.0f : MIN(h0, hend);
    float vmax = is_angle ? 90.0f : MAX(h0, hend);
    float vstep = is_angle ? 1.0f : (dh != 0.0f ? fabsf(dh) : 1.0f);
    int n_given = 0;
    int nval = is_angle ? 91 : nh;

    if (se_have_par("cig_min")) vmin = se_get_par_float("cig_min");
    if (se_have_par("cig_max")) vmax = se_get_par_float("cig_max");
    if (se_have_par("cig_step")) vstep = se_get_par_float("cig_step");
    if (se_have_par("cig_n")) {
        nval = se_get_par_int("cig_n");
        n_given = 1;
    }

    if (is_angle) {
        if (se_have_par("angle_min")) vmin = se_get_par_float("angle_min");
        if (se_have_par("angle_max")) vmax = se_get_par_float("angle_max");
        if (se_have_par("angle_step")) vstep = se_get_par_float("angle_step");
        if (se_have_par("angle_n")) {
            nval = se_get_par_int("angle_n");
            n_given = 1;
        }
    } else {
        if (se_have_par("offset_min")) vmin = se_get_par_float("offset_min");
        if (se_have_par("offset_max")) vmax = se_get_par_float("offset_max");
        if (se_have_par("offset_step")) vstep = se_get_par_float("offset_step");
        if (se_have_par("offset_n")) {
            nval = se_get_par_int("offset_n");
            n_given = 1;
        }
    }

    if (vstep <= 0.0f) ERROR(("CIG axis step must be positive."));

    if (!n_given) {
        nval = (int)floorf((vmax - vmin) / vstep + 0.5f) + 1;
    }
    if (nval <= 0) ERROR(("CIG axis sample number must be positive."));

    *g0 = vmin;
    *dg = vstep;
    *ng = nval;
    return 0;
}

static sep_t *open_grad_volume(const char *name, char *file,
                               int nz, int nx, int n3)
{
    sep_t *fp = sep_open(file, SEP_READ, 0);
    if (fp->headers->ndim < 3) {
        ERROR(("Need 3D data volume for %s= (z,x,source_or_receiver)", name));
    }
    if (fp->headers->n[0] != nz || fp->headers->n[1] != nx || fp->headers->n[2] != n3) {
        ERROR(("Dimension mismatch in %s=. Expected n1=%d n2=%d n3=%d", name, nz, nx, n3));
    }
    return fp;
}

void diff2(float *trace, int n, float d)
{
    int i;
    float *tmp = alloc1float(n);
    for (i = 0; i < n; i++) {
        tmp[i] = trace[i];
    }
    for (i = 1; i < n - 1; i++) {
        trace[i] = (tmp[i + 1] - tmp[i - 1]) / (2.0f * d);
    }
    trace[0] = (tmp[1] - tmp[0]) / d;
    trace[n - 1] = (tmp[n - 1] - tmp[n - 2]) / d;
    free1float(tmp);
}

int main(int argc, char *argv[])
{
    se_par_init(argc, argv);

    char *unit;
    const char *type = NULL;
    const char *mode = NULL;
    int isdiff;
    int adj = 0, cig = 0, cmp = 0;
    int mode_given = 0;
    int cig_is_angle = 0;
    int use_input_angle_grad = 0;
    off_t nzx = 0;
    int nt = 0, nx = 0, sny = 0, rny = 0, ns = 0, nh = 0, nz = 0;
    int i = 0, ix = 0, iz = 0, ih = 0, is = 0, ist = 0, iht = 0, ng = 0, ithr = 0, nthr = 0;
    int ig = 0, ig_offset = 0;
    float *trace = NULL, **traces = NULL, **out = NULL;
    float *stable = NULL, *rtable = NULL, *stablex = NULL, *rtablex = NULL;
    float *stbl0 = NULL, *stbl1 = NULL, *stblx0 = NULL, *stblx1 = NULL;
    float *rtbl0 = NULL, *rtbl1 = NULL, *rtblx0 = NULL, *rtblx1 = NULL;
    float *grad_slice0 = NULL, *grad_slice1 = NULL;
    int stbl0_idx = -1, stbl1_idx = -1, stblx0_idx = -1, stblx1_idx = -1;
    int rtbl0_idx = -1, rtbl1_idx = -1, rtblx0_idx = -1, rtblx1_idx = -1;
    float *sgradx = NULL, *sgradz = NULL, *rgradx = NULL, *rgradz = NULL;
    float ds = 0.0f, s0 = 0.0f, x0 = 0.0f, sy0 = 0.0f, sdy = 0.0f, ry0 = 0.0f, rdy = 0.0f;
    float s = 0.0f, h = 0.0f, h0 = 0.0f, dh = 0.0f, dx = 0.0f;
    float ti = 0.0f, t0 = 0.0f, t1 = 0.0f, t2 = 0.0f, dt = 0.0f, z0 = 0.0f, dz = 0.0f, tau = 0.0f;
    float aal = 0.0f, tx = 0.0f, aper = 0.0f;
    float cig0 = 0.0f, dcig = 1.0f, angle = 0.0f;
    float angle_taper_start = -1.0f, angle_taper_end = -1.0f;
    int angle_interp = 1;
    int angle_grad_radius = 2;
    int aperture_trace_taper = 20;
    int direct_mute = 0;
    int direct_mute_invert = 0;
    int direct_mute_qc_shot = -1;
    float direct_mute_slope = 0.0f;
    float direct_mute_intercept = 0.0f;
    float direct_mute_taper = 0.05f;
    sep_t *dat = NULL, *mig = NULL, *stim = NULL, *sder = NULL, *rtim = NULL, *rder = NULL;
    sep_t *mute_before = NULL, *mute_after = NULL;
    sep_t *sgradx_fp = NULL, *sgradz_fp = NULL, *rgradx_fp = NULL, *rgradz_fp = NULL;
    char *dat_f = NULL, *mig_f = NULL, *stim_f = NULL, *sder_f = NULL, *rtim_f = NULL, *rder_f = NULL;
    char *sgradx_f = NULL, *sgradz_f = NULL, *rgradx_f = NULL, *rgradz_f = NULL;
    char *mute_before_f = NULL, *mute_after_f = NULL;

    int aperture_trace = -1; /* limit the number of traces used for imaging; -1 means no limit */

    if (!se_have_par("aperture_trace"))
        aperture_trace = -1;
    else
        aperture_trace = se_get_par_int("aperture_trace");

    INFO(("aperture_trace = %d\n", aperture_trace));

    mode_given = se_have_par("mode") || se_have_par("cig_mode");
    if (se_have_par("mode"))
        mode = se_get_par_str("mode");
    else if (se_have_par("cig_mode"))
        mode = se_get_par_str("cig_mode");
    else
        mode = "offset";

    if (str_ieq(mode, "offset")) {
        cig_is_angle = 0;
    } else if (str_ieq(mode, "angle")) {
        cig_is_angle = 1;
    } else {
        ERROR(("Unknown CIG mode: %s. Supported modes are offset and angle.", mode));
    }

    if (!se_have_par("adj"))
        adj = 1;
    else
        adj = se_get_par_int("adj");
    if (!se_have_par("cig"))
        cig = mode_given ? 1 : 0;
    else
        cig = se_get_par_int("cig");
    if (!se_have_par("cmp"))
        cmp = 1;
    else
        cmp = se_get_par_int("cmp");

    if (!se_have_par("diff"))
        isdiff = 0;
    else
        isdiff = se_get_par_int("diff");

    if (!se_have_par("seismic_data"))
        ERROR(("Need seismic_data="));
    else
        dat_f = se_get_par_str("seismic_data");
    if (!se_have_par("migration"))
        ERROR(("Need migration="));
    else
        mig_f = se_get_par_str("migration");
    if (!se_have_par("stable"))
        ERROR(("Need stable="));
    else
        stim_f = se_get_par_str("stable");
    if (!se_have_par("sderiv"))
        ERROR(("Need sderiv="));
    else
        sder_f = se_get_par_str("sderiv");
    if (!se_have_par("rtable"))
        ERROR(("Need rtable="));
    else
        rtim_f = se_get_par_str("rtable");
    if (!se_have_par("rderiv"))
        ERROR(("Need rderiv="));
    else
        rder_f = se_get_par_str("rderiv");

    if (cig && cig_is_angle) {
        int has_sgx = se_have_par("sgradx");
        int has_sgz = se_have_par("sgradz");
        int has_rgx = se_have_par("rgradx");
        int has_rgz = se_have_par("rgradz");
        if (has_sgx || has_sgz || has_rgx || has_rgz) {
            if (!(has_sgx && has_sgz && has_rgx && has_rgz)) {
                ERROR(("For mode=angle, please provide all of sgradx=, sgradz=, rgradx=, and rgradz=, or provide none of them."));
            }
            use_input_angle_grad = 1;
            sgradx_f = se_get_par_str("sgradx");
            sgradz_f = se_get_par_str("sgradz");
            rgradx_f = se_get_par_str("rgradx");
            rgradz_f = se_get_par_str("rgradz");
        }
    }

    if (adj) {
        dat = sep_open(dat_f, SEP_READ, 0);
        mig = sep_open(mig_f, SEP_WRITE, 0);
    } else {
        mig = sep_open(mig_f, SEP_READ, 0);
        dat = sep_open(dat_f, SEP_WRITE, 0);
    }

    if (adj) {
        if (dat->headers->ndim < 3)
            ERROR(("Need 3D data volume for seismic_data= (t,h/r,s)"));
        nt = dat->headers->n[0]; /* time samples */
        nh = dat->headers->n[1]; /* offset samples */
        ns = dat->headers->n[2]; /* shot samples */

        t0 = (float)dat->headers->o[0]; /* time origin */
        dt = (float)dat->headers->d[0]; /* time sampling */
        h0 = (float)dat->headers->o[1]; /* offset origin */
        dh = (float)dat->headers->d[1]; /* offset sampling */
        s0 = (float)dat->headers->o[2]; /* shot origin */
        ds = (float)dat->headers->d[2]; /* shot sampling */
    } else {
        if (!se_have_par("nt"))
            ERROR(("Need nt="));
        else
            nt = se_get_par_int("nt");
        if (!se_have_par("nh"))
            nh = 1;
        else
            nh = se_get_par_int("nh");
        if (!se_have_par("ns"))
            ns = 1;
        else
            ns = se_get_par_int("ns");
        if (!se_have_par("t0"))
            t0 = 0.0f;
        else
            t0 = se_get_par_float("t0");
        if (!se_have_par("dt"))
            ERROR(("Need dt="));
        else
            dt = se_get_par_float("dt");
        if (!se_have_par("h0"))
            h0 = 0.0f;
        else
            h0 = se_get_par_float("h0");
        if (!se_have_par("dh"))
            ERROR(("Need dh="));
        else
            dh = se_get_par_float("dh");
        if (!se_have_par("s0"))
            s0 = 0.0f;
        else
            s0 = se_get_par_float("s0");
        if (!se_have_par("ds"))
            ERROR(("Need ds="));
        else
            ds = se_get_par_float("ds");
    }

    if (1 == nh) dh = 0.0f;
    if (1 == ns) ds = 0.0f;

    /*
     * Optional smooth direct-wave/side-noise mute for input shot gathers.
     * Receiver position is s+h for CMP/offset input and h for receiver input.
     * direct_mute_slope has units of seconds per lateral-coordinate unit.
     */
    if (se_have_par("direct_mute"))
        direct_mute = se_get_par_int("direct_mute");
    if (se_have_par("direct_mute_slope"))
        direct_mute_slope = se_get_par_float("direct_mute_slope");
    if (se_have_par("direct_mute_invert"))
        direct_mute_invert = se_get_par_int("direct_mute_invert");
    if (se_have_par("direct_mute_intercept"))
        direct_mute_intercept = se_get_par_float("direct_mute_intercept");
    if (se_have_par("direct_mute_taper"))
        direct_mute_taper = se_get_par_float("direct_mute_taper");
    if (se_have_par("direct_mute_qc_shot"))
        direct_mute_qc_shot = se_get_par_int("direct_mute_qc_shot");

    if (direct_mute) {
        if (!adj)
            ERROR(("direct_mute=1 is an input-data preprocessing option and requires adj=1."));
        if (!isfinite(direct_mute_slope) || direct_mute_slope < 0.0f)
            ERROR(("direct_mute_slope must be finite and non-negative."));
        if (!isfinite(direct_mute_intercept))
            ERROR(("direct_mute_intercept must be finite."));
        if (!isfinite(direct_mute_taper) || direct_mute_taper <= 0.0f)
            ERROR(("direct_mute_taper must be finite and positive to avoid a hard cutoff."));
        if (direct_mute_invert != 0 && direct_mute_invert != 1)
            ERROR(("direct_mute_invert must be 0 or 1."));
        if (direct_mute_qc_shot < -1)
            ERROR(("direct_mute_qc_shot must be -1 (disabled) or a zero-based shot index."));
        if (direct_mute_qc_shot >= ns)
            ERROR(("direct_mute_qc_shot=%d is outside the input shot range [0,%d].",
                   direct_mute_qc_shot, ns - 1));
        if (direct_mute_qc_shot >= 0) {
            if (!se_have_par("direct_mute_qc_before") ||
                !se_have_par("direct_mute_qc_after"))
                ERROR(("direct_mute_qc_shot requires direct_mute_qc_before= and direct_mute_qc_after=."));
            mute_before_f = se_get_par_str("direct_mute_qc_before");
            mute_after_f = se_get_par_str("direct_mute_qc_after");
            mute_before = sep_open(mute_before_f, SEP_WRITE, 0);
            mute_after = sep_open(mute_after_f, SEP_WRITE, 0);
            sep_t *qc_files[2] = {mute_before, mute_after};
            for (int iqc = 0; iqc < 2; iqc++) {
                qc_files[iqc]->headers->ndim = 2;
                qc_files[iqc]->headers->n[0] = nt;
                qc_files[iqc]->headers->n[1] = nh;
                qc_files[iqc]->headers->o[0] = t0;
                qc_files[iqc]->headers->d[0] = dt;
                qc_files[iqc]->headers->o[1] = h0;
                qc_files[iqc]->headers->d[1] = dh;
                sep_set_header(qc_files[iqc], "label1", "Time");
                sep_set_header(qc_files[iqc], "unit1", "s");
                sep_set_header(qc_files[iqc], "label2", cmp ? "Offset" : "Receiver");
            }
            INFO(("direct mute QC: write zero-based shot %d before/after mute to %s and %s\n",
                  direct_mute_qc_shot, mute_before_f, mute_after_f));
        }

        INFO(("direct mute enabled: boundary = %g %c %g*abs(receiver-source) s, cosine taper = %g s\n",
              direct_mute_intercept, direct_mute_invert ? '-' : '+',
              direct_mute_slope, direct_mute_taper));
    }

    stim = sep_open(stim_f, SEP_READ, 0);
    sder = sep_open(sder_f, SEP_READ, 0);

    if (stim->headers->ndim < 3)
        ERROR(("Need 3D data volume for stable= (z,x,s)"));
    if (sder->headers->ndim < 3)
        ERROR(("Need 3D data volume for sderiv= (z,x,s)"));

    nz = stim->headers->n[0];
    nx = stim->headers->n[1];
    sny = stim->headers->n[2];

    z0 = (float)stim->headers->o[0];
    dz = (float)stim->headers->d[0];
    x0 = (float)stim->headers->o[1];
    dx = (float)stim->headers->d[1];
    sy0 = (float)stim->headers->o[2];
    sdy = (float)stim->headers->d[2];

    nzx = (off_t)nz * (off_t)nx;

    rtim = sep_open(rtim_f, SEP_READ, 0);
    rder = sep_open(rder_f, SEP_READ, 0);

    if (rtim->headers->ndim < 3)
        ERROR(("Need 3D data volume for rtable= (z,x,r)"));
    if (rder->headers->ndim < 3)
        ERROR(("Need 3D data volume for rderiv= (z,x,r)"));

    rny = rtim->headers->n[2];
    ry0 = (float)rtim->headers->o[2];
    rdy = (float)rtim->headers->d[2];

    if (use_input_angle_grad) {
        sgradx_fp = open_grad_volume("sgradx", sgradx_f, nz, nx, sny);
        sgradz_fp = open_grad_volume("sgradz", sgradz_f, nz, nx, sny);
        rgradx_fp = open_grad_volume("rgradx", rgradx_f, nz, nx, rny);
        rgradz_fp = open_grad_volume("rgradz", rgradz_f, nz, nx, rny);
        sgradx = alloc1float(nzx);
        sgradz = alloc1float(nzx);
        rgradx = alloc1float(nzx);
        rgradz = alloc1float(nzx);
        INFO(("mode=angle: stream input spatial traveltime gradients sgradx/sgradz/rgradx/rgradz.\n"));
    } else if (cig && cig_is_angle) {
        INFO(("mode=angle: no spatial gradient files are provided; compute dT/dx and dT/dz by finite differences of traveltime tables.\n"));
    }

    if (!se_have_par("tau"))
        tau = 0.0f;
    else
        tau = se_get_par_float("tau"); /* static time-shift, in second */
    INFO(("tau = %f\n", tau));

    if (!se_have_par("aperture"))
        aper = 90.0f;
    else
        aper = se_get_par_float("aperture"); /* migration aperture, in degree */
    if (!se_have_par("antialias"))
        aal = 1.0f;
    else
        aal = se_get_par_float("antialias"); /* antialiasing */

    /*
     * Angle-domain cleanup controls.
     *
     * angle_interp=1:
     *   linearly distribute one migration contribution to two neighboring
     *   angle bins instead of nearest-neighbor binning.
     *
     * angle_taper_start / angle_taper_end:
     *   optional cosine taper in degrees. A negative start disables tapering.
     */
    if (se_have_par("angle_interp"))
        angle_interp = se_get_par_int("angle_interp");
    if (se_have_par("angle_taper_start"))
        angle_taper_start = se_get_par_float("angle_taper_start");
    if (se_have_par("angle_taper_end"))
        angle_taper_end = se_get_par_float("angle_taper_end");
    if (se_have_par("angle_grad_radius"))
        angle_grad_radius = se_get_par_int("angle_grad_radius");
    if (angle_grad_radius < 1)
        angle_grad_radius = 1;

    if (se_have_par("aperture_trace_taper"))
        aperture_trace_taper = se_get_par_int("aperture_trace_taper");
    if (aperture_trace_taper < 0)
        aperture_trace_taper = 0;

    if (cig && cig_is_angle) {
        INFO(("angle_interp = %d, angle_taper_start = %g, angle_taper_end = %g, angle_grad_radius = %d, aperture_trace_taper = %d\n",
              angle_interp, angle_taper_start, angle_taper_end,
              angle_grad_radius, aperture_trace_taper));
    }

    if (cig) {
        parse_output_axis(cig_is_angle, nh, h0, dh, &cig0, &dcig, &ng);
    } else {
        ng = 1;
        cig0 = 0.0f;
        dcig = 1.0f;
    }

    INFO(("cig = %d, mode = %s, n3 = %d, o3 = %g, d3 = %g\n", cig, cig_is_angle ? "angle" : "offset", ng, cig0, dcig));

    mig->headers->ndim = 3;

    if (adj) {
        mig->headers->n[0] = nz;
        mig->headers->n[1] = nx;
        mig->headers->o[0] = z0;
        mig->headers->d[0] = dz;
        mig->headers->o[1] = x0;
        mig->headers->d[1] = dx;
        sep_set_header(mig, "label1", "Depth");
        sep_set_header(mig, "label2", "Lateral");
        unit = sep_get_hdr(dat, "unit1", NULL);

        if (NULL != unit)
            sep_set_header(mig, "unit1", unit);
        if (cig) {
            mig->headers->n[2] = ng;
            mig->headers->o[2] = cig0;
            mig->headers->d[2] = dcig;
            if (cig_is_angle) {
                sep_set_header(mig, "label3", "Half-opening angle");
                sep_set_header(mig, "unit3", "degree");
            } else {
                sep_set_header(mig, "label3", cmp ? "Offset" : "Receiver");
                if (NULL != unit)
                    sep_set_header(mig, "unit3", unit);
            }
        } else {
            mig->headers->n[2] = 1;
            mig->headers->o[2] = 0.0f;
            mig->headers->d[2] = 1.0f;
        }
    } else {
        dat->headers->n[0] = nt;
        dat->headers->n[1] = nh;
        dat->headers->n[2] = ns;

        dat->headers->o[0] = t0;
        dat->headers->d[0] = dt;
        dat->headers->o[1] = h0;
        dat->headers->d[1] = dh;
        dat->headers->o[2] = s0;
        dat->headers->d[2] = ds;

        sep_set_header(dat, "label1", "Time");
        sep_set_header(dat, "unit1", "s");

        if (cmp)
            sep_set_header(dat, "label2", "Offset");
        else
            sep_set_header(dat, "label2", "Receiver");
        sep_set_header(dat, "label3", "Shot");
    }

    /* allocate temporary memory */
    out = alloc2float(nzx, ng);
    trace = alloc1float(nt);

    stable = alloc1float(nzx);
    stablex = alloc1float(nzx);
    rtable = alloc1float(nzx);
    rtablex = alloc1float(nzx);
    stbl0 = alloc1float(nzx);
    stbl1 = alloc1float(nzx);
    stblx0 = alloc1float(nzx);
    stblx1 = alloc1float(nzx);
    rtbl0 = alloc1float(nzx);
    rtbl1 = alloc1float(nzx);
    rtblx0 = alloc1float(nzx);
    rtblx1 = alloc1float(nzx);
    if (use_input_angle_grad) {
        grad_slice0 = alloc1float(nzx);
        grad_slice1 = alloc1float(nzx);
    }

    /* type of interpolation, default Hermite */
    if (!se_have_par("type"))
        type = "hermit";
    else
        type = se_get_par_str("type");

    if (type[0] != 'l' && type[0] != 'p' && type[0] != 'h') {
        ERROR(("Unknown interpolation type: %s. Supported types are linear (l), partial (p), and hermite (h)", type));
    }

    /* initialize interpolation */
    tinterp_init(nzx, sdy, rdy);

    /* initialize summation */
    kirmig_init(nt, dt, t0);

    if (adj) {
        memset(out[0], 0, nzx * ng * sizeof(float));
    } else {
        sep_read_fromsep(mig, nzx * ng, 0, &out[0], NULL);
    }

    auto t_start = std::chrono::steady_clock::now();

    /* determine number of threads without a data race */
#ifdef _OPENMP
    nthr = omp_get_max_threads();
#else
    nthr = 1;
#endif

    INFO((">>Using %d threads<<\n", nthr));

    if (!adj) {
        traces = alloc2float(nt, nthr);
    }

    for (is = 0; is < ns; is++) { /* shot */
        s = s0 + is * ds;
        INFO(("shot %d of %d;", is + 1, ns));

        /* cubic Hermite spline interpolation of source traveltime */
        ist = (int)((s - sy0) / sdy);
        if (ist <= 0) {
            load_table_slice(stim, 0, stbl0, nzx, &stbl0_idx, "stable");
            load_table_slice(sder, 0, stblx0, nzx, &stblx0_idx, "sderiv");
            copy_n(stable, stbl0, nzx);
            copy_n(stablex, stblx0, nzx);
        } else if (ist >= sny - 1) {
            load_table_slice(stim, sny - 1, stbl0, nzx, &stbl0_idx, "stable");
            load_table_slice(sder, sny - 1, stblx0, nzx, &stblx0_idx, "sderiv");
            copy_n(stable, stbl0, nzx);
            copy_n(stablex, stblx0, nzx);
        } else {
            load_table_slice(stim, ist, stbl0, nzx, &stbl0_idx, "stable");
            load_table_slice(stim, ist + 1, stbl1, nzx, &stbl1_idx, "stable");
            load_table_slice(sder, ist, stblx0, nzx, &stblx0_idx, "sderiv");
            load_table_slice(sder, ist + 1, stblx1, nzx, &stblx1_idx, "sderiv");
            switch (type[0]) {
            case 'l': /* linear */
                tinterp_linear(true, stable, s - ist * sdy - sy0, stbl0, stbl1);
                dinterp_linear(true, stablex, s - ist * sdy - sy0, stbl0, stbl1);
                break;

            case 'p': /* partial */
                tinterp_partial(true, stable, s - ist * sdy - sy0, nz, nx, dx, stbl0, stbl1);
                dinterp_partial(true, stablex, s - ist * sdy - sy0, nz, nx, dx, stbl0, stbl1);
                break;

            case 'h': /* Hermite */
                tinterp_hermite(true, stable, s - ist * sdy - sy0, stbl0, stbl1, stblx0, stblx1);
                dinterp_hermite(true, stablex, s - ist * sdy - sy0, stbl0, stbl1, stblx0, stblx1);
                break;
            }
        }

        if (use_input_angle_grad) {
            interp_volume_linear_stream(sgradx, sgradx_fp, sny, nzx, s, sy0, sdy,
                                        grad_slice0, grad_slice1, "sgradx");
            interp_volume_linear_stream(sgradz, sgradz_fp, sny, nzx, s, sy0, sdy,
                                        grad_slice0, grad_slice1, "sgradz");
        }

        for (ih = 0; ih < nh; ih++) { /* offset or receiver */
            h = h0 + ih * dh;

            /* cubic Hermite spline interpolation of receiver traveltime */
            iht = cmp ? (int)((s + h - ry0) / rdy) : (int)((h - ry0) / rdy);

            float trace_ap_weight = 1.0f;

            if (adj) {
                /* read trace first to keep the input I/O position correct */
                se_fsio_read_float(dat->data->io, trace, nt);

                if (direct_mute && is == direct_mute_qc_shot)
                    se_fsio_write_float(mute_before->data->io, trace, nt);

                if (direct_mute) {
                    const float receiver_x = cmp ? (s + h) : h;
                    apply_direct_wave_mute(
                        trace, nt, t0, dt, s, receiver_x,
                        direct_mute_slope, direct_mute_intercept,
                        direct_mute_taper, direct_mute_invert);
                }

                if (direct_mute && is == direct_mute_qc_shot)
                    se_fsio_write_float(mute_after->data->io, trace, nt);

                trace_ap_weight = 1.0f;
                if (aperture_trace != -1) {
                    const float ad = fabsf((float)(iht - ist));

                    if (ad >= (float)aperture_trace)
                        continue;

                    if (aperture_trace_taper > 0) {
                        const float taper_begin =
                            MAX(0.0f, (float)(aperture_trace - aperture_trace_taper));

                        if (ad > taper_begin) {
                            const float denom =
                                (float)aperture_trace - taper_begin;

                            if (denom > 0.0f) {
                                const float x = (ad - taper_begin) / denom;
                                trace_ap_weight =
                                    0.5f * (1.0f + cosf((float)M_PI * x));
                            }
                        }
                    }
                }

                if (cig && !cig_is_angle) {
                    ig_offset = value_to_nearest_bin(h, cig0, dcig, ng);
                    if (ig_offset < 0) continue;
                } else {
                    ig_offset = 0;
                }

                doubint(nt, trace);

                if (isdiff == 1) {
                    diff2(trace, nt, dt);
                }
            } else {
                if (aperture_trace != -1) {
                    const float ad = fabsf((float)(iht - ist));

                    if (ad >= (float)aperture_trace) {
                        for (i = 0; i < nt; i++) trace[i] = 0.0f;
                        se_fsio_write_float(dat->data->io, trace, nt);
                        continue;
                    }

                    if (aperture_trace_taper > 0) {
                        const float taper_begin =
                            MAX(0.0f, (float)(aperture_trace - aperture_trace_taper));

                        if (ad > taper_begin) {
                            const float denom =
                                (float)aperture_trace - taper_begin;

                            if (denom > 0.0f) {
                                const float x = (ad - taper_begin) / denom;
                                trace_ap_weight =
                                    0.5f * (1.0f + cosf((float)M_PI * x));
                            }
                        }
                    }
                }

                if (cig && !cig_is_angle) {
                    ig_offset = value_to_nearest_bin(h, cig0, dcig, ng);
                    if (ig_offset < 0) {
                        for (i = 0; i < nt; i++) trace[i] = 0.0f;
                        se_fsio_write_float(dat->data->io, trace, nt);
                        continue;
                    }
                } else {
                    ig_offset = 0;
                }
                for (ithr = 0; ithr < nthr; ithr++) {
                    for (i = 0; i < nt; i++) {
                        traces[ithr][i] = 0.0f;
                    }
                }
            }

            if (iht <= 0) {
                load_table_slice(rtim, 0, rtbl0, nzx, &rtbl0_idx, "rtable");
                load_table_slice(rder, 0, rtblx0, nzx, &rtblx0_idx, "rderiv");
                copy_n(rtable, rtbl0, nzx);
                copy_n(rtablex, rtblx0, nzx);
            } else if (iht >= rny - 1) {
                load_table_slice(rtim, rny - 1, rtbl0, nzx, &rtbl0_idx, "rtable");
                load_table_slice(rder, rny - 1, rtblx0, nzx, &rtblx0_idx, "rderiv");
                copy_n(rtable, rtbl0, nzx);
                copy_n(rtablex, rtblx0, nzx);
            } else {
                float rcoord = cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0;
                load_table_slice(rtim, iht, rtbl0, nzx, &rtbl0_idx, "rtable");
                load_table_slice(rtim, iht + 1, rtbl1, nzx, &rtbl1_idx, "rtable");
                load_table_slice(rder, iht, rtblx0, nzx, &rtblx0_idx, "rderiv");
                load_table_slice(rder, iht + 1, rtblx1, nzx, &rtblx1_idx, "rderiv");
                switch (type[0]) {
                case 'l': /* linear */
                    tinterp_linear(false, rtable, rcoord, rtbl0, rtbl1);
                    dinterp_linear(false, rtablex, rcoord, rtbl0, rtbl1);
                    break;

                case 'p': /* partial */
                    tinterp_partial(false, rtable, rcoord, nz, nx, dx, rtbl0, rtbl1);
                    dinterp_partial(false, rtablex, rcoord, nz, nx, dx, rtbl0, rtbl1);
                    break;

                case 'h': /* Hermite */
                    tinterp_hermite(false, rtable, rcoord, rtbl0, rtbl1, rtblx0, rtblx1);
                    dinterp_hermite(false, rtablex, rcoord, rtbl0, rtbl1, rtblx0, rtblx1);
                    break;
                }
            }

            if (use_input_angle_grad) {
                float rpos = cmp ? (s + h) : h;
                interp_volume_linear_stream(rgradx, rgradx_fp, rny, nzx, rpos, ry0, rdy,
                                            grad_slice0, grad_slice1, "rgradx");
                interp_volume_linear_stream(rgradz, rgradz_fp, rny, nzx, rpos, ry0, rdy,
                                            grad_slice0, grad_slice1, "rgradz");
            }

#ifdef _OPENMP
#pragma omp parallel for private(iz, ix, t1, t2, ti, tx, ithr, ig, angle)
#endif
            for (i = 0; i < nzx; i++) {
#ifdef _OPENMP
                ithr = omp_get_thread_num();
#else
                ithr = 0;
#endif
                iz = i % nz;
                ix = (i - iz) / nz;

                /* aperture, cone angle */
                if (cmp) {
                    if (h >= 0.0f) {
                        if (atanf((s - x0 - ix * dx) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - s - h) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                    } else {
                        if (atanf((s + h - x0 - ix * dx) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - s) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                    }
                } else {
                    if (h - s >= 0.0f) {
                        if (atanf((s - x0 - ix * dx) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - h) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                    } else {
                        if (atanf((h - x0 - ix * dx) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - s) / (iz * dz)) * 180.0f / (float)M_PI > aper)
                            continue;
                    }
                }

                t1 = stable[i];
                t2 = rtable[i];
                ti = t1 + t2 + tau;

                tx = MAX(fabsf(stablex[i] * ds), fabsf(rtablex[i] * dh));

                if (cig && cig_is_angle) {
                    if (use_input_angle_grad) {
                        angle = half_opening_angle_deg_from_grad(
                            sgradx[i], sgradz[i], rgradx[i], rgradz[i]);
                    } else {
                        angle = half_opening_angle_deg_from_tables(
                            stable, rtable, nz, nx, iz, ix, dx, dz,
                            angle_grad_radius);
                    }

                    const float aw =
                        angle_taper_weight(
                            angle, angle_taper_start, angle_taper_end) *
                        trace_ap_weight;

                    if (aw <= 0.0f) continue;

                    if (angle_interp) {
                        int ig0_local = -1, ig1_local = -1;
                        float w0_local = 0.0f, w1_local = 0.0f;

                        if (!value_to_linear_bins(
                                angle, cig0, dcig, ng,
                                &ig0_local, &ig1_local,
                                &w0_local, &w1_local)) {
                            continue;
                        }

                        if (adj) {
                            /*
                             * First evaluate the Kirchhoff sample once, then
                             * distribute it between the two neighboring angle
                             * bins. Each OpenMP iteration owns one spatial
                             * sample i, so these two updates are race-free.
                             */
                            float picked = 0.0f;
                            kirmig_pick(1, ti, tx * aal, &picked, trace);
                            picked *= aw;

                            out[ig0_local][i] += w0_local * picked;
                            if (ig1_local != ig0_local && w1_local > 0.0f)
                                out[ig1_local][i] += w1_local * picked;
                        } else {
                            /*
                             * Forward/adjoint consistency: the same two-bin
                             * interpolation and taper are applied when modeling
                             * data from an angle gather.
                             */
                            float v0 = aw * w0_local * out[ig0_local][i];
                            kirmig_pick(0, ti, tx * aal, &v0, traces[ithr]);

                            if (ig1_local != ig0_local && w1_local > 0.0f) {
                                float v1 = aw * w1_local * out[ig1_local][i];
                                kirmig_pick(0, ti, tx * aal, &v1, traces[ithr]);
                            }
                        }
                    } else {
                        ig = value_to_nearest_bin(angle, cig0, dcig, ng);
                        if (ig < 0) continue;

                        if (adj) {
                            float picked = 0.0f;
                            kirmig_pick(1, ti, tx * aal, &picked, trace);
                            out[ig][i] += aw * picked;
                        } else {
                            float v = aw * out[ig][i];
                            kirmig_pick(0, ti, tx * aal, &v, traces[ithr]);
                        }
                    }
                } else {
                    if (cig)
                        ig = ig_offset;
                    else
                        ig = 0;

                    if (adj) {
                        float picked = 0.0f;
                        kirmig_pick(1, ti, tx * aal, &picked, trace);
                        out[ig][i] += trace_ap_weight * picked;
                    } else {
                        float v = trace_ap_weight * out[ig][i];
                        kirmig_pick(0, ti, tx * aal, &v, traces[ithr]);
                    }
                }
            }

            if (!adj) {
                for (i = 0; i < nt; i++) {
                    trace[i] = 0.0f;
                }
                for (ithr = 0; ithr < nthr; ithr++) {
                    for (i = 0; i < nt; i++) {
                        trace[i] += traces[ithr][i];
                    }
                }
                doubint(nt, trace);
                se_fsio_write_float(dat->data->io, trace, nt);
            }
        } /* ih */
    }

    INFO(("FINISH."));
    auto t_end = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration<double>(t_end - t_start).count();
    INFO(("Done. Elapsed time: %.3f s.", elapsed_seconds));

    if (adj)
        se_fsio_write_float(mig->data->io, out[0], nzx * ng);

    se_par_destroy();
    sep_close(dat);
    sep_close(mig);
    sep_close(stim);
    sep_close(sder);
    sep_close(rtim);
    sep_close(rder);
    if (sgradx_fp) sep_close(sgradx_fp);
    if (sgradz_fp) sep_close(sgradz_fp);
    if (rgradx_fp) sep_close(rgradx_fp);
    if (rgradz_fp) sep_close(rgradz_fp);
    if (mute_before) sep_close(mute_before);
    if (mute_after) sep_close(mute_after);

    if (traces) free2float(traces);
    if (out) free2float(out);
    if (trace) free1float(trace);
    if (stable) free1float(stable);
    if (stablex) free1float(stablex);
    if (rtable) free1float(rtable);
    if (rtablex) free1float(rtablex);
    if (stbl0) free1float(stbl0);
    if (stbl1) free1float(stbl1);
    if (stblx0) free1float(stblx0);
    if (stblx1) free1float(stblx1);
    if (rtbl0) free1float(rtbl0);
    if (rtbl1) free1float(rtbl1);
    if (rtblx0) free1float(rtblx0);
    if (rtblx1) free1float(rtblx1);
    if (grad_slice0) free1float(grad_slice0);
    if (grad_slice1) free1float(grad_slice1);

    if (sgradx) free1float(sgradx);
    if (sgradz) free1float(sgradz);
    if (rgradx) free1float(rgradx);
    if (rgradz) free1float(rgradz);

    return 0;
}
