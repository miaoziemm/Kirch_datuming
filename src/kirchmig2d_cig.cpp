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

static float finite_diff_x(const float *tab, int nz, int nx, int iz, int ix, float dx)
{
    if (nx <= 1 || dx == 0.0f) return 0.0f;
    if (ix <= 0) {
        return (tab[idx2d(iz, 1, nz)] - tab[idx2d(iz, 0, nz)]) / dx;
    }
    if (ix >= nx - 1) {
        return (tab[idx2d(iz, nx - 1, nz)] - tab[idx2d(iz, nx - 2, nz)]) / dx;
    }
    return (tab[idx2d(iz, ix + 1, nz)] - tab[idx2d(iz, ix - 1, nz)]) / (2.0f * dx);
}

static float finite_diff_z(const float *tab, int nz, int nx, int iz, int ix, float dz)
{
    (void)nx;
    if (nz <= 1 || dz == 0.0f) return 0.0f;
    if (iz <= 0) {
        return (tab[idx2d(1, ix, nz)] - tab[idx2d(0, ix, nz)]) / dz;
    }
    if (iz >= nz - 1) {
        return (tab[idx2d(nz - 1, ix, nz)] - tab[idx2d(nz - 2, ix, nz)]) / dz;
    }
    return (tab[idx2d(iz + 1, ix, nz)] - tab[idx2d(iz - 1, ix, nz)]) / (2.0f * dz);
}

static float half_opening_angle_deg_from_grad(float sx, float sz, float rx, float rz)
{
    float ns = sqrtf(sx * sx + sz * sz);
    float nr = sqrtf(rx * rx + rz * rz);
    if (ns <= 1.0e-12f || nr <= 1.0e-12f) return NAN;

    /*
       Here sx,sz and rx,rz are the spatial gradients of the source and
       receiver traveltimes at the image point. In 2D, the ADCIG angle is
       taken as the half-opening angle:

           alpha = 0.5 * acos( grad(Ts) dot grad(Tr) / |grad(Ts)||grad(Tr)| ).
    */
    float c = (sx * rx + sz * rz) / (ns * nr);
    c = clampf_local(c, -1.0f, 1.0f);
    return 0.5f * acosf(c) * 180.0f / (float)M_PI;
}

static float half_opening_angle_deg_from_tables(const float *stable,
                                                const float *rtable,
                                                int nz, int nx,
                                                int iz, int ix,
                                                float dx, float dz)
{
    float sx = finite_diff_x(stable, nz, nx, iz, ix, dx);
    float sz = finite_diff_z(stable, nz, nx, iz, ix, dz);
    float rx = finite_diff_x(rtable, nz, nx, iz, ix, dx);
    float rz = finite_diff_z(rtable, nz, nx, iz, ix, dz);
    return half_opening_angle_deg_from_grad(sx, sz, rx, rz);
}

static void copy_n(float *dst, const float *src, off_t n)
{
    for (off_t i = 0; i < n; i++) dst[i] = src[i];
}

static void interp_volume_linear(float *out, float **tbl, int npos, off_t nxy,
                                 float coord, float o, float d)
{
    if (npos <= 1 || d == 0.0f) {
        copy_n(out, tbl[0], nxy);
        return;
    }

    float f = (coord - o) / d;
    int i0 = (int)floorf(f);

    if (i0 <= 0) {
        copy_n(out, tbl[0], nxy);
        return;
    }
    if (i0 >= npos - 1) {
        copy_n(out, tbl[npos - 1], nxy);
        return;
    }

    float w = f - (float)i0;
    const float *a = tbl[i0];
    const float *b = tbl[i0 + 1];
    for (off_t i = 0; i < nxy; i++) {
        out[i] = (1.0f - w) * a[i] + w * b[i];
    }
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

static void read_grad_volume(const char *name, char *file, int nz, int nx, int n3,
                             off_t nzx, float ***tbl_out)
{
    sep_t *fp = sep_open(file, SEP_READ, 0);
    if (fp->headers->ndim < 3) {
        ERROR(("Need 3D data volume for %s= (z,x,source_or_receiver)", name));
    }
    if (fp->headers->n[0] != nz || fp->headers->n[1] != nx || fp->headers->n[2] != n3) {
        ERROR(("Dimension mismatch in %s=. Expected n1=%d n2=%d n3=%d", name, nz, nx, n3));
    }
    *tbl_out = alloc2float(nzx, n3);
    se_fsio_read_float(fp->data->io, (*tbl_out)[0], nzx * n3);
    sep_close(fp);
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
    float **stbl = NULL, **rtbl = NULL, *stable = NULL, *rtable = NULL;
    float **stblx = NULL, **rtblx = NULL, *stablex = NULL, *rtablex = NULL;
    float **sgradx_tbl = NULL, **sgradz_tbl = NULL, **rgradx_tbl = NULL, **rgradz_tbl = NULL;
    float *sgradx = NULL, *sgradz = NULL, *rgradx = NULL, *rgradz = NULL;
    float ds = 0.0f, s0 = 0.0f, x0 = 0.0f, sy0 = 0.0f, sdy = 0.0f, ry0 = 0.0f, rdy = 0.0f;
    float s = 0.0f, h = 0.0f, h0 = 0.0f, dh = 0.0f, dx = 0.0f;
    float ti = 0.0f, t0 = 0.0f, t1 = 0.0f, t2 = 0.0f, dt = 0.0f, z0 = 0.0f, dz = 0.0f, tau = 0.0f;
    float aal = 0.0f, tx = 0.0f, aper = 0.0f;
    float cig0 = 0.0f, dcig = 1.0f, angle = 0.0f;
    sep_t *dat = NULL, *mig = NULL, *stim = NULL, *sder = NULL, *rtim = NULL, *rder = NULL;
    char *dat_f = NULL, *mig_f = NULL, *stim_f = NULL, *sder_f = NULL, *rtim_f = NULL, *rder_f = NULL;
    char *sgradx_f = NULL, *sgradz_f = NULL, *rgradx_f = NULL, *rgradz_f = NULL;

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

    stbl = alloc2float(nzx, sny);
    se_fsio_read_float(stim->data->io, stbl[0], nzx * sny);
    sep_close(stim);

    stblx = alloc2float(nzx, sny);
    se_fsio_read_float(sder->data->io, stblx[0], nzx * sny);
    sep_close(sder);

    rtim = sep_open(rtim_f, SEP_READ, 0);
    rder = sep_open(rder_f, SEP_READ, 0);

    if (rtim->headers->ndim < 3)
        ERROR(("Need 3D data volume for rtable= (z,x,r)"));
    if (rder->headers->ndim < 3)
        ERROR(("Need 3D data volume for rderiv= (z,x,r)"));

    rny = rtim->headers->n[2];
    ry0 = (float)rtim->headers->o[2];
    rdy = (float)rtim->headers->d[2];

    rtbl = alloc2float(nzx, rny);
    se_fsio_read_float(rtim->data->io, rtbl[0], nzx * rny);
    sep_close(rtim);

    rtblx = alloc2float(nzx, rny);
    se_fsio_read_float(rder->data->io, rtblx[0], nzx * rny);
    sep_close(rder);

    if (use_input_angle_grad) {
        read_grad_volume("sgradx", sgradx_f, nz, nx, sny, nzx, &sgradx_tbl);
        read_grad_volume("sgradz", sgradz_f, nz, nx, sny, nzx, &sgradz_tbl);
        read_grad_volume("rgradx", rgradx_f, nz, nx, rny, nzx, &rgradx_tbl);
        read_grad_volume("rgradz", rgradz_f, nz, nx, rny, nzx, &rgradz_tbl);
        sgradx = alloc1float(nzx);
        sgradz = alloc1float(nzx);
        rgradx = alloc1float(nzx);
        rgradz = alloc1float(nzx);
        INFO(("mode=angle: use input spatial traveltime gradients sgradx/sgradz/rgradx/rgradz.\n"));
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
            for (i = 0; i < nzx; i++) {
                stable[i] = stbl[0][i];
                stablex[i] = stblx[0][i];
            }
        } else if (ist >= sny - 1) {
            for (i = 0; i < nzx; i++) {
                stable[i] = stbl[sny - 1][i];
                stablex[i] = stblx[sny - 1][i];
            }
        } else {
            switch (type[0]) {
            case 'l': /* linear */
                tinterp_linear(true, stable, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1]);
                dinterp_linear(true, stablex, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1]);
                break;

            case 'p': /* partial */
                tinterp_partial(true, stable, s - ist * sdy - sy0, nz, nx, dx, stbl[ist], stbl[ist + 1]);
                dinterp_partial(true, stablex, s - ist * sdy - sy0, nz, nx, dx, stbl[ist], stbl[ist + 1]);
                break;

            case 'h': /* Hermite */
                tinterp_hermite(true, stable, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1], stblx[ist], stblx[ist + 1]);
                dinterp_hermite(true, stablex, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1], stblx[ist], stblx[ist + 1]);
                break;
            }
        }

        if (use_input_angle_grad) {
            interp_volume_linear(sgradx, sgradx_tbl, sny, nzx, s, sy0, sdy);
            interp_volume_linear(sgradz, sgradz_tbl, sny, nzx, s, sy0, sdy);
        }

        for (ih = 0; ih < nh; ih++) { /* offset or receiver */
            h = h0 + ih * dh;

            /* cubic Hermite spline interpolation of receiver traveltime */
            iht = cmp ? (int)((s + h - ry0) / rdy) : (int)((h - ry0) / rdy);

            if (adj) {
                /* read trace first to keep the input I/O position correct */
                se_fsio_read_float(dat->data->io, trace, nt);

                if (aperture_trace != -1 && fabsf((float)(iht - ist)) >= (float)aperture_trace) continue;

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
                for (i = 0; i < nzx; i++) {
                    rtable[i] = rtbl[0][i];
                    rtablex[i] = rtblx[0][i];
                }
            } else if (iht >= rny - 1) {
                for (i = 0; i < nzx; i++) {
                    rtable[i] = rtbl[rny - 1][i];
                    rtablex[i] = rtblx[rny - 1][i];
                }
            } else {
                float rcoord = cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0;
                switch (type[0]) {
                case 'l': /* linear */
                    tinterp_linear(false, rtable, rcoord, rtbl[iht], rtbl[iht + 1]);
                    dinterp_linear(false, rtablex, rcoord, rtbl[iht], rtbl[iht + 1]);
                    break;

                case 'p': /* partial */
                    tinterp_partial(false, rtable, rcoord, nz, nx, dx, rtbl[iht], rtbl[iht + 1]);
                    dinterp_partial(false, rtablex, rcoord, nz, nx, dx, rtbl[iht], rtbl[iht + 1]);
                    break;

                case 'h': /* Hermite */
                    tinterp_hermite(false, rtable, rcoord, rtbl[iht], rtbl[iht + 1], rtblx[iht], rtblx[iht + 1]);
                    dinterp_hermite(false, rtablex, rcoord, rtbl[iht], rtbl[iht + 1], rtblx[iht], rtblx[iht + 1]);
                    break;
                }
            }

            if (use_input_angle_grad) {
                float rpos = cmp ? (s + h) : h;
                interp_volume_linear(rgradx, rgradx_tbl, rny, nzx, rpos, ry0, rdy);
                interp_volume_linear(rgradz, rgradz_tbl, rny, nzx, rpos, ry0, rdy);
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

                if (cig) {
                    if (cig_is_angle) {
                        if (use_input_angle_grad) {
                            angle = half_opening_angle_deg_from_grad(sgradx[i], sgradz[i], rgradx[i], rgradz[i]);
                        } else {
                            angle = half_opening_angle_deg_from_tables(stable, rtable, nz, nx, iz, ix, dx, dz);
                        }
                        ig = value_to_nearest_bin(angle, cig0, dcig, ng);
                        if (ig < 0) continue;
                    } else {
                        ig = ig_offset;
                    }
                } else {
                    ig = 0;
                }

                t1 = stable[i];
                t2 = rtable[i];
                ti = t1 + t2 + tau;

                tx = MAX(fabsf(stablex[i] * ds), fabsf(rtablex[i] * dh));

                kirmig_pick(adj, ti, tx * aal, out[ig] + i, adj ? trace : traces[ithr]);
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

    if (traces) free2float(traces);
    if (out) free2float(out);
    if (trace) free1float(trace);
    if (stable) free1float(stable);
    if (stablex) free1float(stablex);
    if (rtable) free1float(rtable);
    if (rtablex) free1float(rtablex);
    if (stbl) free2float(stbl);
    if (stblx) free2float(stblx);
    if (rtbl) free2float(rtbl);
    if (rtblx) free2float(rtblx);

    if (sgradx) free1float(sgradx);
    if (sgradz) free1float(sgradz);
    if (rgradx) free1float(rgradx);
    if (rgradz) free1float(rgradz);
    if (sgradx_tbl) free2float(sgradx_tbl);
    if (sgradz_tbl) free2float(sgradz_tbl);
    if (rgradx_tbl) free2float(rgradx_tbl);
    if (rgradz_tbl) free2float(rgradz_tbl);

    return 0;
}
