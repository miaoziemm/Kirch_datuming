#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>
#include "./help/eikods2d_help.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <limits>
#include <vector>
#include <unistd.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

/*
 * Branch-aware spatial traveltime gradient along one axis.
 *
 * A centered finite difference, (T[i+1]-T[i-1])/(2 d), straddles the
 * first-arrival crease and averages two different ray branches, which makes
 * the scattering angle jump.  This routine instead back-traces the ray: it
 * uses only the axis neighbor with the SMALLER traveltime (an upwind,
 * one-sided stencil; a 3-point upwind stencil when it stays on the same
 * branch).  At a crease it therefore follows one ray leg, exactly like the
 * ray-vector angle used by Madagascar cram.
 *
 * Flat indexing is axis0-fast (z): ss0=1, ss1=n1, ss2=n1*n2.
 */
static float upwind_grad_axis(const float *t, int i0, int axis,
                              const int nn[3], const int ss[3],
                              const float dd[3])
{
    int ii[3], i1, i2;

    ii[0] = i0 % nn[0];
    ii[1] = (i0 / ss[1]) % nn[1];
    ii[2] = i0 / ss[2];

    const int a = i0 - ss[axis];
    const int b = i0 + ss[axis];

    if (ii[axis] == 0) {
        i1 = b;                                   /* boundary: only + side */
    } else if (ii[axis] != nn[axis] - 1 && t[b] < t[a]) {
        i1 = b;                                   /* upwind is the + side */
    } else {
        i1 = a;                                   /* upwind is the - side */
    }

    if (!(t[i1] < t[i0])) return 0.0f;            /* no upwind along axis */

    float g;
    if (i1 == b) {
        i2 = i1 + ss[axis];
        if (ii[axis] < nn[axis] - 2 && t[i2] < t[i1])
            g = (-t[i2] + 4.0f * t[i1] - 3.0f * t[i0]) / (2.0f * dd[axis]);
        else
            g = (t[i1] - t[i0]) / dd[axis];
    } else {
        i2 = i1 - ss[axis];
        if (ii[axis] > 1 && t[i2] < t[i1])
            g = (3.0f * t[i0] - 4.0f * t[i1] + t[i2]) / (2.0f * dd[axis]);
        else
            g = (t[i0] - t[i1]) / dd[axis];
    }
    return g;
}

static void write_float_shot_at(sep_t *sep, const float *buf, size_t nsample, int ishot, const char *label)
{
    if (sep == NULL || sep->data == NULL || sep->data->io == NULL)
        ERROR(("Invalid SEP output handle for %s", label));

    se_fsio *io = sep->data->io;
    if (io->io_data_type != FIO_DATA_TYPE_NATIVE_FLOAT)
        ERROR(("Parallel random-access output for %s requires native_float data format", label));

    if (nsample > std::numeric_limits<size_t>::max() / sizeof(float))
        ERROR(("Byte-count overflow while writing %s", label));
    const size_t nbytes = nsample * sizeof(float);
    if ((size_t)ishot > std::numeric_limits<size_t>::max() / nbytes)
        ERROR(("Offset overflow while writing %s", label));

    const off_t offset = (off_t)((size_t)ishot * nbytes);
    const char *bytes = reinterpret_cast<const char *>(buf);
    size_t done = 0;
    while (done < nbytes)
    {
        ssize_t written = pwrite(io->fd, bytes + done, nbytes - done, offset + (off_t)done);
        if (written < 0)
        {
            if (errno == EINTR)
                continue;
            ERROR(("pwrite failed for %s shot %d file %s: errno=%d - %s",
                   label, ishot + 1, io->file_name, errno, strerror(errno)));
        }
        if (written == 0)
            ERROR(("pwrite wrote zero bytes for %s shot %d file %s", label, ishot + 1, io->file_name));
        done += (size_t)written;
    }
}

static void write_float_all(sep_t *sep, const float *buf, size_t nsample, int nshot, const char *label)
{
    if (sep == NULL || sep->data == NULL || sep->data->io == NULL)
        ERROR(("Invalid SEP output handle for %s", label));

    if (nsample > std::numeric_limits<size_t>::max() / (sizeof(float) * (size_t)nshot))
        ERROR(("Byte-count overflow while writing %s", label));

    const size_t nbytes = nsample * (size_t)nshot * sizeof(float);
    const char *bytes = reinterpret_cast<const char *>(buf);
    size_t done = 0;
    while (done < nbytes)
    {
        ssize_t written = pwrite(sep->data->io->fd, bytes + done, nbytes - done, (off_t)done);
        if (written < 0)
        {
            if (errno == EINTR)
                continue;
            ERROR(("pwrite failed for %s file %s: errno=%d - %s",
                   label, sep->data->io->file_name, errno, strerror(errno)));
        }
        if (written == 0)
            ERROR(("pwrite wrote zero bytes for %s file %s", label, sep->data->io->file_name));
        done += (size_t)written;
    }
}

static void preallocate_float_output(sep_t *sep, size_t nsample, int nshot, const char *label)
{
    if (sep == NULL || sep->data == NULL || sep->data->io == NULL)
        return;

    if (nsample > std::numeric_limits<size_t>::max() / sizeof(float))
        ERROR(("Byte-count overflow while sizing %s", label));
    const size_t shot_bytes = nsample * sizeof(float);

    if ((size_t)nshot > std::numeric_limits<size_t>::max() / shot_bytes)
        ERROR(("Total byte-count overflow while sizing %s", label));
    const size_t total_bytes = shot_bytes * (size_t)nshot;

    if (ftruncate(sep->data->io->fd, (off_t)total_bytes) != 0)
        ERROR(("ftruncate failed for %s file %s: errno=%d - %s",
               label, sep->data->io->file_name, errno, strerror(errno)));
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        eikods2d_help::print_help();
        return 1;
    }

    int b1 = 1, b2 = 1, b3 = 1, n1 = 0, n2 = 0, n3 = 1, i, nshot = 1, ndim = 3, is, order = 2, l = 1;
    size_t n123 = 0;
    float br1 = 0.0f, br2 = 0.0f, br3 = 0.0f, o1 = 0.0f, o2 = 0.0f, o3 = 0.0f, d1 = 0.0f, d2 = 0.0f, d3 = 0.0f;
    float **s, *v;
    const char *sfile = NULL, *in_f = NULL, *out_f = NULL, *gradz_f = NULL, *gradx_f = NULL, *tdl1_f = NULL, *tds1_f = NULL, *tdl2_f = NULL, *tds2_f = NULL;
    bool isvel = true, plane[3] = {false, false, false};
    int efmm = 0;
    sep_t *vel = NULL, *time = NULL, *shots = NULL, *gradz = NULL, *gradx = NULL, *tdl1 = NULL, *tds1 = NULL, *tdl2 = NULL, *tds2 = NULL;

    se_par_init(argc, argv);

    if (!se_have_par("in"))
        ERROR(("Need in="));
    else
        in_f = se_get_par_str("in");
    if (!se_have_par("out"))
        ERROR(("Need out="));
    else
        out_f = se_get_par_str("out");

    vel = sep_open(in_f, SEP_READ, 0);
    time = sep_open(out_f, SEP_WRITE, 0);

    time->headers->ndim = 3;


    printf("Reading velocity model from %s\n", in_f);
    printf("Writing traveltime to %s\n", out_f);
    printf("n1=%d, n2=%d, n3=%d\n", vel->headers->n[0], vel->headers->n[1], vel->headers->n[2]);

    n1 = vel->headers->n[0];
    n2 = vel->headers->n[1];
    n3 = vel->headers->n[2];
    d1 = (float)vel->headers->d[0];
    d2 = (float)vel->headers->d[1];
    d3 = (float)vel->headers->d[2];
    o1 = (float)vel->headers->o[0];
    o2 = (float)vel->headers->o[1];
    o3 = (float)vel->headers->o[2];

    // if(!sep_have_hdr_int(vel,"n1")) ERROR(("No n1= in input")); else n1 = sep_get_hdr_int(vel,"n1",0);
    // if(!sep_have_hdr_int(vel,"n2")) ERROR(("No n2= in input")); else n2 = sep_get_hdr_int(vel,"n2",0);
    // if(!sep_have_hdr_int(vel,"n3")) n3=1; else n3=sep_get_hdr_int(vel,"n3",1);

    // if(!sep_have_hdr_float(vel,"d1")) ERROR(("No d1= in input")); else d1=(float)sep_get_hdr_float(vel,"d1",0.0);
    // if(!sep_have_hdr_float(vel,"d2")) ERROR(("No d2= in input")); else d2=(float)sep_get_hdr_float(vel,"d2",0.0);
    // if(!sep_have_hdr_float(vel,"d3")) d3=d2; else d3=(float)sep_get_hdr_float(vel,"d3",d2);

    // if(!sep_have_hdr_float(vel,"o1")) o1=0.; else o1=(float)sep_get_hdr_float(vel,"o1",0.0);
    // if(!sep_have_hdr_float(vel,"o2")) o2=0.; else o2=(float)sep_get_hdr_float(vel,"o2",0.0);
    // if(!sep_have_hdr_float(vel,"o3")) o3=0.; else o3=(float)sep_get_hdr_float(vel,"o3",0.0);

    if (!se_have_par("vel"))
        isvel = true;
    else
        isvel = se_get_par_int("vel");
    /* if y, the input is velocity; n, slowness squared */

    if (!se_have_par("order"))
        order = 2;
    else
        order = se_get_par_int("order");
    /* [1,2] Accuracy order */

    if (!se_have_par("efmm"))
        efmm = 0;
    else
        efmm = se_get_par_int("efmm");
    /* if y, use fast efmming instead of fast marching */

    if (!se_have_par("br1"))
        br1 = d1;
    else
        br1 = se_get_par_float("br1");
    if (!se_have_par("br2"))
        br2 = d2;
    else
        br2 = se_get_par_float("br2");
    if (!se_have_par("br3"))
        br3 = d3;
    else
        br3 = se_get_par_float("br3");
    /* Constant-velocity box around the source (in physical dimensions) */

    if (!se_have_par("plane1"))
        plane[2] = false;
    else
        plane[2] = se_get_par_int("plane1");
    if (!se_have_par("plane2"))
        plane[1] = false;
    else
        plane[1] = se_get_par_int("plane2");
    if (!se_have_par("plane3"))
        plane[0] = false;
    else
        plane[0] = se_get_par_int("plane3");
    /* plane-wave source */

    if (!se_have_par("b1"))
        b1 = plane[2] ? n1 : (int)(br1 / d1 + 0.5);
    else
        b1 = se_get_par_int("b1");
    if (!se_have_par("b2"))
        b2 = plane[1] ? n2 : (int)(br2 / d2 + 0.5);
    else
        b2 = se_get_par_int("b2");
    if (!se_have_par("b3"))
        b3 = plane[0] ? n3 : (int)(br3 / d3 + 0.5);
    else
        b3 = se_get_par_int("b3");
    /* Constant-velocity box around the source (in samples) */

    if (b1 < 1)
        b1 = 1;
    if (b2 < 1)
        b2 = 1;
    if (b3 < 1)
        b3 = 1;

    if (se_have_par("shotfile"))
        sfile = se_get_par_str("shotfile");
    /* File with shot locations (n2=number of shots, n1=3) */

    float dshot=0.0;
    float oshot=0.0;

    if (NULL != sfile)
    {
        shots = sep_open(sfile, SEP_READ, 0);
        nshot = shots->headers->n[1];
        ndim = shots->headers->n[0];

        // if(!sep_have_hdr_int(shots,"n2"))
        //     ERROR(("No n2= in shotfile")); else nshot=sep_get_hdr_int(shots,"n2",0);
        // if(!sep_have_hdr_int(shots,"n1"))
        //     ERROR(("Need n1=3 in shotfile")); else ndim=sep_get_hdr_int(shots,"n1",0);
        if (ndim != 3)
            ERROR(("Need n1=3 in shotfile"));

        s = alloc2float(ndim, nshot);
        se_fsio_read_float(shots->data->io, s[0], nshot * ndim);
        sep_close(shots);
        oshot = s[0][1];
        dshot = -s[0][1] + s[1][1];

        // sep_copy_headers(time, vel);
        // sep_set_axis(time, 3, nshot, 0.0, 1.0, "Shot");
    }
    else
    {
        nshot = 1;
        ndim = 3;

        s = alloc2float(ndim, nshot);

        if (!se_have_par("zshot"))
            s[0][0] = 0.;
        else
            s[0][0] = se_get_par_float("zshot");
        /* Shot location (used if no shotfile) */
        if (!se_have_par("yshot"))
            s[0][1] = o2 + 0.5 * (n2 - 1) * d2;
        else
            s[0][1] = se_get_par_float("yshot");
        if (!se_have_par("xshot"))
            s[0][2] = o3 + 0.5 * (n3 - 1) * d3;
        else
            s[0][2] = se_get_par_float("xshot");

        oshot = s[0][1];
        dshot = -s[0][1] + (o2 + 0.5 * (n2 - 1) * d2);
        INFO(("Shooting from zshot=%g yshot=%g xshot=%g",
              s[0][0], s[0][1], s[0][2]));
    }

    // sep_copy_headers(time, vel);

    if (n1 <= 0 || n2 <= 0 || n3 <= 0)
        ERROR(("Invalid input dimensions n1=%d n2=%d n3=%d", n1, n2, n3));

    n123 = (size_t)n1 * (size_t)n2 * (size_t)n3;
    if (n123 > (size_t)std::numeric_limits<int>::max())
        ERROR(("Grid is too large for the eikods integer indexing: n1*n2*n3=%zu", n123));

    v = alloc1float(n123);

    se_fsio_read_float(vel->data->io, v, n123);
    if (isvel)
    {
        /* transform velocity to slowness squared */
#ifdef SE_USE_OMP
#pragma omp parallel for
#endif
        for (i = 0; i < (int)n123; i++)
        {
            float slow_i = v[i];
            v[i] = 1. / (slow_i * slow_i);
        }
    }

    /* first-order derivative */
    if (!se_have_par("tdl1"))
        ERROR(("Output derivative tdl1= missing."));
    else
        tdl1_f = se_get_par_str("tdl1");
    if (!se_have_par("tds1"))
        ERROR(("Output derivative tds1= missing."));
    else
        tds1_f = se_get_par_str("tds1");
    tdl1 = sep_open(tdl1_f, SEP_WRITE, 0);
    tds1 = sep_open(tds1_f, SEP_WRITE, 0);
    tdl1->headers->ndim = 3;
    tds1->headers->ndim = 3;

    /* Optional branch-aware spatial gradients of the traveltime, computed by
       an upwind (ray-backtracing) stencil: gradz=dT/dz, gradx=dT/dx. These
       feed the scattering-angle classification directly. */
    if (se_have_par("gradz"))
    {
        gradz_f = se_get_par_str("gradz");
        gradz = sep_open(gradz_f, SEP_WRITE, 0);
        gradz->headers->ndim = 3;
    }
    if (se_have_par("gradx"))
    {
        gradx_f = se_get_par_str("gradx");
        gradx = sep_open(gradx_f, SEP_WRITE, 0);
        gradx->headers->ndim = 3;
    }

    /* second-order derivative */
    if (se_have_par("tdl2"))
    {
        tdl2_f = se_get_par_str("tdl2");
        tdl2 = sep_open(tdl2_f, SEP_WRITE, 0);
        tdl2->headers->ndim = 3;
    }
    if (se_have_par("tds2"))
    {
        tds2_f = se_get_par_str("tds2");
        tds2 = sep_open(tds2_f, SEP_WRITE, 0);
        tds2->headers->ndim = 3;
    }

    if (!se_have_par("l"))
        l = 1;
    else
        l = se_get_par_int("l");
    /* source perturbation direction */

    time->headers->n[0] = n1;
    time->headers->d[0] = d1;
    time->headers->o[0] = o1;
    time->headers->n[1] = n2;
    time->headers->d[1] = d2;
    time->headers->o[1] = o2;
    time->headers->n[2] = nshot;
    time->headers->d[2] = dshot;
    time->headers->o[2] = oshot;

    tdl1->headers->n[0] = n1;
    tdl1->headers->d[0] = d1;
    tdl1->headers->o[0] = o1;
    tdl1->headers->n[1] = n2;
    tdl1->headers->d[1] = d2;
    tdl1->headers->o[1] = o2;
    tdl1->headers->n[2] = nshot;
    tdl1->headers->d[2] = dshot;
    tdl1->headers->o[2] = oshot;

    tds1->headers->n[0] = n1;
    tds1->headers->d[0] = d1;
    tds1->headers->o[0] = o1;
    tds1->headers->n[1] = n2;
    tds1->headers->d[1] = d2;
    tds1->headers->o[1] = o2;
    tds1->headers->n[2] = nshot;
    tds1->headers->d[2] = dshot;
    tds1->headers->o[2] = oshot;

    if (gradz != NULL)
    {
        gradz->headers->n[0] = n1;
        gradz->headers->d[0] = d1;
        gradz->headers->o[0] = o1;
        gradz->headers->n[1] = n2;
        gradz->headers->d[1] = d2;
        gradz->headers->o[1] = o2;
        gradz->headers->n[2] = nshot;
        gradz->headers->d[2] = dshot;
        gradz->headers->o[2] = oshot;
    }
    if (gradx != NULL)
    {
        gradx->headers->n[0] = n1;
        gradx->headers->d[0] = d1;
        gradx->headers->o[0] = o1;
        gradx->headers->n[1] = n2;
        gradx->headers->d[1] = d2;
        gradx->headers->o[1] = o2;
        gradx->headers->n[2] = nshot;
        gradx->headers->d[2] = dshot;
        gradx->headers->o[2] = oshot;
    }

    if (tdl2 != NULL)
    {
        tdl2->headers->n[0] = n1;
        tdl2->headers->d[0] = d1;
        tdl2->headers->o[0] = o1;
        tdl2->headers->n[1] = n2;
        tdl2->headers->d[1] = d2;
        tdl2->headers->o[1] = o2;
        tdl2->headers->n[2] = nshot;
        tdl2->headers->d[2] = dshot;
        tdl2->headers->o[2] = oshot;
    }
    if (tds2 != NULL)
    {
        tds2->headers->n[0] = n1;
        tds2->headers->d[0] = d1;
        tds2->headers->o[0] = o1;
        tds2->headers->n[1] = n2;
        tds2->headers->d[1] = d2;
        tds2->headers->o[1] = o2;
        tds2->headers->n[2] = nshot;
        tds2->headers->d[2] = dshot;
        tds2->headers->o[2] = oshot;
    }

    preallocate_float_output(time, n123, nshot, "time");
    preallocate_float_output(tdl1, n123, nshot, "tdl1");
    preallocate_float_output(tds1, n123, nshot, "tds1");
    if (gradz != NULL)
        preallocate_float_output(gradz, n123, nshot, "gradz");
    if (gradx != NULL)
        preallocate_float_output(gradx, n123, nshot, "gradx");
    if (tdl2 != NULL)
        preallocate_float_output(tdl2, n123, nshot, "tdl2");
    if (tds2 != NULL)
        preallocate_float_output(tds2, n123, nshot, "tds2");

    const size_t nout = 3 + (tdl2 != NULL ? 1 : 0) + (tds2 != NULL ? 1 : 0);
    if ((size_t)nshot > std::numeric_limits<size_t>::max() / n123)
        ERROR(("Total output sample-count overflow: nshot=%d n123=%zu", nshot, n123));
    const size_t samples_total = n123 * (size_t)nshot;
    if (samples_total > std::numeric_limits<size_t>::max() / sizeof(float))
        ERROR(("Output byte-count overflow: samples=%zu", samples_total));
    const size_t bytes_per_output = samples_total * sizeof(float);
    if (nout > std::numeric_limits<size_t>::max() / bytes_per_output)
        ERROR(("Buffered output byte-count overflow: outputs=%zu bytes_per_output=%zu", nout, bytes_per_output));
    const size_t buffered_bytes = bytes_per_output * nout;
    int buffer_mb = se_have_par("eikods_buffer_mb") ? se_get_par_int("eikods_buffer_mb") : 4096;
    if (buffer_mb < 0)
        buffer_mb = 0;
    size_t buffer_limit = (size_t)buffer_mb * 1024u * 1024u;
    const bool buffered_output = (buffer_limit > 0 && buffered_bytes <= buffer_limit);

    std::vector<float> time_all, gzall, gxall, dl1_all, ds1_all, dl2_all, ds2_all;
    if (buffered_output)
    {
        time_all.resize(samples_total);
        dl1_all.resize(samples_total);
        ds1_all.resize(samples_total);
        if (gradz != NULL)
            gzall.resize(samples_total);
        if (gradx != NULL)
            gxall.resize(samples_total);
        if (tdl2 != NULL)
            dl2_all.resize(samples_total);
        if (tds2 != NULL)
            ds2_all.resize(samples_total);
        INFO(("buffered run: %d shots, %zu samples/shot, %.1f MiB output buffer",
              nshot, n123, buffered_bytes / (1024.0 * 1024.0)));
    }
    else
    {
        INFO(("streaming run: %d shots, %zu samples/shot (requires %.1f MiB > eikods_buffer_mb=%d)",
              nshot, n123, buffered_bytes / (1024.0 * 1024.0), buffer_mb));
    }

    auto t_start = std::chrono::steady_clock::now();
#ifdef SE_USE_OMP
    INFO(("OpenMP enabled: computing shots with up to %d threads%s",
          omp_get_max_threads(), buffered_output ? " and deferred output writes" : " and random-access writes"));
#pragma omp parallel
    {
        if (!efmm)
            eikods_init(n3, n2, n1);
        else
            efmm_eikods_init(n3, n2, n1);

        std::vector<float> t(n123), dl1(n123), ds1(n123);
        std::vector<float> gz, gx;         /* branch-aware spatial gradients */
        if (gradz != NULL)
            gz.resize(n123);
        if (gradx != NULL)
            gx.resize(n123);
        std::vector<float> dl2_buf, ds2_buf;
        if (tdl2 != NULL || tds2 != NULL)
        {
            dl2_buf.resize(n123);
            ds2_buf.resize(n123);
        }
        std::vector<int> p(n123);

#pragma omp for schedule(dynamic)
        for (is = 0; is < nshot; is++)
        {
            int tid = omp_get_thread_num();
            if (is == 0 || is % 50 == 0 || is == nshot - 1)
                INFO(("shot %d/%d start on thread %d", is + 1, nshot, tid));

            if (efmm)
                efmm_eikods(t.data(), v, p.data(), plane,
                            n3, n2, n1, o3, o2, o1, d3, d2, d1,
                            s[is][2], s[is][1], s[is][0], b3, b2, b1,
                            order, l, dl1.data(), ds1.data(),
                            dl2_buf.empty() ? NULL : dl2_buf.data(),
                            ds2_buf.empty() ? NULL : ds2_buf.data());
            else
                eikods(t.data(), v, p.data(), plane,
                       n3, n2, n1, o3, o2, o1, d3, d2, d1,
                       s[is][2], s[is][1], s[is][0], b3, b2, b1,
                       order, l, dl1.data(), ds1.data(),
                       dl2_buf.empty() ? NULL : dl2_buf.data(),
                       ds2_buf.empty() ? NULL : ds2_buf.data());

            /* Branch-aware spatial gradients via upwind (ray-backtracing)
               stencils. At a first-arrival crease each component uses only
               the smaller-traveltime neighbor, so it follows the correct ray
               leg instead of straddling it like a centered FD. */
            if (gradz != NULL || gradx != NULL)
            {
                const int nn[3] = {n1, n2, n3};
                const int ss[3] = {1, n1, n1 * n2};
                const float dd[3] = {d1, d2, d3};
                for (int q = 0; q < (int)n123; q++)
                {
                    if (gradz != NULL)
                        gz[q] = upwind_grad_axis(t.data(), q, 0, nn, ss, dd);
                    if (gradx != NULL)
                        gx[q] = upwind_grad_axis(t.data(), q, 1, nn, ss, dd);
                }
            }

            if (buffered_output)
            {
                const size_t off = (size_t)is * n123;
                std::copy(t.begin(), t.end(), time_all.begin() + off);
                std::copy(dl1.begin(), dl1.end(), dl1_all.begin() + off);
                std::copy(ds1.begin(), ds1.end(), ds1_all.begin() + off);
                if (gradz != NULL)
                    std::copy(gz.begin(), gz.end(), gzall.begin() + off);
                if (gradx != NULL)
                    std::copy(gx.begin(), gx.end(), gxall.begin() + off);
                if (tdl2 != NULL)
                    std::copy(dl2_buf.begin(), dl2_buf.end(), dl2_all.begin() + off);
                if (tds2 != NULL)
                    std::copy(ds2_buf.begin(), ds2_buf.end(), ds2_all.begin() + off);
            }
            else
            {
                write_float_shot_at(time, t.data(), n123, is, "time");
                write_float_shot_at(tdl1, dl1.data(), n123, is, "tdl1");
                write_float_shot_at(tds1, ds1.data(), n123, is, "tds1");
                if (gradz != NULL)
                    write_float_shot_at(gradz, gz.data(), n123, is, "gradz");
                if (gradx != NULL)
                    write_float_shot_at(gradx, gx.data(), n123, is, "gradx");
                if (tdl2 != NULL)
                    write_float_shot_at(tdl2, dl2_buf.data(), n123, is, "tdl2");
                if (tds2 != NULL)
                    write_float_shot_at(tds2, ds2_buf.data(), n123, is, "tds2");
            }

            if (is == 0 || is % 50 == 0 || is == nshot - 1)
                INFO(("shot %d/%d done on thread %d", is + 1, nshot, tid));
        }
        if (!efmm)
            eikods_close();
    }
#else
    if (!efmm)
        eikods_init(n3, n2, n1);
    else
        efmm_eikods_init(n3, n2, n1);

    std::vector<float> t(n123), dl1(n123), ds1(n123);
    std::vector<float> gz, gx;         /* branch-aware spatial gradients */
    if (gradz != NULL)
        gz.resize(n123);
    if (gradx != NULL)
        gx.resize(n123);
    std::vector<float> dl2_buf, ds2_buf;
    if (tdl2 != NULL || tds2 != NULL)
    {
        dl2_buf.resize(n123);
        ds2_buf.resize(n123);
    }
    std::vector<int> p(n123);

    for (is = 0; is < nshot; is++)
    {
        if (is == 0 || is % 50 == 0 || is == nshot - 1)
            INFO(("shot %d/%d start", is + 1, nshot));

        if (efmm)
            efmm_eikods(t.data(), v, p.data(), plane,
                        n3, n2, n1, o3, o2, o1, d3, d2, d1,
                        s[is][2], s[is][1], s[is][0], b3, b2, b1,
                        order, l, dl1.data(), ds1.data(),
                        dl2_buf.empty() ? NULL : dl2_buf.data(),
                        ds2_buf.empty() ? NULL : ds2_buf.data());
        else
            eikods(t.data(), v, p.data(), plane,
                   n3, n2, n1, o3, o2, o1, d3, d2, d1,
                   s[is][2], s[is][1], s[is][0], b3, b2, b1,
                   order, l, dl1.data(), ds1.data(),
                   dl2_buf.empty() ? NULL : dl2_buf.data(),
                   ds2_buf.empty() ? NULL : ds2_buf.data());

        /* Branch-aware spatial gradients via upwind (ray-backtracing). */
        if (gradz != NULL || gradx != NULL)
        {
            const int nn[3] = {n1, n2, n3};
            const int ss[3] = {1, n1, n1 * n2};
            const float dd[3] = {d1, d2, d3};
            for (int q = 0; q < (int)n123; q++)
            {
                if (gradz != NULL)
                    gz[q] = upwind_grad_axis(t.data(), q, 0, nn, ss, dd);
                if (gradx != NULL)
                    gx[q] = upwind_grad_axis(t.data(), q, 1, nn, ss, dd);
            }
        }

        if (buffered_output)
        {
            const size_t off = (size_t)is * n123;
            std::copy(t.begin(), t.end(), time_all.begin() + off);
            std::copy(dl1.begin(), dl1.end(), dl1_all.begin() + off);
            std::copy(ds1.begin(), ds1.end(), ds1_all.begin() + off);
            if (gradz != NULL)
                std::copy(gz.begin(), gz.end(), gzall.begin() + off);
            if (gradx != NULL)
                std::copy(gx.begin(), gx.end(), gxall.begin() + off);
            if (tdl2 != NULL)
                std::copy(dl2_buf.begin(), dl2_buf.end(), dl2_all.begin() + off);
            if (tds2 != NULL)
                std::copy(ds2_buf.begin(), ds2_buf.end(), ds2_all.begin() + off);
        }
        else
        {
            write_float_shot_at(time, t.data(), n123, is, "time");
            write_float_shot_at(tdl1, dl1.data(), n123, is, "tdl1");
            write_float_shot_at(tds1, ds1.data(), n123, is, "tds1");
            if (gradz != NULL)
                write_float_shot_at(gradz, gz.data(), n123, is, "gradz");
            if (gradx != NULL)
                write_float_shot_at(gradx, gx.data(), n123, is, "gradx");
            if (tdl2 != NULL)
                write_float_shot_at(tdl2, dl2_buf.data(), n123, is, "tdl2");
            if (tds2 != NULL)
                write_float_shot_at(tds2, ds2_buf.data(), n123, is, "tds2");
        }
    }
    if (!efmm)
        eikods_close();
#endif

    if (buffered_output)
    {
        INFO(("writing buffered output"));
        write_float_all(time, time_all.data(), n123, nshot, "time");
        write_float_all(tdl1, dl1_all.data(), n123, nshot, "tdl1");
        write_float_all(tds1, ds1_all.data(), n123, nshot, "tds1");
        if (gradz != NULL)
            write_float_all(gradz, gzall.data(), n123, nshot, "gradz");
        if (gradx != NULL)
            write_float_all(gradx, gxall.data(), n123, nshot, "gradx");
        if (tdl2 != NULL)
            write_float_all(tdl2, dl2_all.data(), n123, nshot, "tdl2");
        if (tds2 != NULL)
            write_float_all(tds2, ds2_all.data(), n123, nshot, "tds2");
    }

    INFO(("FINISH."));
    auto t_end = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration<double>(t_end - t_start).count();
    INFO(("Done. Elapsed time: %.3f s.", elapsed_seconds));

    sep_close(vel);
    sep_close(time);
    sep_close(tdl1);
    sep_close(tds1);
    if (gradz != NULL)
        sep_close(gradz);
    if (gradx != NULL)
        sep_close(gradx);
    if (tdl2 != NULL)
        sep_close(tdl2);
    if (tds2 != NULL)
        sep_close(tds2);

    return 0;
}

/* 	$Id: Meikonal.c 7107 2011-04-10 02:04:14Z ivlad $	 */
