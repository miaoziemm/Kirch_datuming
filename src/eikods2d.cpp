#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>
#include "./help/eikods2d_help.h"

#include <limits>
#include <vector>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

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
    const char *sfile = NULL, *in_f = NULL, *out_f = NULL, *tdl1_f = NULL, *tds1_f = NULL, *tdl2_f = NULL, *tds2_f = NULL;
    bool isvel = true, plane[3] = {false, false, false};
    int efmm = 0;
    sep_t *vel = NULL, *time = NULL, *shots = NULL, *tdl1 = NULL, *tds1 = NULL, *tdl2 = NULL, *tds2 = NULL;

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

auto t_start = std::chrono::steady_clock::now();

    INFO(("streaming run: %d shots, %zu samples/shot", nshot, n123));
#ifdef SE_USE_OMP
    INFO(("OpenMP enabled: computing shots with up to %d threads and ordered streaming writes", omp_get_max_threads()));
#pragma omp parallel
    {
        if (!efmm)
            eikods_init(n3, n2, n1);
        else
            efmm_eikods_init(n3, n2, n1);

        std::vector<float> t(n123), dl1(n123), ds1(n123);
        std::vector<float> dl2_buf, ds2_buf;
        if (tdl2 != NULL || tds2 != NULL)
        {
            dl2_buf.resize(n123);
            ds2_buf.resize(n123);
        }
        std::vector<int> p(n123);

#pragma omp for schedule(dynamic) ordered
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

#pragma omp ordered
            {
                se_fsio_write_float(time->data->io, t.data(), n123);
                se_fsio_write_float(tdl1->data->io, dl1.data(), n123);
                se_fsio_write_float(tds1->data->io, ds1.data(), n123);
                if (tdl2 != NULL)
                    se_fsio_write_float(tdl2->data->io, dl2_buf.data(), n123);
                if (tds2 != NULL)
                    se_fsio_write_float(tds2->data->io, ds2_buf.data(), n123);
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

        se_fsio_write_float(time->data->io, t.data(), n123);
        se_fsio_write_float(tdl1->data->io, dl1.data(), n123);
        se_fsio_write_float(tds1->data->io, ds1.data(), n123);
        if (tdl2 != NULL)
            se_fsio_write_float(tdl2->data->io, dl2_buf.data(), n123);
        if (tds2 != NULL)
            se_fsio_write_float(tds2->data->io, ds2_buf.data(), n123);
    }
    if (!efmm)
        eikods_close();
#endif

    INFO(("FINISH."));
    auto t_end = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration<double>(t_end - t_start).count();
    INFO(("Done. Elapsed time: %.3f s.", elapsed_seconds));

    sep_close(vel);
    sep_close(time);
    sep_close(tdl1);
    sep_close(tds1);
    if (tdl2 != NULL)
        sep_close(tdl2);
    if (tds2 != NULL)
        sep_close(tds2);

    return 0;
}

/* 	$Id: Meikonal.c 7107 2011-04-10 02:04:14Z ivlad $	 */
