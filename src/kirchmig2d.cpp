#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>
#include "./help/kirchmig2d_help.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ----------------------------------------------------------------------
   Early-time mute with cosine ramp.

   tmute = tmute0 + offset / vmute

   Notes:
   1. If tmute0 < 0 and vmute <= 0, mute is disabled.
   2. off, vmute, coordinate units must be consistent.
      For example:
        x in m  -> vmute in m/s
        x in km -> vmute in km/s
---------------------------------------------------------------------- */
static void apply_early_mute(
    float *trace,
    int nt,
    float dt,
    float t0,
    float off,
    float tmute0,
    float vmute,
    float mutewidth)
{
    if (trace == NULL || nt <= 0 || dt <= 0.0f) return;

    if (tmute0 < 0.0f && vmute <= 0.0f) return;

    float tmute = 0.0f;

    if (tmute0 > 0.0f)
        tmute += tmute0;

    if (vmute > 0.0f)
        tmute += off / vmute;

    if (tmute <= t0) return;

    for (int it = 0; it < nt; it++)
    {
        float tt = t0 + it * dt;

        if (tt <= tmute)
        {
            trace[it] = 0.0f;
        }
        else if (mutewidth > 0.0f && tt < tmute + mutewidth)
        {
            float u = (tt - tmute) / mutewidth;
            if (u < 0.0f) u = 0.0f;
            if (u > 1.0f) u = 1.0f;

            /* 0 -> 1 cosine ramp */
            float w = 0.5f * (1.0f - cosf((float)M_PI * u));
            trace[it] *= w;
        }
    }
}

/* ----------------------------------------------------------------------
   One-sided cosine taper near the aperture boundary.

   dist <= (1 - taper_frac) * aper_halfwidth : weight = 1
   dist >= aper_halfwidth                    : weight = 0
   middle zone                               : cosine taper

   taper_frac = 0.2 means the outer 20% of the aperture is tapered.
---------------------------------------------------------------------- */
static float one_aperture_taper(
    float dist,
    float aper_halfwidth,
    float taper_frac)
{
    if (aper_halfwidth <= 0.0f) return 0.0f;

    if (dist >= aper_halfwidth) return 0.0f;

    if (taper_frac <= 0.0f) return 1.0f;

    if (taper_frac > 0.95f) taper_frac = 0.95f;

    float aper_inner = (1.0f - taper_frac) * aper_halfwidth;

    if (dist <= aper_inner) return 1.0f;

    float denom = aper_halfwidth - aper_inner;
    if (denom <= 0.0f) return 1.0f;

    float u = (dist - aper_inner) / denom;
    if (u < 0.0f) u = 0.0f;
    if (u > 1.0f) u = 1.0f;

    return 0.5f * (1.0f + cosf((float)M_PI * u));
}

/* ----------------------------------------------------------------------
   Source-side and receiver-side aperture weight.

   If aperture >= 89.9 deg, aperture checking is effectively disabled.
---------------------------------------------------------------------- */
static float aperture_weight(
    float ximg,
    float zimg,
    float xsrc,
    float xrec,
    float aper_deg,
    float taper_frac)
{
    if (zimg <= 0.0f) return 0.0f;

    if (aper_deg >= 89.9f)
        return 1.0f;

    if (aper_deg <= 0.0f)
        return 0.0f;

    float aper_rad = aper_deg * (float)M_PI / 180.0f;
    float aper_halfwidth = zimg * tanf(aper_rad);

    if (aper_halfwidth <= 0.0f)
        return 0.0f;

    float dxs = fabsf(ximg - xsrc);
    float dxr = fabsf(ximg - xrec);

    float ws = one_aperture_taper(dxs, aper_halfwidth, taper_frac);
    if (ws <= 0.0f) return 0.0f;

    float wr = one_aperture_taper(dxr, aper_halfwidth, taper_frac);
    if (wr <= 0.0f) return 0.0f;

    return ws * wr;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        kirchmig2d_help::print_help();
        return 1;
    }

    const char *unit = NULL, *type = "hermit";

    int cig = 0;
    int cmp = 0;
    int normalize = 0;

    off_t nzx;

    int nt = 0, nx = 0, sny = 0, rny = 0, ns = 0, nh = 0, nz = 0;
    int i, ix, iz, ih, is, ist, iht, ng, ithr, nthr;

    float *trace = NULL;
    float **traces = NULL;
    float **out = NULL;
    float **wsum = NULL;

    float **stbl = NULL;
    float **rtbl = NULL;
    float *stable = NULL;
    float *rtable = NULL;

    float **stblx = NULL;
    float **rtblx = NULL;
    float *stablex = NULL;
    float *rtablex = NULL;

    float ds = 0.0f, s0 = 0.0f, x0 = 0.0f;
    float sy0 = 0.0f, sdy = 0.0f;
    float ry0 = 0.0f, rdy = 0.0f;
    float s, h, h0 = 0.0f, dh = 0.0f;
    float dx = 0.0f, dz = 0.0f, z0 = 0.0f;
    float ti, t0 = 0.0f, t1, t2, dt = 0.0f, tau = 0.0f;
    float aal = 1.0f, tx, aper = 90.0f;

    /* New artifact-suppression parameters */
    float offmax = -1.0f;        /* <=0: disabled */
    float tmute0 = -1.0f;        /* <0 and vmute<=0: disabled */
    float vmute = 0.0f;          /* <=0: no offset-dependent mute */
    float mutewidth = 0.02f;     /* cosine ramp width, in seconds */
    float aper_taper = 0.20f;    /* outer aperture taper fraction */

    const char *data_f = NULL, *mig_f = NULL;
    const char *stim_f = NULL, *sder_f = NULL;
    const char *rtim_f = NULL, *rder_f = NULL;

    sep_t *dat = NULL, *mig = NULL;
    sep_t *stim = NULL, *sder = NULL;
    sep_t *rtim = NULL, *rder = NULL;

    long long n_s_clamp = 0;
    long long n_r_clamp = 0;
    long long n_trace_offmute = 0;
    long long n_aper_skip = 0;
    long long n_t_out = 0;
    long long n_t_in = 0;

    se_par_init(argc, argv);

    if (!se_have_par("cmp"))
        cmp = 0;
    else
        cmp = se_get_par_int("cmp");

    if (!se_have_par("data"))
        ERROR(("Need data="));
    else
        data_f = se_get_par_str("data");

    if (!se_have_par("mig"))
        ERROR(("Need mig="));
    else
        mig_f = se_get_par_str("mig");

    dat = sep_open(data_f, SEP_READ, 0);
    mig = sep_open(mig_f, SEP_WRITE, 0);

    nt = dat->headers->n[0];
    t0 = dat->headers->o[0];
    dt = dat->headers->d[0];

    nh = dat->headers->n[1];
    h0 = dat->headers->o[1];
    dh = dat->headers->d[1];

    ns = dat->headers->n[2];
    s0 = dat->headers->o[2];
    ds = dat->headers->d[2];

    printf("data: nt=%d, dt=%g, t0=%g; nh=%d, dh=%g, h0=%g; ns=%d, ds=%g, s0=%g\n",
           nt, dt, t0, nh, dh, h0, ns, ds, s0);

    if (1 == nh) dh = 0.0f;
    if (1 == ns) ds = 0.0f;

    if (!se_have_par("stable"))
        ERROR(("Need stable="));
    else
    {
        stim_f = se_get_par_str("stable");
        stim = sep_open(stim_f, SEP_READ, 0);
    }

    if (!se_have_par("sderiv"))
        ERROR(("Need sderiv="));
    else
    {
        sder_f = se_get_par_str("sderiv");
        sder = sep_open(sder_f, SEP_READ, 0);
    }

    nz = stim->headers->n[0];
    z0 = (float)stim->headers->o[0];
    dz = (float)stim->headers->d[0];

    nx = stim->headers->n[1];
    x0 = (float)stim->headers->o[1];
    dx = (float)stim->headers->d[1];

    nzx = (off_t)nz * (off_t)nx;

    sny = stim->headers->n[2];
    sy0 = (float)stim->headers->o[2];
    sdy = (float)stim->headers->d[2];

    stbl = alloc2float(nzx, sny);
    se_fsio_read_float(stim->data->io, stbl[0], nzx * sny);
    sep_close(stim);

    stblx = alloc2float(nzx, sny);
    se_fsio_read_float(sder->data->io, stblx[0], nzx * sny);
    sep_close(sder);

    printf("sy0 = %f, sdy = %f\n", sy0, sdy);
    printf("source time table: nx = %d, nz = %d, ns = %d\n", nx, nz, sny);
    printf("dsource time table: nx = %d, nz = %d, ds = %f\n", nx, nz, sdy);

    if (!se_have_par("rtable"))
        ERROR(("Need rtable="));
    else
    {
        rtim_f = se_get_par_str("rtable");
        rtim = sep_open(rtim_f, SEP_READ, 0);
    }

    if (!se_have_par("rderiv"))
        ERROR(("Need rderiv="));
    else
    {
        rder_f = se_get_par_str("rderiv");
        rder = sep_open(rder_f, SEP_READ, 0);
    }

    rny = rtim->headers->n[2];
    ry0 = (float)rtim->headers->o[2];
    rdy = (float)rtim->headers->d[2];

    rtbl = alloc2float(nzx, rny);
    se_fsio_read_float(rtim->data->io, rtbl[0], nzx * rny);
    sep_close(rtim);

    rtblx = alloc2float(nzx, rny);
    se_fsio_read_float(rder->data->io, rtblx[0], nzx * rny);
    sep_close(rder);

    printf("ry0 = %f, rdy = %f\n", ry0, rdy);
    printf("receiver time table: nx = %d, nz = %d, nr = %d\n", nx, nz, rny);
    printf("dreceiver time table: nx = %d, nz = %d, dr = %f\n", nx, nz, rdy);

    if (!se_have_par("tau"))
        tau = 0.0f;
    else
        tau = se_get_par_float("tau");
    /* static time-shift, in seconds */

    if (!se_have_par("aperture"))
        aper = 90.0f;
    else
        aper = se_get_par_float("aperture");
    /* migration aperture, in degree */

    if (!se_have_par("antialias"))
        aal = 1.0f;
    else
        aal = se_get_par_float("antialias");
    /* antialiasing coefficient */

    if (!se_have_par("cig"))
        cig = 0;
    else
        cig = se_get_par_int("cig");
    /* y - output common offset gathers */

    if (!se_have_par("type"))
        type = "hermit";
    else
        type = se_get_par_str("type");
    /* interpolation type: linear / partial / hermit */

    /* New parameters */
    if (!se_have_par("offmax"))
        offmax = -1.0f;
    else
        offmax = se_get_par_float("offmax");
    /* maximum absolute offset. offmax<=0 means disabled */

    if (!se_have_par("tmute0"))
        tmute0 = -1.0f;
    else
        tmute0 = se_get_par_float("tmute0");
    /* zero-offset early mute time, in seconds. tmute0<0 and vmute<=0 means disabled */

    if (!se_have_par("vmute"))
        vmute = 0.0f;
    else
        vmute = se_get_par_float("vmute");
    /* mute velocity. Units must be consistent with coordinates */

    if (!se_have_par("mutewidth"))
        mutewidth = 0.02f;
    else
        mutewidth = se_get_par_float("mutewidth");
    /* mute cosine ramp width, in seconds */

    if (!se_have_par("aper_taper"))
        aper_taper = 0.20f;
    else
        aper_taper = se_get_par_float("aper_taper");
    /* aperture taper fraction */

    if (!se_have_par("normalize"))
        normalize = 0;
    else
        normalize = se_get_par_int("normalize");
    /* normalize image by accumulated aperture weights */

    printf("parameters: cmp=%d, cig=%d, type=%s\n", cmp, cig, type);
    printf("parameters: aperture=%g deg, aper_taper=%g, antialias=%g\n",
           aper, aper_taper, aal);
    printf("parameters: offmax=%g, tmute0=%g, vmute=%g, mutewidth=%g, normalize=%d\n",
           offmax, tmute0, vmute, mutewidth, normalize);

    ng = cig ? nh : 1;

    sep_set_axis(mig, 0, nz, z0, dz, "Depth");
    sep_set_axis(mig, 1, nx, x0, dx, "Lateral");

    unit = sep_get_hdr(dat, "unit2", NULL);
    if (NULL != unit)
        sep_set_header(mig, "unit1", unit);

    if (cig)
    {
        sep_set_axis(mig, 2, nh, h0, dh, cmp ? "Offset" : "Receiver");

        if (NULL != unit)
            sep_set_header(mig, "unit3", unit);
    }
    else
    {
        sep_set_header_int(mig, "n3", 1);
    }

    sep_write_headers(mig);

    /* allocate temporary memory */
    out = alloc2float(nzx, ng);
    wsum = alloc2float(nzx, ng);

    trace = alloc1float(nt);

    stable = alloc1float(nzx);
    stablex = alloc1float(nzx);
    rtable = alloc1float(nzx);
    rtablex = alloc1float(nzx);

    memset(out[0], 0, nzx * ng * sizeof(float));
    memset(wsum[0], 0, nzx * ng * sizeof(float));

    /* initialize interpolation */
    tinterp_init(nzx, sdy, rdy);

    /* initialize summation */
    kirmig_init(nt, dt, t0);

    /* get number of threads */
#ifdef SE_USE_OMP
    nthr = 1;
#pragma omp parallel
    {
#pragma omp single
        {
            nthr = omp_get_num_threads();
        }
    }
#else
    nthr = 1;
#endif

    WARN((">>Using %d threads<<\n", nthr));

    /* Correct allocation: one pointer per thread */
    traces = (float **)calloc(nthr, sizeof(float *));

    for (is = 0; is < ns; is++)
    {
        /* shot coordinate */
        s = s0 + is * ds;
        WARN(("shot %d of %d, s=%g", is + 1, ns, s));

        /* source traveltime interpolation */
        if (sny <= 1 || fabsf(sdy) <= 0.0f)
        {
            for (i = 0; i < nzx; i++)
            {
                stable[i] = stbl[0][i];
                stablex[i] = stblx[0][i];
            }
        }
        else
        {
            float fsidx = (s - sy0) / sdy;
            ist = (int)floorf(fsidx);

            if (ist < 0)
            {
                n_s_clamp++;
                for (i = 0; i < nzx; i++)
                {
                    stable[i] = stbl[0][i];
                    stablex[i] = stblx[0][i];
                }
            }
            else if (ist >= sny - 1)
            {
                n_s_clamp++;
                for (i = 0; i < nzx; i++)
                {
                    stable[i] = stbl[sny - 1][i];
                    stablex[i] = stblx[sny - 1][i];
                }
            }
            else
            {
                float soff = s - (sy0 + ist * sdy);

                switch (type[0])
                {
                case 'l':
                    tinterp_linear(true, stable, soff, stbl[ist], stbl[ist + 1]);
                    dinterp_linear(true, stablex, soff, stbl[ist], stbl[ist + 1]);
                    break;

                case 'p':
                    tinterp_partial(true, stable, soff, nz, nx, dx, stbl[ist], stbl[ist + 1]);
                    dinterp_partial(true, stablex, soff, nz, nx, dx, stbl[ist], stbl[ist + 1]);
                    break;

                case 'h':
                default:
                    tinterp_hermite(true, stable, soff, stbl[ist], stbl[ist + 1], stblx[ist], stblx[ist + 1]);
                    dinterp_hermite(true, stablex, soff, stbl[ist], stbl[ist + 1], stblx[ist], stblx[ist + 1]);
                    break;
                }
            }
        }

        for (ih = 0; ih < nh; ih++)
        {
            h = h0 + ih * dh;

            /* Receiver coordinate.
               cmp=0: h is absolute receiver coordinate.
               cmp=1: h is offset, receiver = source + offset.
            */
            float xsrc = s;
            float xrec = cmp ? (s + h) : h;
            float off = fabsf(xrec - xsrc);

            /* receiver traveltime interpolation */
            if (rny <= 1 || fabsf(rdy) <= 0.0f)
            {
                for (i = 0; i < nzx; i++)
                {
                    rtable[i] = rtbl[0][i];
                    rtablex[i] = rtblx[0][i];
                }
            }
            else
            {
                float rcoord = xrec;
                float fridx = (rcoord - ry0) / rdy;
                iht = (int)floorf(fridx);

                if (iht < 0)
                {
                    n_r_clamp++;
                    for (i = 0; i < nzx; i++)
                    {
                        rtable[i] = rtbl[0][i];
                        rtablex[i] = rtblx[0][i];
                    }
                }
                else if (iht >= rny - 1)
                {
                    n_r_clamp++;
                    for (i = 0; i < nzx; i++)
                    {
                        rtable[i] = rtbl[rny - 1][i];
                        rtablex[i] = rtblx[rny - 1][i];
                    }
                }
                else
                {
                    float roff = rcoord - (ry0 + iht * rdy);

                    switch (type[0])
                    {
                    case 'l':
                        tinterp_linear(false, rtable, roff, rtbl[iht], rtbl[iht + 1]);
                        dinterp_linear(false, rtablex, roff, rtbl[iht], rtbl[iht + 1]);
                        break;

                    case 'p':
                        tinterp_partial(false, rtable, roff, nz, nx, dx, rtbl[iht], rtbl[iht + 1]);
                        dinterp_partial(false, rtablex, roff, nz, nx, dx, rtbl[iht], rtbl[iht + 1]);
                        break;

                    case 'h':
                    default:
                        tinterp_hermite(false, rtable, roff, rtbl[iht], rtbl[iht + 1], rtblx[iht], rtblx[iht + 1]);
                        dinterp_hermite(false, rtablex, roff, rtbl[iht], rtbl[iht + 1], rtblx[iht], rtblx[iht + 1]);
                        break;
                    }
                }
            }

            /* read trace.
               Important: even if this trace is muted by offmax, we must read it
               because data are read sequentially as [nt, nh, ns].
            */
            se_fsio_read_float(dat->data->io, trace, nt);

            /* offset mute */
            if (offmax > 0.0f && off > offmax)
            {
                n_trace_offmute++;
                continue;
            }

            /* early-time mute before double integration */
            apply_early_mute(trace, nt, dt, t0, off, tmute0, vmute, mutewidth);

            /* double integration for antialiasing */
            doubint(nt, trace);

            for (ithr = 0; ithr < nthr; ithr++)
                traces[ithr] = trace;

#ifdef SE_USE_OMP
#pragma omp parallel for private(iz, ix, t1, t2, ti, tx, ithr) reduction(+:n_aper_skip,n_t_out,n_t_in)
#endif
            for (i = 0; i < nzx; i++)
            {
#ifdef SE_USE_OMP
                ithr = omp_get_thread_num();
#else
                ithr = 0;
#endif
                iz = i % nz;
                ix = (i - iz) / nz;

                float ximg = x0 + ix * dx;
                float zimg = z0 + iz * dz;

                /* corrected aperture with true depth z0 + iz*dz */
                float waper = aperture_weight(ximg, zimg, xsrc, xrec, aper, aper_taper);

                if (waper <= 0.0f)
                {
                    n_aper_skip++;
                    continue;
                }

                t1 = stable[i];
                t2 = rtable[i];
                ti = t1 + t2 + tau;

                /* skip invalid sampling time */
                if (ti < t0 || ti > t0 + (nt - 1) * dt)
                {
                    n_t_out++;
                    continue;
                }

                n_t_in++;

                tx = fmaxf(fabsf(stablex[i] * ds), fabsf(rtablex[i] * dh));

                int ig = cig ? ih : 0;

                /* kirmig_pick directly accumulates into out.
                   To apply aperture taper without modifying library code:
                   1. save old value
                   2. call kirmig_pick
                   3. multiply the increment by waper
                */
                float before = out[ig][i];

                kirmig_pick(1, ti, tx * aal, out[ig] + i, traces[ithr]);

                out[ig][i] = before + waper * (out[ig][i] - before);

                wsum[ig][i] += waper;
            }
        }
    }

    WARN(("."));

    printf("\n---- migration statistics ----\n");
    printf("source table clamp count        = %lld\n", n_s_clamp);
    printf("receiver table clamp count      = %lld\n", n_r_clamp);
    printf("offset-muted trace count        = %lld\n", n_trace_offmute);
    printf("aperture-skipped sample count   = %lld\n", n_aper_skip);
    printf("valid time sample count         = %lld\n", n_t_in);
    printf("out-of-record time sample count = %lld\n", n_t_out);
    printf("--------------------------------\n");

    if (normalize)
    {
        for (int ig = 0; ig < ng; ig++)
        {
            for (i = 0; i < nzx; i++)
            {
                if (wsum[ig][i] > 1.0e-6f)
                    out[ig][i] /= wsum[ig][i];
            }
        }
    }

    se_fsio_write_float(mig->data->io, out[0], nzx * ng);

    sep_close(mig);
    sep_close(dat);

    return 0;
}