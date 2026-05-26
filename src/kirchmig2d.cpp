#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

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

void diff2(float *trace,int n, float d)
{
    int i;
    float *tmp = alloc1float(n);
    for (i=0; i < n; i++) {
    tmp[i] = trace[i];
    }
    for (i=1; i < n-1; i++) {
    trace[i] = (tmp[i+1]-tmp[i-1])/(2.*d);
    }
    trace[0] = (tmp[1]-tmp[0])/d;
    trace[n-1] = (tmp[n-1]-tmp[n-2])/d;
    free1float(tmp);
}


int main(int argc, char *argv[])
{
    se_par_init(argc, argv);
    char *unit;
    const char *type = NULL;
    int adj = 0, cig = 0, cmp = 0;
    off_t nzx = 0;
    int nt = 0, nx = 0, sny = 0, rny = 0, ns = 0, nh = 0, nz = 0, i = 0, ix = 0, iz = 0, ih = 0, is = 0, ist = 0, iht = 0, ng = 0, ithr = 0, nthr = 0;
    float *trace = NULL, **traces = NULL, **out = NULL, **stbl = NULL, **rtbl = NULL, *stable = NULL, *rtable = NULL, **stblx = NULL, **rtblx = NULL, *stablex = NULL, *rtablex = NULL;
    float ds = 0.0f, s0 = 0.0f, x0 = 0.0f, sy0 = 0.0f, sdy = 0.0f, ry0 = 0.0f, rdy = 0.0f, s = 0.0f, h = 0.0f, h0 = 0.0f, dh = 0.0f, dx = 0.0f, ti = 0.0f, t0 = 0.0f, t1 = 0.0f, t2 = 0.0f, dt = 0.0f, z0 = 0.0f, dz = 0.0f, tau = 0.0f;
    float aal = 0.0f, tx = 0.0f, aper = 0.0f;
    sep_t *dat = NULL, *mig = NULL, *stim = NULL, *sder = NULL, *rtim = NULL, *rder = NULL;
    char *dat_f = NULL, *mig_f = NULL, *stim_f = NULL, *sder_f = NULL, *rtim_f = NULL, *rder_f = NULL;

    int aperture_trace = -1; //用于限定成像道数，如果为-1则不限定

    if(!se_have_par("aperture_trace"))
        aperture_trace = -1;
    else
        aperture_trace = se_get_par_int("aperture_trace");

    INFO(("aperture_trace = %d\n", aperture_trace));

    if (!se_have_par("adj"))
        adj = 1;
    else
        adj = se_get_par_int("adj");
    if (!se_have_par("cig"))
        cig = 0;
    else
        cig = se_get_par_int("cig");
    if (!se_have_par("cmp"))
        cmp = 1;
    else
        cmp = se_get_par_int("cmp");

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

    if (adj)
    {
        dat = sep_open(dat_f, SEP_READ, 0);
        mig = sep_open(mig_f, SEP_WRITE, 0);
    }
    else
    {
        mig = sep_open(mig_f, SEP_READ, 0);
        dat = sep_open(dat_f, SEP_WRITE, 0);
    }


    if (adj)
    {
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
    }
    else
    {
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
            t0 = 0.0;
        else
            t0 = se_get_par_float("t0");
        if (!se_have_par("dt"))
            ERROR(("Need dt="));
        else
            dt = se_get_par_float("dt");
        if (!se_have_par("h0"))
            h0 = 0.0;
        else
            h0 = se_get_par_float("h0");
        if (!se_have_par("dh"))
            ERROR(("Need dh="));
        else
            dh = se_get_par_float("dh");
        if (!se_have_par("s0"))
            s0 = 0.0;
        else
            s0 = se_get_par_float("s0");
        if (!se_have_par("ds"))
            ERROR(("Need ds="));
        else
            ds = se_get_par_float("ds");
    }

    if (1 == nh)
        dh = 0.0;
    if (1 == ns)
        ds = 0.0;

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

    if (!se_have_par("tau"))
        tau = 0.0;
    else
        tau = se_get_par_float("tau"); /* static time-shift (in second) */


    // 打印tau
    INFO(("tau = %f\n", tau));
    
    if (!se_have_par("aperture"))
        aper = 90.0;
    else
        aper = se_get_par_float("aperture"); /* migration aperture (in degree) */
    if (!se_have_par("antialias"))
        aal = 1.0;
    else
        aal = se_get_par_float("antialias"); /* antialiasing */
    if (!se_have_par("cig"))
        cig = 0;
    else
        cig = se_get_par_int("cig"); /* y - output common offset gathers */

    ng = cig ? nh : 1;
    mig->headers->ndim = 3;

    if (adj)
    {
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
        if (cig)
        {
            mig->headers->n[2] = nh;
            mig->headers->o[2] = h0;
            mig->headers->d[2] = dh;
            sep_set_header(mig, "label3", cmp ? "Offset" : "Receiver");

            if (NULL != unit)
                sep_set_header(mig, "unit3", unit);
        }
        else
        {
            mig->headers->n[2] = 1;
        }
    }
    else
    {
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

    /* type of interpolation (default Hermit) */
    if (!se_have_par("type"))
        type = "hermit";
    else
        type = se_get_par_str("type");

    if(type[0] != 'l' && type[0] != 'p' && type[0] != 'h') {
        ERROR(("Unknown interpolation type: %s. Supported types are: linear (l), partial (p) and hermite (h)", type));
    }

    /* initialize interpolation */
    tinterp_init(nzx, sdy, rdy);

    /* initialize summation */
    kirmig_init(nt, dt, t0);

    if (adj)
    {
        memset(out[0], 0, nzx * ng * sizeof(float));
    }
    else
    {
        sep_read_fromsep(mig, nzx * ng, 0, &out[0], NULL);
    }

    /* determine number of threads without a data race */
#ifdef _OPENMP
    nthr = omp_get_max_threads();
#else
    nthr = 1;
#endif

    INFO((">>Using %d threads<<\n", nthr));

    if (!adj)
    {
        traces = alloc2float(nt, nthr);
    }

    for (is = 0; is < ns; is++)
    { /* shot */
        s = s0 + is * ds;
        INFO(("shot %d of %d;", is + 1, ns));

        /* cubic Hermite spline interpolation */
        ist = (s - sy0) / sdy;
        if (ist <= 0)
        {
            for (i = 0; i < nzx; i++)
            {
                stable[i] = stbl[0][i];
                stablex[i] = stblx[0][i];
            }
        }
        else if (ist >= sny - 1)
        {
            for (i = 0; i < nzx; i++)
            {
                stable[i] = stbl[sny - 1][i];
                stablex[i] = stblx[sny - 1][i];
            }
        }
        else
        {
            switch (type[0])
            {
            case 'l': /* linear */
                tinterp_linear(true, stable, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1]);
                dinterp_linear(true, stablex, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1]);
                break;

            case 'p': /* partial */
                tinterp_partial(true, stable, s - ist * sdy - sy0, nz, nx, dx, stbl[ist], stbl[ist + 1]);
                dinterp_partial(true, stablex, s - ist * sdy - sy0, nz, nx, dx, stbl[ist], stbl[ist + 1]);
                break;

            case 'h': /* hermit */
                tinterp_hermite(true, stable, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1], stblx[ist], stblx[ist + 1]);
                dinterp_hermite(true, stablex, s - ist * sdy - sy0, stbl[ist], stbl[ist + 1], stblx[ist], stblx[ist + 1]);
                break;
            }
        }

        for (ih = 0; ih < nh; ih++)
        { /* offset */
            h = h0 + ih * dh;
            /* cubic Hermite spline interpolation */
            iht = cmp ? (s + h - ry0) / rdy : (h - ry0) / rdy;


             if (adj)
            {
                /* read trace */
                se_fsio_read_float(dat->data->io, trace, nt);
                // 根据aperture_trace参数限定成像道数
                if (aperture_trace != -1 && fabs(iht-ist) >= aperture_trace) continue;
                doubint(nt, trace);

                // diff2(trace, nt, dt);

            }
            else
            {
                for (ithr = 0; ithr < nthr; ithr++)
                {
                    for (i = 0; i < nt; i++)
                    {
                        traces[ithr][i] = 0.;
                    }
                }
            }

            if (iht <= 0)
            {
                for (i = 0; i < nzx; i++)
                {
                    rtable[i] = rtbl[0][i];
                    rtablex[i] = rtblx[0][i];
                }
            }
            else if (iht >= rny - 1)
            {
                for (i = 0; i < nzx; i++)
                {
                    rtable[i] = rtbl[rny - 1][i];
                    rtablex[i] = rtblx[rny - 1][i];
                }
            }
            else
            {
                switch (type[0])
                {
                case 'l': /* linear */
                    tinterp_linear(false, rtable, cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0, rtbl[iht], rtbl[iht + 1]);
                    dinterp_linear(false, rtablex, cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0, rtbl[iht], rtbl[iht + 1]);
                    break;

                case 'p': /* partial */
                    tinterp_partial(false, rtable, cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0, nz, nx, dx, rtbl[iht], rtbl[iht + 1]);
                    dinterp_partial(false, rtablex, cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0, nz, nx, dx, rtbl[iht], rtbl[iht + 1]);
                    break;

                case 'h': /* hermit */
                    tinterp_hermite(false, rtable, cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0, rtbl[iht], rtbl[iht + 1], rtblx[iht], rtblx[iht + 1]);
                    dinterp_hermite(false, rtablex, cmp ? s + h - iht * rdy - ry0 : h - iht * rdy - ry0, rtbl[iht], rtbl[iht + 1], rtblx[iht], rtblx[iht + 1]);
                    break;
                }
            }

           

#ifdef _OPENMP
#pragma omp parallel for private(iz, ix, t1, t2, ti, tx, ithr)
#endif
            for (i = 0; i < nzx; i++)
            {
#ifdef _OPENMP
                ithr = omp_get_thread_num();
#else
                ithr = 0;
#endif
                iz = i % nz;
                ix = (i - iz) / nz;

                /* aperture (cone angle) */
                if (cmp)
                {
                    if (h >= 0.)
                    {
                        if (atanf((s - x0 - ix * dx) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - s - h) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                    }
                    else
                    {
                        if (atanf((s + h - x0 - ix * dx) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - s) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                    }
                }
                else
                {
                    if (h - s >= 0.)
                    {
                        if (atanf((s - x0 - ix * dx) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - h) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                    }
                    else
                    {
                        if (atanf((h - x0 - ix * dx) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                        if (atanf((x0 + ix * dx - s) / (iz * dz)) * 180. / M_PI > aper)
                            continue;
                    }
                }

                t1 = stable[i];
                t2 = rtable[i];
                ti = t1 + t2 + tau;

                 tx = MAX(fabsf(stablex[i] * ds), fabsf(rtablex[i] * dh));

                kirmig_pick(adj, ti, tx * aal, out[cig ? ih : 0] + i, adj ? trace : traces[ithr]);
            }

            if (!adj)
            {
                for (i = 0; i < nt; i++)
                {
                    trace[i] = 0.0f;
                }
                for (ithr = 0; ithr < nthr; ithr++)
                {
                    for (i = 0; i < nt; i++)
                    {
                        trace[i] += traces[ithr][i];
                    }
                }
                doubint(nt, trace);
                se_fsio_write_float(dat->data->io, trace, nt);
            }
        } /* ih */
    }
    INFO(("FINISH."));

    if (adj)
        se_fsio_write_float(mig->data->io, out[0], nzx * ng);

    se_par_destroy();
    sep_close(dat);
    sep_close(mig);

    if (traces)
        free2float(traces);
    if (out)
        free2float(out);
    if (trace)
        free1float(trace);
    if (stable)
        free1float(stable);
    if (stablex)
        free1float(stablex);
    if (rtable)
        free1float(rtable);
    if (rtablex)
        free1float(rtablex);
    if (stbl)
        free2float(stbl);
    if (stblx)
        free2float(stblx);
    if (rtbl)
        free2float(rtbl);
    if (rtblx)
        free2float(rtblx);
        
    return 0;
}
