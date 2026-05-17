#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void reverse_trace(int nt, int nh, int ns, float ***tr)
{
    INFO(("Reversing traces..."));
    for (int is = 0; is < ns; is++) {
        for (int ih = 0; ih < nh; ih++) {
            for (int it = 0; it < nt / 2; it++) {
                float tmp = tr[is][ih][it];
                tr[is][ih][it] = tr[is][ih][nt - 1 - it];
                tr[is][ih][nt - 1 - it] = tmp;
            }
        }
    }
}

static int next_pow2(int n)
{
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

static void build_filter_F(float tau, float dt, int nsam, float *F)
{
    float *G = (float *)malloc((size_t)nsam * sizeof(float));
    G[0] = 0.0f;
    for (int k = 1; k < nsam; k++) {
        float t = (tau + k * dt) / tau;
        G[k] = sqrtf(t * t - 1.0f);
    }
    for (int k = 0; k < nsam - 1; k++) G[k] = G[k + 1] - G[k];
    for (int k = nsam - 2; k > 0; k--) G[k] = G[k] - G[k - 1];
    for (int k = 0; k < nsam - 1; k++) F[k] = G[k];
    free(G);
}

static void accumulate_pair_fft(
    const float *x, float *out, int nt, int nfft, float dt, int nsam,
    float tau, float weight,
    float *xpad, float *hpad, float *ypad,
    fftwf_complex *X, fftwf_complex *H, fftwf_complex *Y,
    fftwf_plan plan_x, fftwf_plan plan_h, fftwf_plan plan_y,
    float *Fbuf)
{
    if (tau <= 0.0f || !isfinite(tau) || fabsf(weight) < 1e-20f) return;

    const float q = tau / dt;
    const int m = (int)ceilf(q);
    const float delta = (float)m - q;

    memset(xpad, 0, (size_t)nfft * sizeof(float));
    memset(hpad, 0, (size_t)nfft * sizeof(float));
    memcpy(xpad, x, (size_t)nt * sizeof(float));

    build_filter_F(tau, dt, nsam, Fbuf);

    for (int k = 0; k <= nsam - 2; k++) {
        int i1 = m + k;
        int i2 = m + k - 1;
        float fk = Fbuf[k] / dt;
        if (i1 >= 0 && i1 < nfft) hpad[i1] += (1.0f - delta) * fk;
        if (i2 >= 0 && i2 < nfft) hpad[i2] += delta * fk;
    }

    fftwf_execute(plan_x);
    fftwf_execute(plan_h);

    for (int iw = 0; iw <= nfft / 2; iw++) {
        float ar = X[iw][0], ai = X[iw][1];
        float br = H[iw][0], bi = H[iw][1];
        Y[iw][0] = ar * br - ai * bi;
        Y[iw][1] = ar * bi + ai * br;
    }

    fftwf_execute(plan_y);

    for (int it = 0; it < nt; it++) {
        out[it] += weight * (ypad[it] / nfft);
    }

    /* 边界修正，保证与 kirdat_pick 完全一致 */
    const float x0 = x[0];
    for (int k = 0; k <= nsam - 2; k++) {
        int nedge = m + k - 1;
        if (nedge >= 0 && nedge < nt) {
            out[nedge] -= weight * (delta * Fbuf[k] * x0 / dt);
        }
    }
}

int main(int argc, char *argv[])
{
    se_par_init(argc, argv);

    int verb, nt, nh, ns, nsg, nrg, aper, tap;
    int cmp = 1;
    float sdatum = 0.0f, rdatum = 0.0f;
    float length, dt, h0, dh, s0, ds, sg0, dsg, rg0, drg;
    float ***tr_in, ***tr_out, **stable, **rtable;
    sep_t *in, *out, *sgreen, *rgreen, *interm = NULL, *model;
    int should_datum = 0;
    char *in_f = NULL, *out_f = NULL, *sgreen_f = NULL, *rgreen_f = NULL, *interm_f = NULL, *model_f = NULL;

    if(!se_have_par("input_file")) ERROR(("Need input_file=")); else in_f = se_get_par_str("input_file");
    if(!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");
    if(!se_have_par("sgreen_file")) ERROR(("Need sgreen_file=")); else sgreen_f = se_get_par_str("sgreen_file");
    if(!se_have_par("rgreen_file")) ERROR(("Need rgreen_file=")); else rgreen_f = se_get_par_str("rgreen_file");
    if(!se_have_par("model_file")) ERROR(("Need model_file=")); else model_f = se_get_par_str("model_file");

    in = sep_open(in_f, SEP_READ, 0);
    out = sep_open(out_f, SEP_WRITE, 0);
    model = sep_open(model_f, SEP_READ, 0);

    verb = se_have_par("verb") ? se_get_par_int("verb") : 1;
    cmp = se_have_par("cmp") ? se_get_par_int("cmp") : 1;
    if (cmp != 0 && cmp != 1) ERROR(("cmp must be 0 or 1."));

    if (sep_have_hdr_int(model, "should_datum")) should_datum = sep_get_hdr_int(model, "should_datum", 0);
    else ERROR(("Need should_datum in model_file header"));

    aper = se_have_par("aperture") ? se_get_par_int("aperture") : 50;
    tap = se_have_par("taper") ? se_get_par_int("taper") : 10;
    length = se_have_par("length") ? se_get_par_float("length") : 0.025f;

    if(in->headers->ndim < 3) ERROR(("Input must be 3D."));
    nt = in->headers->n[0]; nh = in->headers->n[1]; ns = in->headers->n[2];
    dt = in->headers->d[0]; h0 = in->headers->o[1]; dh = in->headers->d[1];
    s0 = in->headers->o[2]; ds = in->headers->d[2];

    sdatum = (float)(model->headers->o[0] + should_datum * model->headers->d[0]);
    rdatum = sdatum;

    tr_in = alloc3float(nt, nh, ns);
    tr_out = alloc3float(nt, nh, ns);
    se_fsio_read_float(in->data->io, tr_in[0][0], nt * nh * ns);

    reverse_trace(nt, nh, ns, tr_in);

    sgreen = sep_open(sgreen_f, SEP_READ, 0);
    nsg = sgreen->headers->n[0]; sg0 = sgreen->headers->o[0]; dsg = sgreen->headers->d[0];
    stable = alloc2float(nsg, nsg);
    se_fsio_read_float(sgreen->data->io, stable[0], nsg * nsg);
    sep_close(sgreen);

    rgreen = sep_open(rgreen_f, SEP_READ, 0);
    nrg = rgreen->headers->n[0]; rg0 = rgreen->headers->o[0]; drg = rgreen->headers->d[0];
    rtable = alloc2float(nrg, nrg);
    se_fsio_read_float(rgreen->data->io, rtable[0], nrg * nrg);
    sep_close(rgreen);

    if (se_have_par("interm")) {
        interm_f = se_get_par_str("interm");
        interm = sep_open(interm_f, SEP_WRITE, 0);
    }

    int nsam = (int)(length / dt) + 2;
    float max_tau = 0.0f;
    for (int i = 0; i < nrg * nrg; i++) if (rtable[0][i] > max_tau) max_tau = rtable[0][i];
    for (int i = 0; i < nsg * nsg; i++) if (stable[0][i] > max_tau) max_tau = stable[0][i];
    int pmax = (int)ceilf(max_tau / dt) + (nsam - 2);
    int nfft = next_pow2(nt + pmax + 2);

    float *xpad = (float *)fftwf_malloc((size_t)nfft * sizeof(float));
    float *hpad = (float *)fftwf_malloc((size_t)nfft * sizeof(float));
    float *ypad = (float *)fftwf_malloc((size_t)nfft * sizeof(float));
    fftwf_complex *X = (fftwf_complex *)fftwf_malloc((size_t)(nfft / 2 + 1) * sizeof(fftwf_complex));
    fftwf_complex *H = (fftwf_complex *)fftwf_malloc((size_t)(nfft / 2 + 1) * sizeof(fftwf_complex));
    fftwf_complex *Y = (fftwf_complex *)fftwf_malloc((size_t)(nfft / 2 + 1) * sizeof(fftwf_complex));
    float *Fbuf = (float *)malloc((size_t)(nsam - 1) * sizeof(float));

    fftwf_plan plan_x = fftwf_plan_dft_r2c_1d(nfft, xpad, X, FFTW_ESTIMATE);
    fftwf_plan plan_h = fftwf_plan_dft_r2c_1d(nfft, hpad, H, FFTW_ESTIMATE);
    fftwf_plan plan_y = fftwf_plan_dft_c2r_1d(nfft, Y, ypad, FFTW_ESTIMATE);

    for (int is = 0; is < ns; is++) {
        if (verb) INFO(("Processing common-shot gather %d of %d.", is + 1, ns));
        for (int ih = 0; ih < nh; ih++) {
            int c = (int)((((cmp ? (s0 + is * ds) : 0.0f) + h0 + ih * dh - rg0) / drg) + 0.5f);
            int left = (ih - aper < 0) ? 0 : ih - aper;
            int right = (ih + aper > nh - 1) ? nh - 1 : ih + aper;
            for (int ic = left; ic <= right; ic++) {
                int cc = (int)((((cmp ? (s0 + is * ds) : 0.0f) + h0 + ic * dh - rg0) / drg) + 0.5f);
                if (c < 0 || c >= nrg || cc < 0 || cc >= nrg) ERROR(("Receiver table too small."));
                float coef = 1.0f;
                coef *= (ic - left >= tap) ? 1.0f : (float)((ic - left) / tap);
                coef *= (right - ic >= tap) ? 1.0f : (float)((right - ic) / tap);
                float tau = rtable[cc][c];
                float dist = rdatum * rdatum + (ic - ih) * dh * (ic - ih) * dh;
                float w = coef / M_PI * dh * rdatum * tau / dist;
                accumulate_pair_fft(tr_in[is][ic], tr_out[is][ih], nt, nfft, dt, nsam, tau, w,
                                    xpad, hpad, ypad, X, H, Y, plan_x, plan_h, plan_y, Fbuf);
            }
        }
    }

    if (interm) se_fsio_write_float(interm->data->io, tr_out[0][0], nt * nh * ns);
    memset(tr_in[0][0], 0, (size_t)nt * nh * ns * sizeof(float));

    if (cmp == 1) {
        float s = fabsf((ns - 1) * ds), h = fabsf((nh - 1) * dh), dr;
        int jump;
        if (fabsf(ds) >= fabsf(dh)) { dr = fabsf(dh); jump = 1; }
        else { dr = fabsf(ds); jump = (int)(dh / ds + 0.5f); }
        int nr = (int)((s + h) / dr + 1.5f);
        for (int ir = 0; ir < nr; ir++) {
            if (verb) INFO(("Processing common-receiver gather %d of %d.", ir + 1, nr));
            float r = ir * dr + ((ds <= 0.f) ? -1.f : 0.f) * s + ((dh <= 0.f) ? -1.f : 0.f) * h;
            int sleft = (int)((ir * dr + ((ds <= 0.f) ? -1.f : 0.f) * s + ((ds <= 0.f) ? 0.f : -1.f) * h) / ds + 0.5f);
            int sright = (int)((ir * dr + ((ds <= 0.f) ? -1.f : 0.f) * s + ((ds <= 0.f) ? -1.f : 0.f) * h) / ds + 0.5f);
            if (sleft < 0) sleft = 0;
            if (sright > ns - 1) sright = ns - 1;
            int left = (int)((r - sleft * ds) / dh + 0.5f); if (left < 0 || left > nh - 1) sleft++;
            int right = (int)((r - sright * ds) / dh + 0.5f); if (right < 0 || right > nh - 1) sright--;
            for (int is = sleft; is <= sright; is += jump) {
                int c = (int)((s0 + is * ds - sg0) / dsg + 0.5f);
                if (c < 0 || c >= nsg) ERROR(("Source table too small."));
                int ih = (int)((r - is * ds) / dh + 0.5f);
                left = (is - jump * aper < sleft) ? sleft : is - jump * aper;
                right = (is + jump * aper > sright) ? sright : is + jump * aper;
                for (int ic = left; ic <= right; ic += jump) {
                    int cc = (int)((s0 + ic * ds - sg0) / dsg + 0.5f);
                    if (cc < 0 || cc >= nsg) ERROR(("Source table too small."));
                    int hh = (int)((r - ic * ds) / dh + 0.5f);
                    float coef = 1.0f;
                    coef *= (ic - left >= tap) ? 1.0f : (float)((ic - left) / jump / tap);
                    coef *= (right - ic >= tap) ? 1.0f : (float)((right - ic) / jump / tap);
                    float tau = stable[cc][c];
                    float dist = sdatum * sdatum + (ic - is) * ds * (ic - is) * ds;
                    float w = coef / M_PI * ds * sdatum * tau / dist;
                    accumulate_pair_fft(tr_out[ic][hh], tr_in[is][ih], nt, nfft, dt, nsam, tau, w,
                                        xpad, hpad, ypad, X, H, Y, plan_x, plan_h, plan_y, Fbuf);
                }
            }
        }
    } else {
        for (int ih = 0; ih < nh; ih++) {
            if (verb) INFO(("Processing absolute-coordinate receiver gather %d of %d.", ih + 1, nh));
            for (int is = 0; is < ns; is++) {
                int c = (int)((s0 + is * ds - sg0) / dsg + 0.5f);
                if (c < 0 || c >= nsg) ERROR(("Source table too small."));
                int left = (is - aper < 0) ? 0 : is - aper;
                int right = (is + aper > ns - 1) ? ns - 1 : is + aper;
                for (int ic = left; ic <= right; ic++) {
                    int cc = (int)((s0 + ic * ds - sg0) / dsg + 0.5f);
                    if (cc < 0 || cc >= nsg) ERROR(("Source table too small."));
                    float coef = 1.0f;
                    coef *= (ic - left >= tap) ? 1.0f : (float)((ic - left) / tap);
                    coef *= (right - ic >= tap) ? 1.0f : (float)((right - ic) / tap);
                    float tau = stable[cc][c];
                    float dist = sdatum * sdatum + (ic - is) * ds * (ic - is) * ds;
                    float w = coef / M_PI * ds * sdatum * tau / dist;
                    accumulate_pair_fft(tr_out[ic][ih], tr_in[is][ih], nt, nfft, dt, nsam, tau, w,
                                        xpad, hpad, ypad, X, H, Y, plan_x, plan_h, plan_y, Fbuf);
                }
            }
        }
    }

    out->headers->ndim = 3;
    out->headers->n[0] = nt; out->headers->n[1] = nh; out->headers->n[2] = ns;
    out->headers->d[0] = dt; out->headers->d[1] = dh; out->headers->d[2] = ds;
    out->headers->o[0] = 0.; out->headers->o[1] = h0; out->headers->o[2] = s0;

    reverse_trace(nt, nh, ns, tr_in);
    se_fsio_write_float(out->data->io, tr_in[0][0], nt * nh * ns);

    fftwf_destroy_plan(plan_x); fftwf_destroy_plan(plan_h); fftwf_destroy_plan(plan_y);
    fftwf_free(xpad); fftwf_free(hpad); fftwf_free(ypad);
    fftwf_free(X); fftwf_free(H); fftwf_free(Y); free(Fbuf);

    sep_close(in); sep_close(out); sep_close(model); if (interm) sep_close(interm);
    free3float(tr_in); free3float(tr_out); free2float(stable); free2float(rtable);

    INFO(("Done."));
    return 0;
}
