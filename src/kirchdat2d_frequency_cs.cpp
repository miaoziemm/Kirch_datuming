#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <chrono>

#include <fftw3.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define IDX3(is, ih, iw, nh, nw) ((((size_t)(is) * (nh) + (ih)) * (nw)) + (iw))

typedef struct
{
    int nt, nh, ns;
    int nsg, nrg;
    int aper, tap, cmp, verb;
    int nfft, nw, nsam;
    float dt, df, h0, dh, s0, ds, sg0, dsg, rg0, drg;
    float sdatum, rdatum;
    float antialias;
    float **stable, **rtable;
    int n_leaf,p;
} FreqContext;

static void reverse_trace_1d(int nt, float *tr)
{
    for (int it = 0; it < nt / 2; it++)
    {
        float t = tr[it];
        tr[it] = tr[nt - 1 - it];
        tr[nt - 1 - it] = t;
    }
}

static void reverse_traces_2d(int nt, int nh, float **tr)
{
    for (int ih = 0; ih < nh; ih++)
    {
        reverse_trace_1d(nt, tr[ih]);
    }
}

static int next_pow2(int n)
{
    int p = 1;
    while (p < n)
        p <<= 1;
    return p;
}

static int choose_nfft(int nt, int nsam, float dt, int nrg, int nsg, float **rtable, float **stable)
{
    float max_tau = 0.0f;
    for (int i = 0; i < nrg * nrg; i++)
        if (rtable[0][i] > max_tau)
            max_tau = rtable[0][i];
    for (int i = 0; i < nsg * nsg; i++)
        if (stable[0][i] > max_tau)
            max_tau = stable[0][i];
    int pmax = (int)ceilf(max_tau / dt) + (nsam - 2);
    return next_pow2(nt + pmax + 2);
}

static void build_filter_F(float tau, float dt, int nsam, float *F)
{
    float *G = (float *)malloc((size_t)nsam * sizeof(float));
    G[0] = 0.0f;
    for (int k = 1; k < nsam; k++)
    {
        float t = (tau + k * dt) / tau;
        G[k] = sqrtf(t * t - 1.0f);
    }
    for (int k = 0; k < nsam - 1; k++)
        G[k] = G[k + 1] - G[k];
    for (int k = nsam - 2; k > 0; k--)
        G[k] = G[k] - G[k - 1];
    for (int k = 0; k < nsam - 1; k++)
        F[k] = G[k];
    free(G);
}

static void build_H_iw(float tau, const FreqContext *ctx, int iw, float *F, float *hr, float *hi)
{
    if (tau <= 0.0f || !isfinite(tau))
    {
        *hr = 0.0f;
        *hi = 0.0f;
        return;
    }

    build_filter_F(tau, ctx->dt, ctx->nsam, F);

    float q = tau / ctx->dt;
    int m = (int)ceilf(q);
    float delta = (float)m - q;

    /*
     * Original form:
     *   H = sum_k F[k]/dt * [
     *          (1-delta) exp(-i*w*(m+k))
     *        + delta     exp(-i*w*(m+k-1)) ]
     *
     * Factorized form:
     *   H = [ (1-delta) exp(-i*w*m) + delta exp(-i*w*(m-1)) ]
     *       * sum_k F[k]/dt * exp(-i*w*k)
     *
     * This is algebraically the same as the original expression, but avoids
     * two sin/cos evaluations for every filter sample k.
     */
    float omega = 2.0f * (float)M_PI * iw / ctx->nfft;
    float theta = -omega;

    float step_r = cosf(theta);
    float step_i = sinf(theta);

    float z_r = 1.0f;
    float z_i = 0.0f;
    float sum_r = 0.0f;
    float sum_i = 0.0f;

    for (int k = 0; k <= ctx->nsam - 2; k++)
    {
        float fk = F[k] / ctx->dt;

        sum_r += fk * z_r;
        sum_i += fk * z_i;

        float next_r = z_r * step_r - z_i * step_i;
        float next_i = z_r * step_i + z_i * step_r;
        z_r = next_r;
        z_i = next_i;
    }

    float phase_m = -omega * (float)m;
    float phase_m1 = -omega * (float)(m - 1);

    float lin_r = (1.0f - delta) * cosf(phase_m) + delta * cosf(phase_m1);
    float lin_i = (1.0f - delta) * sinf(phase_m) + delta * sinf(phase_m1);

    *hr = lin_r * sum_r - lin_i * sum_i;
    *hi = lin_r * sum_i + lin_i * sum_r;
}

static void accumulate_pair_at_freq(float ar, float ai, float br, float bi, float w, float *or_, float *oi)
{
    *or_ += w * (ar * br - ai * bi);
    *oi += w * (ar * bi + ai * br);
}

static inline float table_slope_first_index(float **table, int n, int row, int col, float dcoord)
{
    if (n <= 1 || dcoord == 0.0f)
        return 0.0f;
    if (row <= 0)
        return (table[1][col] - table[0][col]) / dcoord;
    if (row >= n - 1)
        return (table[n - 1][col] - table[n - 2][col]) / dcoord;
    return 0.5f * (table[row + 1][col] - table[row - 1][col]) / dcoord;
}

static inline float antialias_weight_iw(const FreqContext *ctx, int iw, float dtau)
{
    if (ctx->antialias <= 0.0f)
        return 1.0f;
    if (iw <= 0)
        return 1.0f;

    dtau = fabsf(dtau) * ctx->antialias;
    if (dtau <= 1e-12f || !isfinite(dtau))
        return 1.0f;

    float freq = (float)iw / ((float)ctx->nfft * ctx->dt);
    float fnyq = 0.5f / ctx->dt;
    float fmax = 0.5f / dtau;

    if (fmax >= fnyq)
        return 1.0f;
    if (fmax <= 0.0f)
        return 0.0f;

    /* Smooth transition avoids ringing caused by a hard frequency cut. */
    float fpass = 0.8f * fmax;
    if (freq <= fpass)
        return 1.0f;
    if (freq >= fmax)
        return 0.0f;

    float x = (freq - fpass) / (fmax - fpass);
    return 0.5f * (1.0f + cosf((float)M_PI * x));
}

static inline float taper_weight(int left, int center, int right, int tap, int stride)
{
    if (tap <= 0)
        return 1.0f;
    /*
     * Match kirchdat2d_auto_cs.cpp exactly:
     *   cmp=0/receiver pass: (ic-left)/tap and (right-ic)/tap (integer division)
     *   cmp=1/source pass:   (ic-left)/jump/tap and (right-ic)/jump/tap (integer division)
     * We intentionally keep integer-division semantics here for numerical equivalence.
     */
    float wl = (center - left >= tap) ? 1.0f : (float)(((center - left) / stride) / tap);
    float wr = (right - center >= tap) ? 1.0f : (float)(((right - center) / stride) / tap);
    return wl * wr;
}

static void read_traces_to_freq(se_fsio *io, const FreqContext *ctx, float **shot,
                                fftwf_plan p_f, float *pad, fftwf_complex *spec, fftwf_complex *U)
{
    size_t shot_size = (size_t)ctx->nt * ctx->nh;

    /* Read and FFT one shot at a time to avoid buffering the full time cube. */
    se_fsio_seek(io, 0);
    for (int is = 0; is < ctx->ns; is++)
    {
        se_fsio_read_float(io, shot[0], shot_size);
        reverse_traces_2d(ctx->nt, ctx->nh, shot);
        for (int ih = 0; ih < ctx->nh; ih++)
        {
            memset(pad, 0, (size_t)ctx->nfft * sizeof(float));
            memcpy(pad, shot[ih], (size_t)ctx->nt * sizeof(float));
            fftwf_execute(p_f);
            memcpy(&U[IDX3(is, ih, 0, ctx->nh, ctx->nw)], spec, (size_t)ctx->nw * sizeof(fftwf_complex));
        }
    }
}

static void highpassfilt(fftwf_complex *spec, int nw, float df, float fmin)
{
    if (!spec || nw <= 0 || df <= 0.0f)
        return;
    if (fmin <= 0.0f)
        return;

    float trans = fmaxf(4.0f * df, 0.1f * fmin);
    if (trans >= fmin)
        trans = fmin;
    float f1 = fmin + trans;

    for (int iw = 0; iw < nw; iw++)
    {
        float freq = iw * df;
        float w = 1.0f;
        if (freq <= fmin)
        {
            w = 0.0f;
        }
        else if (freq < f1)
        {
            float x = (freq - fmin) / (f1 - fmin);
            w = 0.5f * (1.0f - cosf((float)M_PI * x));
        }
        if (w <= 0.0f)
        {
            spec[iw][0] = 0.0f;
            spec[iw][1] = 0.0f;
        }
        else if (w < 1.0f)
        {
            spec[iw][0] *= w;
            spec[iw][1] *= w;
        }
    }
}

static void lowpassfilt(fftwf_complex *spec, int nw, float df, float fmax)
{
    if (!spec || nw <= 0 || df <= 0.0f)
        return;
    if (fmax <= 0.0f)
    {
        for (int iw = 0; iw < nw; iw++)
        {
            spec[iw][0] = 0.0f;
            spec[iw][1] = 0.0f;
        }
        return;
    }

    float trans = fmaxf(4.0f * df, 0.1f * fmax);
    if (trans >= fmax)
        trans = fmax;
    float f0 = fmax - trans;

    for (int iw = 0; iw < nw; iw++)
    {
        float freq = iw * df;
        float w = 1.0f;
        if (freq >= fmax)
        {
            w = 0.0f;
        }
        else if (freq > f0)
        {
            float x = (freq - f0) / (fmax - f0);
            w = 0.5f * (1.0f + cosf((float)M_PI * x));
        }
        if (w <= 0.0f)
        {
            spec[iw][0] = 0.0f;
            spec[iw][1] = 0.0f;
        }
        else if (w < 1.0f)
        {
            spec[iw][0] *= w;
            spec[iw][1] *= w;
        }
    }
}

static void bandpassfilt(fftwf_complex *spec, int nw, float df, float fmin, float fmax)
{
    if (!spec || nw <= 0 || df <= 0.0f)
        return;

    if (fmax <= 0.0f || fmax < fmin)
    {
        for (int iw = 0; iw < nw; iw++)
        {
            spec[iw][0] = 0.0f;
            spec[iw][1] = 0.0f;
        }
        return;
    }

    if (fmin <= 0.0f)
    {
        lowpassfilt(spec, nw, df, fmax);
        return;
    }

    float lo_trans = fmaxf(4.0f * df, 0.1f * fmin);
    if (lo_trans >= fmin)
        lo_trans = fmin;
    float hi_trans = fmaxf(4.0f * df, 0.1f * fmax);
    if (hi_trans >= fmax)
        hi_trans = fmax;

    float flo0 = fmin;
    float flo1 = fmin + lo_trans;
    float fhi0 = fmax - hi_trans;
    float fhi1 = fmax;

    for (int iw = 0; iw < nw; iw++)
    {
        float freq = iw * df;
        float wl = 1.0f;
        float wh = 1.0f;

        if (freq <= flo0 || freq >= fhi1)
        {
            wl = 0.0f;
            wh = 0.0f;
        }
        else
        {
            if (freq < flo1)
            {
                float x = (freq - flo0) / (flo1 - flo0);
                wl = 0.5f * (1.0f - cosf((float)M_PI * x));
            }
            if (freq > fhi0)
            {
                float x = (freq - fhi0) / (fhi1 - fhi0);
                wh = 0.5f * (1.0f + cosf((float)M_PI * x));
            }
        }

        float w = wl * wh;
        if (w <= 0.0f)
        {
            spec[iw][0] = 0.0f;
            spec[iw][1] = 0.0f;
        }
        else if (w < 1.0f)
        {
            spec[iw][0] *= w;
            spec[iw][1] *= w;
        }
    }
}

static void write_traces_from_freq(fftwf_complex *U, const FreqContext *ctx, fftwf_plan p_b,
                                   fftwf_complex *spec, float *pad, se_fsio *io, int reverse)
{
    for (int is = 0; is < ctx->ns; is++)
    {
        for (int ih = 0; ih < ctx->nh; ih++)
        {
            memcpy(spec, &U[IDX3(is, ih, 0, ctx->nh, ctx->nw)], (size_t)ctx->nw * sizeof(fftwf_complex));
            // bandpassfilt(spec, ctx->nw, ctx->df, 10, 80);
            fftwf_execute(p_b);
            for (int it = 0; it < ctx->nt; it++)
            {
                pad[it] /= ctx->nfft;
            }
            if (reverse)
            {
                reverse_trace_1d(ctx->nt, pad);
            }
            se_fsio_write_float(io, pad, (size_t)ctx->nt);
        }
    }
}

static void apply_receiver_operator_freq_iw(const FreqContext *ctx, fftwf_complex *Uin, fftwf_complex *Utmp, int iw)
{
    float *F = (float *)malloc((size_t)(ctx->nsam - 1) * sizeof(float));
    float hr, hi;
    float *hr_cache = (float *)malloc((size_t)ctx->nrg * ctx->nrg * sizeof(float));
    float *hi_cache = (float *)malloc((size_t)ctx->nrg * ctx->nrg * sizeof(float));
    unsigned char *cached = (unsigned char *)calloc((size_t)ctx->nrg * ctx->nrg, sizeof(unsigned char));
    for (int is = 0; is < ctx->ns; is++)
    {
        for (int ih = 0; ih < ctx->nh; ih++)
        {
            int c = (int)((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f) + ctx->h0 + ih * ctx->dh - ctx->rg0) / ctx->drg) + 0.5f);
            int left = (ih - ctx->aper < 0) ? 0 : ih - ctx->aper;
            int right = (ih + ctx->aper > ctx->nh - 1) ? ctx->nh - 1 : ih + ctx->aper;
            for (int ic = left; ic <= right; ic++)
            {
                int cc = (int)((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f) + ctx->h0 + ic * ctx->dh - ctx->rg0) / ctx->drg) + 0.5f);
                if (c < 0 || c >= ctx->nrg || cc < 0 || cc >= ctx->nrg)
                    ERROR(("Receiver table too small."));
                float coef = taper_weight(left, ic, right, ctx->tap, 1);
                float tau = ctx->rtable[cc][c];
                float dist = ctx->rdatum * ctx->rdatum + (ic - ih) * (ic - ih) * ctx->dh * ctx->dh;
                float w = coef / M_PI * ctx->dh * ctx->rdatum * tau / dist;

                /* Frequency-domain operator anti-aliasing on receiver summation.
                 * The local phase increment is approximated by dT/dr * dh.
                 */
                float slope = table_slope_first_index(ctx->rtable, ctx->nrg, cc, c, ctx->drg);
                float aa = antialias_weight_iw(ctx, iw, slope * fabsf(ctx->dh));
                w *= aa;

                if (fabsf(w) < 1e-20f)
                    continue;
                size_t key = (size_t)cc * ctx->nrg + c;
                if (!cached[key])
                {
                    build_H_iw(tau, ctx, iw, F, &hr_cache[key], &hi_cache[key]);
                    cached[key] = 1;
                }
                hr = hr_cache[key];
                hi = hi_cache[key];
                size_t ia = IDX3(is, ic, 0, ctx->nh, ctx->nw);
                size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);
                float ar = Uin[ia + iw][0], ai = Uin[ia + iw][1];
                float br = hr, bi = hi;
                accumulate_pair_at_freq(ar, ai, br, bi, w, &Utmp[ib + iw][0], &Utmp[ib + iw][1]);
            }
        }
    }
    free(hr_cache);
    free(hi_cache);
    free(cached);
    free(F);
}

static void apply_source_operator_freq_iw(const FreqContext *ctx, fftwf_complex *Utmp, fftwf_complex *Uout, int iw)
{
    float *F = (float *)malloc((size_t)(ctx->nsam - 1) * sizeof(float));
    float hr, hi;
    float *hr_cache = (float *)malloc((size_t)ctx->nsg * ctx->nsg * sizeof(float));
    float *hi_cache = (float *)malloc((size_t)ctx->nsg * ctx->nsg * sizeof(float));
    unsigned char *cached = (unsigned char *)calloc((size_t)ctx->nsg * ctx->nsg, sizeof(unsigned char));
    if (ctx->cmp == 1)
    {
        float s = fabsf((ctx->ns - 1) * ctx->ds), h = fabsf((ctx->nh - 1) * ctx->dh), dr;
        int jump = (fabsf(ctx->ds) >= fabsf(ctx->dh)) ? 1 : (int)(ctx->dh / ctx->ds + 0.5f);
        dr = (fabsf(ctx->ds) >= fabsf(ctx->dh)) ? fabsf(ctx->dh) : fabsf(ctx->ds);
        int nr = (int)((s + h) / dr + 1.5f);
        for (int ir = 0; ir < nr; ir++)
        {
            float r = ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->dh <= 0.f) ? -1.f : 0.f) * h;
            int sleft = (int)((ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->ds <= 0.f) ? 0.f : -1.f) * h) / ctx->ds + 0.5f);
            int sright = (int)((ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->ds <= 0.f) ? -1.f : 0.f) * h) / ctx->ds + 0.5f);
            if (sleft < 0)
                sleft = 0;
            if (sright > ctx->ns - 1)
                sright = ctx->ns - 1;
            int left = (int)((r - sleft * ctx->ds) / ctx->dh + 0.5f);
            if (left < 0 || left > ctx->nh - 1)
                sleft++;
            int right = (int)((r - sright * ctx->ds) / ctx->dh + 0.5f);
            if (right < 0 || right > ctx->nh - 1)
                sright--;
            for (int is = sleft; is <= sright; is += jump)
            {
                int c = (int)((ctx->s0 + is * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
                if (c < 0 || c >= ctx->nsg)
                    ERROR(("Source table too small."));
                int ih = (int)((r - is * ctx->ds) / ctx->dh + 0.5f);
                left = (is - jump * ctx->aper < sleft) ? sleft : is - jump * ctx->aper;
                right = (is + jump * ctx->aper > sright) ? sright : is + jump * ctx->aper;
                for (int ic = left; ic <= right; ic += jump)
                {
                    int cc = (int)((ctx->s0 + ic * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
                    if (cc < 0 || cc >= ctx->nsg)
                        ERROR(("Source table too small."));
                    int hh = (int)((r - ic * ctx->ds) / ctx->dh + 0.5f);
                    float coef = taper_weight(left, ic, right, ctx->tap, jump);
                    float tau = ctx->stable[cc][c];
                    float dist = ctx->sdatum * ctx->sdatum + (ic - is) * (ic - is) * ctx->ds * ctx->ds;
                    float w = coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;

                    /* Frequency-domain operator anti-aliasing on source summation.
                     * In cmp mode the source loop may skip by jump, so the sampled
                     * source interval is jump * ds.
                     */
                    float slope = table_slope_first_index(ctx->stable, ctx->nsg, cc, c, ctx->dsg);
                    float aa = antialias_weight_iw(ctx, iw, slope * fabsf((float)jump * ctx->ds));
                    w *= aa;

                    if (fabsf(w) < 1e-20f)
                        continue;
                    size_t key = (size_t)cc * ctx->nsg + c;
                    if (!cached[key])
                    {
                        build_H_iw(tau, ctx, iw, F, &hr_cache[key], &hi_cache[key]);
                        cached[key] = 1;
                    }
                    hr = hr_cache[key];
                    hi = hi_cache[key];
                    size_t ia = IDX3(ic, hh, 0, ctx->nh, ctx->nw);
                    size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);
                    float ar = Utmp[ia + iw][0], ai = Utmp[ia + iw][1];
                    float br = hr, bi = hi;
                    accumulate_pair_at_freq(ar, ai, br, bi, w, &Uout[ib + iw][0], &Uout[ib + iw][1]);
                }
            }
        }
    }
    else
    {
        for (int ih = 0; ih < ctx->nh; ih++)
        {
            for (int is = 0; is < ctx->ns; is++)
            {
                int c = (int)((ctx->s0 + is * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
                if (c < 0 || c >= ctx->nsg)
                    ERROR(("Source table too small."));
                int left = (is - ctx->aper < 0) ? 0 : is - ctx->aper;
                int right = (is + ctx->aper > ctx->ns - 1) ? ctx->ns - 1 : is + ctx->aper;
                for (int ic = left; ic <= right; ic++)
                {
                    int cc = (int)((ctx->s0 + ic * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
                    if (cc < 0 || cc >= ctx->nsg)
                        ERROR(("Source table too small."));
                    float coef = taper_weight(left, ic, right, ctx->tap, 1);
                    float tau = ctx->stable[cc][c];
                    float dist = ctx->sdatum * ctx->sdatum + (ic - is) * (ic - is) * ctx->ds * ctx->ds;
                    float w = coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;

                    /* Frequency-domain operator anti-aliasing on source summation. */
                    float slope = table_slope_first_index(ctx->stable, ctx->nsg, cc, c, ctx->dsg);
                    float aa = antialias_weight_iw(ctx, iw, slope * fabsf(ctx->ds));
                    w *= aa;

                    if (fabsf(w) < 1e-20f)
                        continue;
                    size_t key = (size_t)cc * ctx->nsg + c;
                    if (!cached[key])
                    {
                        build_H_iw(tau, ctx, iw, F, &hr_cache[key], &hi_cache[key]);
                        cached[key] = 1;
                    }
                    hr = hr_cache[key];
                    hi = hi_cache[key];
                    size_t ia = IDX3(ic, ih, 0, ctx->nh, ctx->nw);
                    size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);
                    float ar = Utmp[ia + iw][0], ai = Utmp[ia + iw][1];
                    float br = hr, bi = hi;
                    accumulate_pair_at_freq(ar, ai, br, bi, w, &Uout[ib + iw][0], &Uout[ib + iw][1]);
                }
            }
        }
    }
    free(hr_cache);
    free(hi_cache);
    free(cached);
    free(F);
}

static void apply_receiver_operator_freq_iw_bf(const FreqContext *ctx,
                                               fftwf_complex *Uin,
                                               fftwf_complex *Utmp,
                                               int iw)
{
    float *F = (float *)malloc((size_t)(ctx->nsam - 1) * sizeof(float));

    float *hr_cache =
        (float *)malloc((size_t)ctx->nrg * ctx->nrg * sizeof(float));

    float *hi_cache =
        (float *)malloc((size_t)ctx->nrg * ctx->nrg * sizeof(float));

    unsigned char *cached =
        (unsigned char *)calloc((size_t)ctx->nrg * ctx->nrg,
                                sizeof(unsigned char));

    /*
     * Traveltime matrix:
     *
     *     tau_mat[ih][ic] = tau from input receiver ic to output receiver ih
     *
     * This is the phase matrix used by the butterfly algorithm.
     */
    float **tau_mat = alloc2float(ctx->nh, ctx->nh);

    /*
     * Amplitude matrix after phase removal:
     *
     *     Amp_mat[ih,ic]
     *       = w(ih,ic) * H(tau,iw) * exp(+i * omega * tau)
     *
     * Therefore, the full kernel is reconstructed by:
     *
     *     K[ih,ic]
     *       = Amp_mat[ih,ic] * exp(-i * omega * tau_mat[ih][ic]).
     *
     * Amp_mat is stored as a 1D complex array:
     *
     *     imat = ih * nh + ic.
     */
    fftwf_complex *Amp_mat =
        (fftwf_complex *)fftwf_malloc((size_t)ctx->nh * (size_t)ctx->nh * sizeof(fftwf_complex));

    /*
     * Input and output vectors for one fixed source is
     * and one fixed frequency iw:
     *
     *     Uin_is[ic]  = Uin[is, ic, iw]
     *     Uout_is[ih] = receiver-side extrapolated result
     */
    fftwf_complex *Uin_is =
        (fftwf_complex *)fftwf_malloc((size_t)ctx->nh * sizeof(fftwf_complex));

    fftwf_complex *Uout_is =
        (fftwf_complex *)fftwf_malloc((size_t)ctx->nh * sizeof(fftwf_complex));

    if (F == NULL || hr_cache == NULL || hi_cache == NULL ||
        cached == NULL || tau_mat == NULL || Amp_mat == NULL ||
        Uin_is == NULL || Uout_is == NULL)
    {
        ERROR(("Out of memory."));
    }

    /*
     * Physical angular frequency.
     *
     * build_H_iw internally uses the digital angular frequency
     *
     *     omega_d = 2*pi*iw/nfft.
     *
     * Since tau is measured in seconds, the corresponding physical
     * angular frequency is:
     *
     *     omega = omega_d / dt = 2*pi*iw/(nfft*dt).
     */
    float omega =
        2.0f * (float)M_PI * (float)iw / ((float)ctx->nfft * ctx->dt);

    for (int is = 0; is < ctx->ns; is++)
    {

        /*
         * Initialize tau_mat, Amp_mat, and Uout_is for the current source.
         *
         * Values outside the aperture remain zero.
         */
        for (int ih = 0; ih < ctx->nh; ih++)
        {

            Uout_is[ih][0] = 0.0f;
            Uout_is[ih][1] = 0.0f;

            for (int ic = 0; ic < ctx->nh; ic++)
            {

                size_t imat =
                    (size_t)ih * (size_t)ctx->nh + (size_t)ic;

                tau_mat[ih][ic] = 0.0f;

                Amp_mat[imat][0] = 0.0f;
                Amp_mat[imat][1] = 0.0f;
            }
        }

        /*
         * Extract input vector at the current source and frequency:
         *
         *     Uin_is[ic] = Uin[is, ic, iw].
         */
        for (int ic = 0; ic < ctx->nh; ic++)
        {

            size_t ia = IDX3(is, ic, 0, ctx->nh, ctx->nw);

            Uin_is[ic][0] = Uin[ia + iw][0];
            Uin_is[ic][1] = Uin[ia + iw][1];
        }

        /*
         * Build tau_mat and Amp_mat.
         *
         * For each output receiver ih, only input receivers inside
         * [left, right] are used. Outside this aperture, Amp_mat remains zero.
         */
        for (int ih = 0; ih < ctx->nh; ih++)
        {

            float s_abs =
                ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f;

            int c = (int)(((s_abs + ctx->h0 + ih * ctx->dh - ctx->rg0) / ctx->drg) + 0.5f);

            int left = (ih - ctx->aper < 0)
                           ? 0
                           : ih - ctx->aper;

            int right = (ih + ctx->aper > ctx->nh - 1)
                            ? ctx->nh - 1
                            : ih + ctx->aper;

            for (int ic = left; ic <= right; ic++)
            {

                int cc = (int)(((s_abs + ctx->h0 + ic * ctx->dh - ctx->rg0) / ctx->drg) + 0.5f);

                if (c < 0 || c >= ctx->nrg ||
                    cc < 0 || cc >= ctx->nrg)
                {
                    ERROR(("Receiver table too small."));
                }

                float coef = taper_weight(left, ic, right, ctx->tap, 1);

                float tau = ctx->rtable[cc][c];

                float dist =
                    ctx->rdatum * ctx->rdatum + (ic - ih) * (ic - ih) * ctx->dh * ctx->dh;

                float w =
                    coef / M_PI * ctx->dh * ctx->rdatum * tau / dist;

                /*
                 * Frequency-domain operator anti-aliasing on receiver summation.
                 * The local phase increment is approximated by dT/dr * dh.
                 */
                float slope =
                    table_slope_first_index(ctx->rtable,
                                            ctx->nrg,
                                            cc,
                                            c,
                                            ctx->drg);

                float aa =
                    antialias_weight_iw(ctx,
                                        iw,
                                        slope * fabsf(ctx->dh));

                w *= aa;

                if (fabsf(w) < 1e-20f)
                    continue;

                /*
                 * Build H(tau, iw) using the same cache strategy as
                 * the original direct-summation code.
                 */
                size_t key =
                    (size_t)cc * (size_t)ctx->nrg + (size_t)c;

                if (!cached[key])
                {
                    build_H_iw(tau,
                               ctx,
                               iw,
                               F,
                               &hr_cache[key],
                               &hi_cache[key]);

                    cached[key] = 1;
                }

                float hr = hr_cache[key];
                float hi = hi_cache[key];

                /*
                 * Phase-amplitude separation:
                 *
                 *     H(tau,iw) = B(tau,iw) * exp(-i * omega * tau)
                 *
                 * Therefore:
                 *
                 *     B(tau,iw) = H(tau,iw) * exp(+i * omega * tau)
                 *
                 * and
                 *
                 *     Amp_mat = w * B
                 *             = w * H(tau,iw) * exp(+i * omega * tau).
                 */
                float phase = omega * tau;
                float cosp = cosf(phase);
                float sinp = sinf(phase);

                /*
                 * H * exp(+i*phase):
                 *
                 *     (hr + i hi) * (cos + i sin)
                 *   = (hr*cos - hi*sin) + i(hr*sin + hi*cos)
                 */
                float amp_r = w * (hr * cosp - hi * sinp);
                float amp_i = w * (hr * sinp + hi * cosp);

                size_t imat =
                    (size_t)ih * (size_t)ctx->nh + (size_t)ic;

                tau_mat[ih][ic] = tau;

                Amp_mat[imat][0] = amp_r;
                Amp_mat[imat][1] = amp_i;
            }
        }

        // /*
        //  * Receiver-side extrapolation in phase-amplitude-separated form:
        //  *
        //  *     Uout_is[ih]
        //  *       = sum_ic Amp_mat[ih,ic]
        //  *                * exp(-i * omega * tau_mat[ih][ic])
        //  *                * Uin_is[ic].
        //  *
        //  * This direct summation is only for verification.
        //  * Later it can be replaced by a butterfly algorithm using
        //  * tau_mat as the phase matrix and Amp_mat as the amplitude matrix.
        //  */
        // for (int ih = 0; ih < ctx->nh; ih++) {
        //     for (int ic = 0; ic < ctx->nh; ic++) {

        //         size_t imat =
        //             (size_t)ih * (size_t)ctx->nh + (size_t)ic;

        //         float amp_r = Amp_mat[imat][0];
        //         float amp_i = Amp_mat[imat][1];

        //         if (fabsf(amp_r) < 1e-20f &&
        //             fabsf(amp_i) < 1e-20f) {
        //             continue;
        //         }

        //         float tau = tau_mat[ih][ic];

        //         /*
        //          * exp(-i * omega * tau)
        //          *   = cos(omega*tau) - i sin(omega*tau).
        //          */
        //         float phase = omega * tau;
        //         float cosp = cosf(phase);
        //         float sinp = sinf(phase);

        //         /*
        //          * Reconstruct the full kernel:
        //          *
        //          *     K = Amp * exp(-i*phase)
        //          *
        //          * If Amp = amp_r + i amp_i, then
        //          *
        //          *     K_r = amp_r*cos + amp_i*sin
        //          *     K_i = amp_i*cos - amp_r*sin
        //          */
        //         float kr = amp_r * cosp + amp_i * sinp;
        //         float ki = amp_i * cosp - amp_r * sinp;

        //         float ar = Uin_is[ic][0];
        //         float ai = Uin_is[ic][1];

        //         /*
        //          * Uout_is[ih] += K * Uin_is[ic].
        //          */
        //         Uout_is[ih][0] += kr * ar - ki * ai;
        //         Uout_is[ih][1] += kr * ai + ki * ar;
        //     }
        // }

        butterfly_apply_1d_phase_amp(ctx->ns, tau_mat, Amp_mat, Uin_is, Uout_is, omega, ctx->p, ctx->n_leaf);

        /*
         * Butterfly replacement interface:
         *
         *     butterfly_receiver_apply(ctx,
         *                              iw,
         *                              tau_mat,
         *                              Amp_mat,
         *                              Uin_is,
         *                              Uout_is);
         *
         * Inside the butterfly algorithm, the kernel should be evaluated as:
         *
         *     K(ih,ic)
         *       = Amp_mat[ih,ic]
         *         * exp(-i * omega * tau_mat[ih][ic]).
         */

        /*
         * Write current-source output vector back to Utmp.
         */
        for (int ih = 0; ih < ctx->nh; ih++)
        {

            size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);

            Utmp[ib + iw][0] += Uout_is[ih][0];
            Utmp[ib + iw][1] += Uout_is[ih][1];
        }
    }

    fftwf_free(Amp_mat);
    fftwf_free(Uin_is);
    fftwf_free(Uout_is);

    free2float(tau_mat);

    free(hr_cache);
    free(hi_cache);
    free(cached);
    free(F);
}

static void apply_source_operator_freq_iw_bf(const FreqContext *ctx,
                                             fftwf_complex *Utmp,
                                             fftwf_complex *Uout,
                                             int iw)
{
    float *F = (float *)malloc((size_t)(ctx->nsam - 1) * sizeof(float));

    float *hr_cache =
        (float *)malloc((size_t)ctx->nsg * ctx->nsg * sizeof(float));

    float *hi_cache =
        (float *)malloc((size_t)ctx->nsg * ctx->nsg * sizeof(float));

    unsigned char *cached =
        (unsigned char *)calloc((size_t)ctx->nsg * ctx->nsg,
                                sizeof(unsigned char));

    /*
     * Source-side traveltime matrix:
     *
     *     tau_mat[is][ic] = tau from input source ic to output source is
     *
     * This is the phase matrix for the butterfly algorithm.
     */
    float **tau_mat = alloc2float(ctx->ns, ctx->ns);

    /*
     * Source-side amplitude matrix after phase removal:
     *
     *     Amp_mat[is,ic]
     *       = w(is,ic) * H(tau,iw) * exp(+i * omega * tau)
     *
     * The full kernel is reconstructed as:
     *
     *     K(is,ic)
     *       = Amp_mat[is,ic] * exp(-i * omega * tau_mat[is][ic]).
     *
     * Stored as a 1D complex array:
     *
     *     imat = is * ns + ic.
     */
    fftwf_complex *Amp_mat =
        (fftwf_complex *)fftwf_malloc((size_t)ctx->ns * (size_t)ctx->ns * sizeof(fftwf_complex));

    /*
     * Input and output vectors for one source-side matrix-vector product:
     *
     *     Uin_is[ic]  = input wavefield at input source index ic
     *     Uout_is[is] = extrapolated output at output source index is
     */
    fftwf_complex *Uin_is =
        (fftwf_complex *)fftwf_malloc((size_t)ctx->ns * sizeof(fftwf_complex));

    fftwf_complex *Uout_is =
        (fftwf_complex *)fftwf_malloc((size_t)ctx->ns * sizeof(fftwf_complex));

    if (F == NULL || hr_cache == NULL || hi_cache == NULL ||
        cached == NULL || tau_mat == NULL || Amp_mat == NULL ||
        Uin_is == NULL || Uout_is == NULL)
    {
        ERROR(("Out of memory."));
    }

    /*
     * Physical angular frequency:
     *
     *     omega = 2*pi*f = 2*pi*iw/(nfft*dt).
     */
    float omega =
        2.0f * (float)M_PI * (float)iw / ((float)ctx->nfft * ctx->dt);

    if (ctx->cmp == 1)
    {

        float s = fabsf((ctx->ns - 1) * ctx->ds);
        float h = fabsf((ctx->nh - 1) * ctx->dh);
        float dr;

        int jump =
            (fabsf(ctx->ds) >= fabsf(ctx->dh))
                ? 1
                : (int)(ctx->dh / ctx->ds + 0.5f);

        dr =
            (fabsf(ctx->ds) >= fabsf(ctx->dh))
                ? fabsf(ctx->dh)
                : fabsf(ctx->ds);

        int nr = (int)((s + h) / dr + 1.5f);

        /*
         * In cmp mode, source-side extrapolation is performed along
         * constant receiver coordinate r.
         */
        for (int ir = 0; ir < nr; ir++)
        {

            float r =
                ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->dh <= 0.f) ? -1.f : 0.f) * h;

            int sleft =
                (int)((ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->ds <= 0.f) ? 0.f : -1.f) * h) / ctx->ds + 0.5f);

            int sright =
                (int)((ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->ds <= 0.f) ? -1.f : 0.f) * h) / ctx->ds + 0.5f);

            if (sleft < 0)
                sleft = 0;
            if (sright > ctx->ns - 1)
                sright = ctx->ns - 1;

            int left =
                (int)((r - sleft * ctx->ds) / ctx->dh + 0.5f);

            if (left < 0 || left > ctx->nh - 1)
                sleft++;

            int right =
                (int)((r - sright * ctx->ds) / ctx->dh + 0.5f);

            if (right < 0 || right > ctx->nh - 1)
                sright--;

            /*
             * Initialize tau_mat, Amp_mat, Uin_is, and Uout_is
             * for the current constant-r gather.
             */
            for (int is = 0; is < ctx->ns; is++)
            {

                Uin_is[is][0] = 0.0f;
                Uin_is[is][1] = 0.0f;

                Uout_is[is][0] = 0.0f;
                Uout_is[is][1] = 0.0f;

                for (int ic = 0; ic < ctx->ns; ic++)
                {

                    size_t imat =
                        (size_t)is * (size_t)ctx->ns + (size_t)ic;

                    tau_mat[is][ic] = 0.0f;

                    Amp_mat[imat][0] = 0.0f;
                    Amp_mat[imat][1] = 0.0f;
                }
            }

            /*
             * Extract input vector along the current constant-r line:
             *
             *     Uin_is[ic] = Utmp[ic, hh, iw],
             *
             * where:
             *
             *     hh = round((r - ic * ds) / dh).
             */
            for (int ic = sleft; ic <= sright; ic += jump)
            {

                int hh =
                    (int)((r - ic * ctx->ds) / ctx->dh + 0.5f);

                if (hh < 0 || hh >= ctx->nh)
                    continue;

                size_t ia = IDX3(ic, hh, 0, ctx->nh, ctx->nw);

                Uin_is[ic][0] = Utmp[ia + iw][0];
                Uin_is[ic][1] = Utmp[ia + iw][1];
            }

            /*
             * Build tau_mat and Amp_mat for the current constant-r line.
             */
            for (int is = sleft; is <= sright; is += jump)
            {

                int c =
                    (int)((ctx->s0 + is * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);

                if (c < 0 || c >= ctx->nsg)
                {
                    ERROR(("Source table too small."));
                }

                int ih =
                    (int)((r - is * ctx->ds) / ctx->dh + 0.5f);

                if (ih < 0 || ih >= ctx->nh)
                    continue;

                int ileft =
                    (is - jump * ctx->aper < sleft)
                        ? sleft
                        : is - jump * ctx->aper;

                int iright =
                    (is + jump * ctx->aper > sright)
                        ? sright
                        : is + jump * ctx->aper;

                for (int ic = ileft; ic <= iright; ic += jump)
                {

                    int cc =
                        (int)((ctx->s0 + ic * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);

                    if (cc < 0 || cc >= ctx->nsg)
                    {
                        ERROR(("Source table too small."));
                    }

                    int hh =
                        (int)((r - ic * ctx->ds) / ctx->dh + 0.5f);

                    if (hh < 0 || hh >= ctx->nh)
                        continue;

                    float coef =
                        taper_weight(ileft, ic, iright, ctx->tap, jump);

                    float tau = ctx->stable[cc][c];

                    float dist =
                        ctx->sdatum * ctx->sdatum + (ic - is) * (ic - is) * ctx->ds * ctx->ds;

                    float w =
                        coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;

                    /*
                     * Frequency-domain operator anti-aliasing on source summation.
                     * In cmp mode the source loop may skip by jump, so the sampled
                     * source interval is jump * ds.
                     */
                    float slope =
                        table_slope_first_index(ctx->stable,
                                                ctx->nsg,
                                                cc,
                                                c,
                                                ctx->dsg);

                    float aa =
                        antialias_weight_iw(ctx,
                                            iw,
                                            slope * fabsf((float)jump * ctx->ds));

                    w *= aa;

                    if (fabsf(w) < 1e-20f)
                        continue;

                    /*
                     * Build H(tau, iw) using the same cache strategy
                     * as the original direct-summation code.
                     */
                    size_t key =
                        (size_t)cc * (size_t)ctx->nsg + (size_t)c;

                    if (!cached[key])
                    {
                        build_H_iw(tau,
                                   ctx,
                                   iw,
                                   F,
                                   &hr_cache[key],
                                   &hi_cache[key]);

                        cached[key] = 1;
                    }

                    float hr = hr_cache[key];
                    float hi = hi_cache[key];

                    /*
                     * Phase-amplitude separation:
                     *
                     *     Amp = w * H(tau,iw) * exp(+i * omega * tau).
                     */
                    float phase = omega * tau;
                    float cosp = cosf(phase);
                    float sinp = sinf(phase);

                    /*
                     * H * exp(+i*phase):
                     *
                     *     (hr + i hi) * (cos + i sin)
                     *   = (hr*cos - hi*sin) + i(hr*sin + hi*cos).
                     */
                    float amp_r = w * (hr * cosp - hi * sinp);
                    float amp_i = w * (hr * sinp + hi * cosp);

                    size_t imat =
                        (size_t)is * (size_t)ctx->ns + (size_t)ic;

                    tau_mat[is][ic] = tau;

                    Amp_mat[imat][0] = amp_r;
                    Amp_mat[imat][1] = amp_i;
                }
            }

            /*
             * Direct summation for verification:
             *
             *     Uout_is[is]
             *       = sum_ic Amp_mat[is,ic]
             *                * exp(-i * omega * tau_mat[is][ic])
             *                * Uin_is[ic].
             *
             * Later this block can be replaced by a butterfly algorithm.
             */
            // for (int is = sleft; is <= sright; is += jump) {
            //     for (int ic = sleft; ic <= sright; ic += jump) {

            //         size_t imat =
            //             (size_t)is * (size_t)ctx->ns + (size_t)ic;

            //         float amp_r = Amp_mat[imat][0];
            //         float amp_i = Amp_mat[imat][1];

            //         if (fabsf(amp_r) < 1e-20f &&
            //             fabsf(amp_i) < 1e-20f) {
            //             continue;
            //         }

            //         float tau = tau_mat[is][ic];

            //         float phase = omega * tau;
            //         float cosp = cosf(phase);
            //         float sinp = sinf(phase);

            //         /*
            //          * K = Amp * exp(-i*phase).
            //          *
            //          * If Amp = amp_r + i amp_i, then:
            //          *
            //          *     K_r = amp_r*cos + amp_i*sin
            //          *     K_i = amp_i*cos - amp_r*sin.
            //          */
            //         float kr = amp_r * cosp + amp_i * sinp;
            //         float ki = amp_i * cosp - amp_r * sinp;

            //         float ar = Uin_is[ic][0];
            //         float ai = Uin_is[ic][1];

            //         Uout_is[is][0] += kr * ar - ki * ai;
            //         Uout_is[is][1] += kr * ai + ki * ar;
            //     }
            // }

            butterfly_apply_1d_phase_amp(ctx->ns,
                                         tau_mat,
                                         Amp_mat,
                                         Uin_is,
                                         Uout_is,
                                         omega, ctx->p, ctx->n_leaf);

            /*
             * Butterfly replacement interface:
             *
             *     butterfly_source_apply(ctx,
             *                            iw,
             *                            tau_mat,
             *                            Amp_mat,
             *                            Uin_is,
             *                            Uout_is);
             *
             * Inside the butterfly algorithm, the kernel should be:
             *
             *     K(is,ic)
             *       = Amp_mat[is,ic]
             *         * exp(-i * omega * tau_mat[is][ic]).
             */

            /*
             * Write the current constant-r output vector back to Uout.
             */
            for (int is = sleft; is <= sright; is += jump)
            {

                int ih =
                    (int)((r - is * ctx->ds) / ctx->dh + 0.5f);

                if (ih < 0 || ih >= ctx->nh)
                    continue;

                size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);

                Uout[ib + iw][0] += Uout_is[is][0];
                Uout[ib + iw][1] += Uout_is[is][1];
            }
        }
    }
    else
    {

        /*
         * In non-cmp mode, source-side extrapolation is performed
         * for each fixed receiver/offset index ih.
         */
        for (int ih = 0; ih < ctx->nh; ih++)
        {

            /*
             * Initialize tau_mat, Amp_mat, Uin_is, and Uout_is
             * for the current ih.
             */
            for (int is = 0; is < ctx->ns; is++)
            {

                Uin_is[is][0] = 0.0f;
                Uin_is[is][1] = 0.0f;

                Uout_is[is][0] = 0.0f;
                Uout_is[is][1] = 0.0f;

                for (int ic = 0; ic < ctx->ns; ic++)
                {

                    size_t imat =
                        (size_t)is * (size_t)ctx->ns + (size_t)ic;

                    tau_mat[is][ic] = 0.0f;

                    Amp_mat[imat][0] = 0.0f;
                    Amp_mat[imat][1] = 0.0f;
                }
            }

            /*
             * Extract input vector for current ih:
             *
             *     Uin_is[ic] = Utmp[ic, ih, iw].
             */
            for (int ic = 0; ic < ctx->ns; ic++)
            {

                size_t ia = IDX3(ic, ih, 0, ctx->nh, ctx->nw);

                Uin_is[ic][0] = Utmp[ia + iw][0];
                Uin_is[ic][1] = Utmp[ia + iw][1];
            }

            /*
             * Build tau_mat and Amp_mat for current ih.
             */
            for (int is = 0; is < ctx->ns; is++)
            {

                int c =
                    (int)((ctx->s0 + is * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);

                if (c < 0 || c >= ctx->nsg)
                {
                    ERROR(("Source table too small."));
                }

                int left =
                    (is - ctx->aper < 0)
                        ? 0
                        : is - ctx->aper;

                int right =
                    (is + ctx->aper > ctx->ns - 1)
                        ? ctx->ns - 1
                        : is + ctx->aper;

                for (int ic = left; ic <= right; ic++)
                {

                    int cc =
                        (int)((ctx->s0 + ic * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);

                    if (cc < 0 || cc >= ctx->nsg)
                    {
                        ERROR(("Source table too small."));
                    }

                    float coef =
                        taper_weight(left, ic, right, ctx->tap, 1);

                    float tau = ctx->stable[cc][c];

                    float dist =
                        ctx->sdatum * ctx->sdatum + (ic - is) * (ic - is) * ctx->ds * ctx->ds;

                    float w =
                        coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;

                    /*
                     * Frequency-domain operator anti-aliasing on source summation.
                     */
                    float slope =
                        table_slope_first_index(ctx->stable,
                                                ctx->nsg,
                                                cc,
                                                c,
                                                ctx->dsg);

                    float aa =
                        antialias_weight_iw(ctx,
                                            iw,
                                            slope * fabsf(ctx->ds));

                    w *= aa;

                    if (fabsf(w) < 1e-20f)
                        continue;

                    /*
                     * Build H(tau, iw) using the same cache strategy
                     * as the original direct-summation code.
                     */
                    size_t key =
                        (size_t)cc * (size_t)ctx->nsg + (size_t)c;

                    if (!cached[key])
                    {
                        build_H_iw(tau,
                                   ctx,
                                   iw,
                                   F,
                                   &hr_cache[key],
                                   &hi_cache[key]);

                        cached[key] = 1;
                    }

                    float hr = hr_cache[key];
                    float hi = hi_cache[key];

                    /*
                     * Phase-amplitude separation:
                     *
                     *     Amp = w * H(tau,iw) * exp(+i * omega * tau).
                     */
                    float phase = omega * tau;
                    float cosp = cosf(phase);
                    float sinp = sinf(phase);

                    float amp_r = w * (hr * cosp - hi * sinp);
                    float amp_i = w * (hr * sinp + hi * cosp);

                    size_t imat =
                        (size_t)is * (size_t)ctx->ns + (size_t)ic;

                    tau_mat[is][ic] = tau;

                    Amp_mat[imat][0] = amp_r;
                    Amp_mat[imat][1] = amp_i;
                }
            }

            /*
             * Direct summation for verification:
             *
             *     Uout_is[is]
             *       = sum_ic Amp_mat[is,ic]
             *                * exp(-i * omega * tau_mat[is][ic])
             *                * Uin_is[ic].
             *
             * Later this block can be replaced by a butterfly algorithm.
             */
            // for (int is = 0; is < ctx->ns; is++) {
            //     for (int ic = 0; ic < ctx->ns; ic++) {

            //         size_t imat =
            //             (size_t)is * (size_t)ctx->ns + (size_t)ic;

            //         float amp_r = Amp_mat[imat][0];
            //         float amp_i = Amp_mat[imat][1];

            //         if (fabsf(amp_r) < 1e-20f &&
            //             fabsf(amp_i) < 1e-20f) {
            //             continue;
            //         }

            //         float tau = tau_mat[is][ic];

            //         float phase = omega * tau;
            //         float cosp = cosf(phase);
            //         float sinp = sinf(phase);

            //         float kr = amp_r * cosp + amp_i * sinp;
            //         float ki = amp_i * cosp - amp_r * sinp;

            //         float ar = Uin_is[ic][0];
            //         float ai = Uin_is[ic][1];

            //         Uout_is[is][0] += kr * ar - ki * ai;
            //         Uout_is[is][1] += kr * ai + ki * ar;
            //     }
            // }
            butterfly_apply_1d_phase_amp(ctx->ns,
                                         tau_mat,
                                         Amp_mat,
                                         Uin_is,
                                         Uout_is,
                                         omega, ctx->p, ctx->n_leaf);

            /*
             * Butterfly replacement interface:
             *
             *     butterfly_source_apply(ctx,
             *                            iw,
             *                            tau_mat,
             *                            Amp_mat,
             *                            Uin_is,
             *                            Uout_is);
             */

            /*
             * Write current-ih output vector back to Uout.
             */
            for (int is = 0; is < ctx->ns; is++)
            {

                size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);

                Uout[ib + iw][0] += Uout_is[is][0];
                Uout[ib + iw][1] += Uout_is[is][1];
            }
        }
    }

    fftwf_free(Amp_mat);
    fftwf_free(Uin_is);
    fftwf_free(Uout_is);

    free2float(tau_mat);

    free(hr_cache);
    free(hi_cache);
    free(cached);
    free(F);
}

static int frequency_progress_step(int nw)
{
    int step = nw / 20; /* approximately 5% */
    if (step < 50)
        step = 50;
    if (nw <= 50)
        step = 1;
    return step;
}

static void print_frequency_progress(const char *stage, int iw, int nw, int verb)
{
    if (!verb)
        return;

    int step = frequency_progress_step(nw);
    if (!(iw == 0 || iw == nw - 1 || ((iw + 1) % step == 0)))
        return;

#ifdef _OPENMP
#pragma omp critical(freq_progress_print)
#endif
    INFO(("Processing %s frequency %d of %d.", stage, iw + 1, nw));
}

int main(int argc, char **argv)
{
    auto t_start = std::chrono::steady_clock::now();

    se_par_init(argc, argv);

    int should_datum = 0;
    int auto_datum = 1;
    int use_bf = 0;
    char *in_f = NULL, *out_f = NULL, *sgreen_f = NULL, *rgreen_f = NULL, *model_f = NULL, *interm_f = NULL;
    sep_t *in, *out, *sgreen, *rgreen, *model, *interm = NULL;
    float **shot = NULL, **stable, **rtable;

    if (!se_have_par("input_file"))
        ERROR(("Need input_file="));
    else
        in_f = se_get_par_str("input_file");
    if (!se_have_par("output_file"))
        ERROR(("Need output_file="));
    else
        out_f = se_get_par_str("output_file");
    if (!se_have_par("sgreen_file"))
        ERROR(("Need sgreen_file="));
    else
        sgreen_f = se_get_par_str("sgreen_file");
    if (!se_have_par("rgreen_file"))
        ERROR(("Need rgreen_file="));
    else
        rgreen_f = se_get_par_str("rgreen_file");
    if (!se_have_par("model_file"))
        ERROR(("Need model_file="));
    else
        model_f = se_get_par_str("model_file");

    in = sep_open(in_f, SEP_READ, 0);
    out = sep_open(out_f, SEP_WRITE, 0);
    model = sep_open(model_f, SEP_READ, 0);

    FreqContext ctx;
    ctx.verb = se_have_par("verb") ? se_get_par_int("verb") : 1;
    ctx.cmp = se_have_par("cmp") ? se_get_par_int("cmp") : 1;
    if (ctx.cmp != 0 && ctx.cmp != 1)
        ERROR(("cmp must be 0 or 1."));

    auto_datum = se_have_par("auto_datum") ? se_get_par_int("auto_datum") : 1;
    if (auto_datum != 0 && auto_datum != 1)
        ERROR(("auto_datum must be 0 or 1."));
    if (auto_datum)
    {
        if (sep_have_hdr_int(model, "should_datum"))
            should_datum = sep_get_hdr_int(model, "should_datum", 0);
        else
            ERROR(("Need should_datum in model_file header when auto_datum=1"));

        ctx.sdatum = (float)(model->headers->o[0] + should_datum * model->headers->d[0]);
        ctx.rdatum = ctx.sdatum;
    }
    else
    {
        if (!se_have_par("sdatum"))
            ERROR(("Need sdatum= when auto_datum=0"));
        if (!se_have_par("rdatum"))
            ERROR(("Need rdatum= when auto_datum=0"));
        ctx.sdatum = se_get_par_float("sdatum");
        ctx.rdatum = se_get_par_float("rdatum");
    }

    use_bf = se_have_par("use_bf") ? se_get_par_int("use_bf") : 1;
    ctx.p = se_have_par("bf_p") ? se_get_par_int("bf_p") : 12;
    ctx.n_leaf = se_have_par("bf_n_leaf") ? se_get_par_int("bf_n_leaf") : 16;

    INFO(("Use butterfly: %s.", use_bf ? "yes" : "no"));
    INFO(("Butterfly parameters: p = %d, n_leaf = %d.", ctx.p, ctx.n_leaf));

    ctx.aper = se_have_par("aperture") ? se_get_par_int("aperture") : 50;
    ctx.tap = se_have_par("taper") ? se_get_par_int("taper") : 10;
    ctx.antialias = se_have_par("antialias") ? se_get_par_float("antialias") : 1.0f;
    if (ctx.antialias < 0.0f)
        ERROR(("antialias must be >= 0."));
    float length = se_have_par("length") ? se_get_par_float("length") : 0.025f;

    if (in->headers->ndim < 3)
        ERROR(("Input must be 3D."));
    ctx.nt = in->headers->n[0];
    ctx.nh = in->headers->n[1];
    ctx.ns = in->headers->n[2];
    ctx.dt = in->headers->d[0];
    ctx.h0 = in->headers->o[1];
    ctx.dh = in->headers->d[1];
    ctx.s0 = in->headers->o[2];
    ctx.ds = in->headers->d[2];

    sgreen = sep_open(sgreen_f, SEP_READ, 0);
    if (sgreen->headers->ndim < 2)
        ERROR(("Source Green's function must be 2D."));
    ctx.nsg = sgreen->headers->n[0];
    ctx.sg0 = sgreen->headers->o[0];
    ctx.dsg = sgreen->headers->d[0];
    stable = alloc2float(ctx.nsg, ctx.nsg);
    se_fsio_read_float(sgreen->data->io, stable[0], ctx.nsg * ctx.nsg);
    sep_close(sgreen);

    rgreen = sep_open(rgreen_f, SEP_READ, 0);
    if (rgreen->headers->ndim < 2)
        ERROR(("Receiver Green's function must be 2D."));
    ctx.nrg = rgreen->headers->n[0];
    ctx.rg0 = rgreen->headers->o[0];
    ctx.drg = rgreen->headers->d[0];
    rtable = alloc2float(ctx.nrg, ctx.nrg);
    se_fsio_read_float(rgreen->data->io, rtable[0], ctx.nrg * ctx.nrg);
    sep_close(rgreen);

    ctx.stable = stable;
    ctx.rtable = rtable;

    if (se_have_par("interm"))
    {
        interm_f = se_get_par_str("interm");
        interm = sep_open(interm_f, SEP_WRITE, 0);
    }

    ctx.nsam = (int)(length / ctx.dt) + 2;
    ctx.nfft = choose_nfft(ctx.nt, ctx.nsam, ctx.dt, ctx.nrg, ctx.nsg, rtable, stable);
    ctx.nw = ctx.nfft / 2 + 1;
    ctx.df = 1.0f / ((float)ctx.nfft * ctx.dt);

    if (ctx.verb)
    {
        INFO(("antialias = %g; set antialias=0 to disable frequency-domain operator anti-aliasing.", ctx.antialias));
    }

    float *rpad = (float *)fftwf_malloc((size_t)ctx.nfft * sizeof(float));
    fftwf_complex *spec = (fftwf_complex *)fftwf_malloc((size_t)ctx.nw * sizeof(fftwf_complex));
    fftwf_plan p_f = fftwf_plan_dft_r2c_1d(ctx.nfft, rpad, spec, FFTW_ESTIMATE);
    fftwf_plan p_b = fftwf_plan_dft_c2r_1d(ctx.nfft, spec, rpad, FFTW_ESTIMATE);

    fftwf_complex *Uin = (fftwf_complex *)fftwf_malloc((size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));
    fftwf_complex *Utmp = (fftwf_complex *)fftwf_malloc((size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));

    memset(Utmp, 0, (size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));

    shot = alloc2float(ctx.nt, ctx.nh);
    read_traces_to_freq(in->data->io, &ctx, shot, p_f, rpad, spec, Uin);
    free2float(shot);
    shot = NULL;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iw = 0; iw < ctx.nw; iw++)
    {
        print_frequency_progress("receiver-side", iw, ctx.nw, ctx.verb);
        if (!use_bf)
            apply_receiver_operator_freq_iw(&ctx, Uin, Utmp, iw);
        else
            apply_receiver_operator_freq_iw_bf(&ctx, Uin, Utmp, iw);
    }

    if (interm)
    {
        write_traces_from_freq(Utmp, &ctx, p_b, spec, rpad, interm->data->io, 0);
    }

    fftwf_free(Uin);
    Uin = NULL;

    fftwf_complex *Uout = (fftwf_complex *)fftwf_malloc((size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));
    memset(Uout, 0, (size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iw = 0; iw < ctx.nw; iw++)
    {
        print_frequency_progress("source-side", iw, ctx.nw, ctx.verb);
        if (!use_bf)
            apply_source_operator_freq_iw(&ctx, Utmp, Uout, iw);
        else
            apply_source_operator_freq_iw_bf(&ctx, Utmp, Uout, iw);
    }

    fftwf_free(Utmp);
    Utmp = NULL;

    out->headers->ndim = 3;
    out->headers->n[0] = ctx.nt;
    out->headers->n[1] = ctx.nh;
    out->headers->n[2] = ctx.ns;
    out->headers->d[0] = ctx.dt;
    out->headers->d[1] = ctx.dh;
    out->headers->d[2] = ctx.ds;
    out->headers->o[0] = 0.0f;
    out->headers->o[1] = ctx.h0;
    out->headers->o[2] = ctx.s0;

    write_traces_from_freq(Uout, &ctx, p_b, spec, rpad, out->data->io, 1);

    fftwf_destroy_plan(p_f);
    fftwf_destroy_plan(p_b);
    fftwf_free(rpad);
    fftwf_free(spec);
    fftwf_free(Uout);

    sep_close(in);
    sep_close(out);
    sep_close(model);
    if (interm)
        sep_close(interm);

    if (shot)
        free2float(shot);
    free2float(stable);
    free2float(rtable);

    auto t_end = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration<double>(t_end - t_start).count();
    INFO(("Done. Elapsed time: %.3f s.", elapsed_seconds));
    return 0;
}
