#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define IDX3(is, ih, iw, nh, nw) ((((size_t)(is) * (nh) + (ih)) * (nw)) + (iw))
typedef struct {
    float *hr, *hi, *er, *ei;
    unsigned char *cached;
    int n;
} HECache;

typedef struct {
    float kr, ki, er, ei, x0;
    int valid;
    size_t ia, ib;
} KernelValue;

typedef struct {
    int nt, nh, ns;
    int nsg, nrg;
    int aper, tap, cmp, verb;
    int nfft, nw, nsam;
    float dt, h0, dh, s0, ds, sg0, dsg, rg0, drg;
    float sdatum, rdatum;
    float **stable, **rtable;
    int use_butterfly;
    int use_edge_correction;
} FreqContext;

static void reverse_trace(int nt, int nh, int ns, float ***tr)
{
    INFO(("Reversing traces..."));
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < ns; is++) {
        for (int ih = 0; ih < nh; ih++) {
            for (int it = 0; it < nt / 2; it++) {
                float t = tr[is][ih][it];
                tr[is][ih][it] = tr[is][ih][nt - 1 - it];
                tr[is][ih][nt - 1 - it] = t;
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

static int choose_nfft(int nt, int nsam, float dt, int nrg, int nsg, float **rtable, float **stable)
{
    float max_tau = 0.0f;
    for (int i = 0; i < nrg * nrg; i++) if (rtable[0][i] > max_tau) max_tau = rtable[0][i];
    for (int i = 0; i < nsg * nsg; i++) if (stable[0][i] > max_tau) max_tau = stable[0][i];
    int pmax = (int)ceilf(max_tau / dt) + (nsam - 2);
    return next_pow2(nt + 2 * pmax + 2);
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

static void build_HE_iw(float tau, const FreqContext *ctx, int iw, float *F, float *hr, float *hi, float *er, float *ei)
{
    if (tau <= 0.0f || !isfinite(tau)) {
        *hr = 0.0f;
        *hi = 0.0f;
        *er = 0.0f;
        *ei = 0.0f;
        return;
    }

    build_filter_F(tau, ctx->dt, ctx->nsam, F);

    float q = tau / ctx->dt;
    int m = (int)ceilf(q);
    float delta = (float)m - q;

    float w = 2.0f * (float)M_PI * iw / ctx->nfft;
    float ar = 0.0f, ai = 0.0f;
    float br = 0.0f, bi = 0.0f;
    for (int k = 0; k <= ctx->nsam - 2; k++) {
        float fk = F[k] / ctx->dt;
        float p1 = (float)(m + k);
        float p2 = (float)(m + k - 1);
        ar += fk * ((1.0f - delta) * cosf(-w * p1) + delta * cosf(-w * p2));
        ai += fk * ((1.0f - delta) * sinf(-w * p1) + delta * sinf(-w * p2));
        br += -delta * fk * cosf(-w * p2);
        bi += -delta * fk * sinf(-w * p2);
    }
    *hr = ar;
    *hi = ai;
    *er = br;
    *ei = bi;
}

static void accumulate_pair_at_freq(float ar, float ai, float br, float bi, float w, float *or_, float *oi)
{
    *or_ += w * (ar * br - ai * bi);
    *oi += w * (ar * bi + ai * br);
}

static void get_he_from_cache_or_build(HECache *cache, size_t key, float tau, const FreqContext *ctx, int iw, float *F,
                                       float *hr, float *hi, float *er, float *ei)
{
    if (!cache->cached[key]) {
        build_HE_iw(tau, ctx, iw, F, &cache->hr[key], &cache->hi[key], &cache->er[key], &cache->ei[key]);
        cache->cached[key] = 1;
    }
    *hr = cache->hr[key];
    *hi = cache->hi[key];
    *er = cache->er[key];
    *ei = cache->ei[key];
}

static void apply_kernel_value(const KernelValue *kv, const fftwf_complex *Uin, fftwf_complex *Uout, int iw, int use_edge_correction)
{
    if (!kv->valid) return;
    float ar = Uin[kv->ia + iw][0], ai = Uin[kv->ia + iw][1];
    accumulate_pair_at_freq(ar, ai, kv->kr, kv->ki, 1.0f, &Uout[kv->ib + iw][0], &Uout[kv->ib + iw][1]);
    if (use_edge_correction) {
        Uout[kv->ib + iw][0] += kv->er * kv->x0;
        Uout[kv->ib + iw][1] += kv->ei * kv->x0;
    }
}
static inline float taper_weight(int left, int center, int right, int tap, int stride);

static KernelValue receiver_kernel(const FreqContext *ctx, int iw, int is, int ih, int ic,
                                   const float *Xin0, HECache *cache, float *F)
{
    KernelValue kv = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0};
    int c = (int)((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f) + ctx->h0 + ih * ctx->dh - ctx->rg0) / ctx->drg) + 0.5f);
    int cc = (int)((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f) + ctx->h0 + ic * ctx->dh - ctx->rg0) / ctx->drg) + 0.5f);
    if (c < 0 || c >= ctx->nrg || cc < 0 || cc >= ctx->nrg) ERROR(("Receiver table too small."));
    int left = (ih - ctx->aper < 0) ? 0 : ih - ctx->aper;
    int right = (ih + ctx->aper > ctx->nh - 1) ? ctx->nh - 1 : ih + ctx->aper;
    float coef = taper_weight(left, ic, right, ctx->tap, 1);
    float tau = ctx->rtable[cc][c];
    float dist = ctx->rdatum * ctx->rdatum + (ic - ih) * (ic - ih) * ctx->dh * ctx->dh;
    float w = coef / M_PI * ctx->dh * ctx->rdatum * tau / dist;
    if (fabsf(w) < 1e-20f) return kv;
    float hr, hi, er, ei;
    size_t key = (size_t)cc * ctx->nrg + c;
    get_he_from_cache_or_build(cache, key, tau, ctx, iw, F, &hr, &hi, &er, &ei);
    kv.kr = w * hr; kv.ki = w * hi; kv.er = w * er; kv.ei = w * ei;
    kv.x0 = Xin0[is * ctx->nh + ic];
    kv.ia = IDX3(is, ic, 0, ctx->nh, ctx->nw);
    kv.ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);
    kv.valid = 1;
    return kv;
}

static inline float taper_weight(int left, int center, int right, int tap, int stride)
{
    if (tap <= 0) return 1.0f;
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

static KernelValue source_kernel(const FreqContext *ctx, int iw, int is, int ih, int ic, int hh,
                                 int left, int right, int jump, const float *Xtmp0, HECache *cache, float *F)
{
    KernelValue kv = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0};
    int c = (int)((ctx->s0 + is * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
    if (c < 0 || c >= ctx->nsg) ERROR(("Source table too small."));
    int cc = (int)((ctx->s0 + ic * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
    if (cc < 0 || cc >= ctx->nsg) ERROR(("Source table too small."));
    float coef = taper_weight(left, ic, right, ctx->tap, jump);
    float tau = ctx->stable[cc][c];
    float dist = ctx->sdatum * ctx->sdatum + (ic - is) * (ic - is) * ctx->ds * ctx->ds;
    float w = coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;
    if (fabsf(w) < 1e-20f) return kv;
    float hr, hi, er, ei;
    size_t key = (size_t)cc * ctx->nsg + c;
    get_he_from_cache_or_build(cache, key, tau, ctx, iw, F, &hr, &hi, &er, &ei);
    kv.kr = w * hr; kv.ki = w * hi; kv.er = w * er; kv.ei = w * ei;
    kv.x0 = Xtmp0[ic * ctx->nh + hh];
    kv.ia = IDX3(ic, hh, 0, ctx->nh, ctx->nw);
    kv.ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);
    kv.valid = 1;
    return kv;
}

typedef KernelValue (*ReceiverKernelCb)(const FreqContext *, int, int, int, int, const float *, HECache *, float *);
typedef KernelValue (*SourceKernelCb)(const FreqContext *, int, int, int, int, int, int, int, int, const float *, HECache *, float *);

static void apply_receiver_direct_operator_common_iw(const FreqContext *ctx, fftwf_complex *Uin, fftwf_complex *Utmp,
                                                     int iw, const float *Xin0, ReceiverKernelCb kcb)
{
    float *F = (float *)malloc((size_t)(ctx->nsam - 1) * sizeof(float));
    HECache cache;
    cache.n = ctx->nrg * ctx->nrg;
    cache.hr = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.hi = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.er = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.ei = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.cached = (unsigned char *)calloc((size_t)cache.n, sizeof(unsigned char));
    for (int is = 0; is < ctx->ns; is++) {
        if (iw == 0 && ctx->verb) INFO(("Processing common-shot gather %d of %d.", is + 1, ctx->ns));
        for (int ih = 0; ih < ctx->nh; ih++) {
            int left = (ih - ctx->aper < 0) ? 0 : ih - ctx->aper;
            int right = (ih + ctx->aper > ctx->nh - 1) ? ctx->nh - 1 : ih + ctx->aper;
            for (int ic = left; ic <= right; ic++) {
                KernelValue kv = kcb(ctx, iw, is, ih, ic, Xin0, &cache, F);
                apply_kernel_value(&kv, Uin, Utmp, iw, ctx->use_edge_correction);
            }
        }
    }
    free(cache.hr); free(cache.hi); free(cache.er); free(cache.ei); free(cache.cached);
    free(F);
}

static void apply_source_direct_operator_common_iw(const FreqContext *ctx, fftwf_complex *Utmp, fftwf_complex *Uout,
                                                   int iw, const float *Xtmp0, SourceKernelCb kcb)
{
    float *F = (float *)malloc((size_t)(ctx->nsam - 1) * sizeof(float));
    HECache cache;
    cache.n = ctx->nsg * ctx->nsg;
    cache.hr = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.hi = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.er = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.ei = (float *)malloc((size_t)cache.n * sizeof(float));
    cache.cached = (unsigned char *)calloc((size_t)cache.n, sizeof(unsigned char));
    if (ctx->cmp == 1) {
        float s = fabsf((ctx->ns - 1) * ctx->ds), h = fabsf((ctx->nh - 1) * ctx->dh), dr;
        int jump = (fabsf(ctx->ds) >= fabsf(ctx->dh)) ? 1 : (int)(ctx->dh / ctx->ds + 0.5f);
        dr = (fabsf(ctx->ds) >= fabsf(ctx->dh)) ? fabsf(ctx->dh) : fabsf(ctx->ds);
        int nr = (int)((s + h) / dr + 1.5f);
        for (int ir = 0; ir < nr; ir++) {
            if (iw == 0 && ctx->verb) INFO(("Processing common-receiver gather %d of %d.", ir + 1, nr));
            float r = ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->dh <= 0.f) ? -1.f : 0.f) * h;
            int sleft = (int)((ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->ds <= 0.f) ? 0.f : -1.f) * h) / ctx->ds + 0.5f);
            int sright = (int)((ir * dr + ((ctx->ds <= 0.f) ? -1.f : 0.f) * s + ((ctx->ds <= 0.f) ? -1.f : 0.f) * h) / ctx->ds + 0.5f);
            if (sleft < 0) sleft = 0;
            if (sright > ctx->ns - 1) sright = ctx->ns - 1;
            int left = (int)((r - sleft * ctx->ds) / ctx->dh + 0.5f);
            if (left < 0 || left > ctx->nh - 1) sleft++;
            int right = (int)((r - sright * ctx->ds) / ctx->dh + 0.5f);
            if (right < 0 || right > ctx->nh - 1) sright--;
            for (int is = sleft; is <= sright; is += jump) {
                int ih = (int)((r - is * ctx->ds) / ctx->dh + 0.5f);
                left = (is - jump * ctx->aper < sleft) ? sleft : is - jump * ctx->aper;
                right = (is + jump * ctx->aper > sright) ? sright : is + jump * ctx->aper;
                for (int ic = left; ic <= right; ic += jump) {
                    int hh = (int)((r - ic * ctx->ds) / ctx->dh + 0.5f);
                    KernelValue kv = kcb(ctx, iw, is, ih, ic, hh, left, right, jump, Xtmp0, &cache, F);
                    apply_kernel_value(&kv, Utmp, Uout, iw, ctx->use_edge_correction);
                }
            }
        }
    } else {
        for (int ih = 0; ih < ctx->nh; ih++) {
            if (iw == 0 && ctx->verb) INFO(("Processing absolute-coordinate receiver gather %d of %d.", ih + 1, ctx->nh));
            for (int is = 0; is < ctx->ns; is++) {
                int left = (is - ctx->aper < 0) ? 0 : is - ctx->aper;
                int right = (is + ctx->aper > ctx->ns - 1) ? ctx->ns - 1 : is + ctx->aper;
                for (int ic = left; ic <= right; ic++) {
                    KernelValue kv = kcb(ctx, iw, is, ih, ic, ih, left, right, 1, Xtmp0, &cache, F);
                    apply_kernel_value(&kv, Utmp, Uout, iw, ctx->use_edge_correction);
                }
            }
        }
    }
    free(cache.hr); free(cache.hi); free(cache.er); free(cache.ei); free(cache.cached);
    free(F);
}

static void fft_traces_to_freq(float ***tr, const FreqContext *ctx, fftwf_plan p_f, float *pad, fftwf_complex *spec, fftwf_complex *U)
{
    for (int is = 0; is < ctx->ns; is++) {
        for (int ih = 0; ih < ctx->nh; ih++) {
            memset(pad, 0, (size_t)ctx->nfft * sizeof(float));
            memcpy(pad, tr[is][ih], (size_t)ctx->nt * sizeof(float));
            fftwf_execute(p_f);
            memcpy(&U[IDX3(is, ih, 0, ctx->nh, ctx->nw)], spec, (size_t)ctx->nw * sizeof(fftwf_complex));
        }
    }
}

static void ifft_traces_from_freq(fftwf_complex *U, const FreqContext *ctx, fftwf_plan p_b, fftwf_complex *spec, float *pad, float ***tr)
{
    for (int is = 0; is < ctx->ns; is++) {
        for (int ih = 0; ih < ctx->nh; ih++) {
            memcpy(spec, &U[IDX3(is, ih, 0, ctx->nh, ctx->nw)], (size_t)ctx->nw * sizeof(fftwf_complex));
            fftwf_execute(p_b);
            for (int it = 0; it < ctx->nt; it++) tr[is][ih][it] = pad[it] / ctx->nfft;
        }
    }
}

static void apply_receiver_direct_operator_freq_iw(const FreqContext *ctx, fftwf_complex *Uin, fftwf_complex *Utmp, int iw, const float *Xin0)
{
    apply_receiver_direct_operator_common_iw(ctx, Uin, Utmp, iw, Xin0, receiver_kernel);
}

static void apply_source_direct_operator_freq_iw(const FreqContext *ctx, fftwf_complex *Utmp, fftwf_complex *Uout, int iw, const float *Xtmp0)
{
    apply_source_direct_operator_common_iw(ctx, Utmp, Uout, iw, Xtmp0, source_kernel);
}

int main(int argc, char **argv)
{
    se_par_init(argc, argv);

    int should_datum = 0;
    char *in_f = NULL, *out_f = NULL, *sgreen_f = NULL, *rgreen_f = NULL, *model_f = NULL, *interm_f = NULL;
    sep_t *in, *out, *sgreen, *rgreen, *model, *interm = NULL;
    float ***tr_in, ***tr_out, ***tr_mid, **stable, **rtable;

    if (!se_have_par("input_file")) ERROR(("Need input_file=")); else in_f = se_get_par_str("input_file");
    if (!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");
    if (!se_have_par("sgreen_file")) ERROR(("Need sgreen_file=")); else sgreen_f = se_get_par_str("sgreen_file");
    if (!se_have_par("rgreen_file")) ERROR(("Need rgreen_file=")); else rgreen_f = se_get_par_str("rgreen_file");
    if (!se_have_par("model_file")) ERROR(("Need model_file=")); else model_f = se_get_par_str("model_file");

    in = sep_open(in_f, SEP_READ, 0);
    out = sep_open(out_f, SEP_WRITE, 0);
    model = sep_open(model_f, SEP_READ, 0);

    FreqContext ctx;
    ctx.verb = se_have_par("verb") ? se_get_par_int("verb") : 1;
    ctx.cmp = se_have_par("cmp") ? se_get_par_int("cmp") : 1;
    ctx.use_butterfly = se_have_par("use_butterfly") ? se_get_par_int("use_butterfly") : 0;
    ctx.use_edge_correction = se_have_par("use_edge_correction") ? se_get_par_int("use_edge_correction") : 1;
    if (ctx.cmp != 0 && ctx.cmp != 1) ERROR(("cmp must be 0 or 1."));
    if (ctx.use_butterfly && ctx.verb) {
        INFO(("use_butterfly=1 requested; butterfly operator placeholder path is reserved and currently falls back to direct frequency-domain summation."));
    }

    if (sep_have_hdr_int(model, "should_datum")) should_datum = sep_get_hdr_int(model, "should_datum", 0);
    else ERROR(("Need should_datum in model_file header"));

    ctx.aper = se_have_par("aperture") ? se_get_par_int("aperture") : 50;
    ctx.tap = se_have_par("taper") ? se_get_par_int("taper") : 10;
    float length = se_have_par("length") ? se_get_par_float("length") : 0.025f;

    if (in->headers->ndim < 3) ERROR(("Input must be 3D."));
    ctx.nt = in->headers->n[0];
    ctx.nh = in->headers->n[1];
    ctx.ns = in->headers->n[2];
    ctx.dt = in->headers->d[0];
    ctx.h0 = in->headers->o[1];
    ctx.dh = in->headers->d[1];
    ctx.s0 = in->headers->o[2];
    ctx.ds = in->headers->d[2];

    ctx.sdatum = (float)(model->headers->o[0] + should_datum * model->headers->d[0]);
    ctx.rdatum = ctx.sdatum;

    tr_in = alloc3float(ctx.nt, ctx.nh, ctx.ns);
    tr_out = alloc3float(ctx.nt, ctx.nh, ctx.ns);
    tr_mid = alloc3float(ctx.nt, ctx.nh, ctx.ns);
    se_fsio_read_float(in->data->io, tr_in[0][0], ctx.nt * ctx.nh * ctx.ns);
    reverse_trace(ctx.nt, ctx.nh, ctx.ns, tr_in);

    sgreen = sep_open(sgreen_f, SEP_READ, 0);
    if (sgreen->headers->ndim < 2) ERROR(("Source Green's function must be 2D."));
    ctx.nsg = sgreen->headers->n[0];
    ctx.sg0 = sgreen->headers->o[0];
    ctx.dsg = sgreen->headers->d[0];
    stable = alloc2float(ctx.nsg, ctx.nsg);
    se_fsio_read_float(sgreen->data->io, stable[0], ctx.nsg * ctx.nsg);
    sep_close(sgreen);

    rgreen = sep_open(rgreen_f, SEP_READ, 0);
    if (rgreen->headers->ndim < 2) ERROR(("Receiver Green's function must be 2D."));
    ctx.nrg = rgreen->headers->n[0];
    ctx.rg0 = rgreen->headers->o[0];
    ctx.drg = rgreen->headers->d[0];
    rtable = alloc2float(ctx.nrg, ctx.nrg);
    se_fsio_read_float(rgreen->data->io, rtable[0], ctx.nrg * ctx.nrg);
    sep_close(rgreen);

    ctx.stable = stable;
    ctx.rtable = rtable;

    if (se_have_par("interm")) {
        interm_f = se_get_par_str("interm");
        interm = sep_open(interm_f, SEP_WRITE, 0);
    }

    ctx.nsam = (int)(length / ctx.dt) + 2;
    ctx.nfft = choose_nfft(ctx.nt, ctx.nsam, ctx.dt, ctx.nrg, ctx.nsg, rtable, stable);
    ctx.nw = ctx.nfft / 2 + 1;

    float *rpad = (float *)fftwf_malloc((size_t)ctx.nfft * sizeof(float));
    fftwf_complex *spec = (fftwf_complex *)fftwf_malloc((size_t)ctx.nw * sizeof(fftwf_complex));
    fftwf_plan p_f = fftwf_plan_dft_r2c_1d(ctx.nfft, rpad, spec, FFTW_ESTIMATE);
    fftwf_plan p_b = fftwf_plan_dft_c2r_1d(ctx.nfft, spec, rpad, FFTW_ESTIMATE);

    fftwf_complex *Uin = (fftwf_complex *)fftwf_malloc((size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));
    fftwf_complex *Utmp = (fftwf_complex *)fftwf_malloc((size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));
    fftwf_complex *Uout = (fftwf_complex *)fftwf_malloc((size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));
    float *Xin0 = (float *)malloc((size_t)ctx.ns * ctx.nh * sizeof(float));
    float *Xtmp0 = (float *)malloc((size_t)ctx.ns * ctx.nh * sizeof(float));

    memset(Utmp, 0, (size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));
    memset(Uout, 0, (size_t)ctx.ns * ctx.nh * ctx.nw * sizeof(fftwf_complex));

    for (int is = 0; is < ctx.ns; is++) for (int ih = 0; ih < ctx.nh; ih++) Xin0[is * ctx.nh + ih] = tr_in[is][ih][0];
    fft_traces_to_freq(tr_in, &ctx, p_f, rpad, spec, Uin);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iw = 0; iw < ctx.nw; iw++) {
        apply_receiver_direct_operator_freq_iw(&ctx, Uin, Utmp, iw, Xin0);
    }

    ifft_traces_from_freq(Utmp, &ctx, p_b, spec, rpad, tr_mid);
    for (int is = 0; is < ctx.ns; is++) for (int ih = 0; ih < ctx.nh; ih++) Xtmp0[is * ctx.nh + ih] = tr_mid[is][ih][0];
    fft_traces_to_freq(tr_mid, &ctx, p_f, rpad, spec, Utmp);

    if (interm) {
        memcpy(tr_out[0][0], tr_mid[0][0], (size_t)ctx.nt * ctx.nh * ctx.ns * sizeof(float));
        se_fsio_write_float(interm->data->io, tr_out[0][0], ctx.nt * ctx.nh * ctx.ns);
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iw = 0; iw < ctx.nw; iw++) {
        apply_source_direct_operator_freq_iw(&ctx, Utmp, Uout, iw, Xtmp0);
    }

    ifft_traces_from_freq(Uout, &ctx, p_b, spec, rpad, tr_in);

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

    reverse_trace(ctx.nt, ctx.nh, ctx.ns, tr_in);
    se_fsio_write_float(out->data->io, tr_in[0][0], ctx.nt * ctx.nh * ctx.ns);

    fftwf_destroy_plan(p_f);
    fftwf_destroy_plan(p_b);
    fftwf_free(rpad);
    fftwf_free(spec);
    fftwf_free(Uin);
    fftwf_free(Utmp);
    fftwf_free(Uout);
    free(Xin0);
    free(Xtmp0);

    sep_close(in);
    sep_close(out);
    sep_close(model);
    if (interm) sep_close(interm);

    free3float(tr_in);
    free3float(tr_out);
    free3float(tr_mid);
    free2float(stable);
    free2float(rtable);

    INFO(("Done."));
    return 0;
}
