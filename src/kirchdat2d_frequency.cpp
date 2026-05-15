#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if defined(SE_USE_OMP) || defined(_OPENMP)
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------------- */
/* This version evaluates the frequency-domain integral/summation directly.   */
/*                                                                            */
/* For each input trace d(t), delay tau, and fractional shift delta, it builds */
/* the same tau-dependent finite filter as the original time-domain code,     */
/* computes spectra by FFTW, and then explicitly evaluates                    */
/*                                                                            */
/*   u[n] = 1/N sum_w D_delta(w) H_tau(w) exp(i*w*n*dt).                      */
/*                                                                            */
/* No inverse FFT is used for the output trace. The last sum over frequencies  */
/* is written explicitly, so the structure can later be replaced by a          */
/* butterfly summation.                                                        */
/* ------------------------------------------------------------------------- */

typedef struct kirdat_integral_workspace_t {
    float *h;              /* zero-padded tau-dependent filter */
    float *x;              /* zero-padded linearly interpolated trace */
    float *out;            /* direct frequency-summation result, length nt */

    fftwf_complex *H;      /* spectrum of h */
    fftwf_complex *X;      /* spectrum of x */

    fftwf_plan plan_h;
    fftwf_plan plan_x;
} kirdat_integral_workspace_t;

static int   g_nt = 0;
static float g_dt = 0.0f;
static int   g_nsam = 0;
static int   g_nfft = 0;
static int   g_nfreq = 0;
static int   g_iw_min = 0;
static int   g_iw_max = 0;
static int   g_mts = 1;
static kirdat_integral_workspace_t *g_ws = NULL;

static int nextpow2_int(int n)
{
    int p = 1;
    while (p < n) p <<= 1;
    return p;
}

static int current_thread_id(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

static float taper_weight(int i, int left, int right, int tap, int jump)
{
    float coef = 1.0f;
    float denom;

    if (tap <= 0) return 1.0f;
    if (jump <= 0) jump = 1;

    denom = (float)tap * (float)jump;

    coef *= (i - left  >= tap * jump) ? 1.0f : ((float)(i - left) / denom);
    coef *= (right - i >= tap * jump) ? 1.0f : ((float)(right - i) / denom);

    if (coef < 0.0f) coef = 0.0f;
    if (coef > 1.0f) coef = 1.0f;
    return coef;
}

static int first_output_index(float tau, float dt, int nt, float *delta)
{
    int it0;

    if (tau <= 0.0f) {
        *delta = 0.0f;
        return 0;
    }

    /* Same logic as the original loop: use the first sample with it*dt >= tau. */
    it0 = (int)ceilf(tau / dt - 1.0e-6f);

    if (it0 < 0) it0 = 0;
    if (it0 >= nt) {
        *delta = 0.0f;
        return nt;
    }

    *delta = (((float)it0) * dt - tau) / dt;
    if (*delta < 0.0f) *delta = 0.0f;
    if (*delta > 1.0f) *delta = 1.0f;

    return it0;
}

static void build_tau_filter(float tau, float *h)
{
    int k;

    memset(h, 0, sizeof(float) * g_nfft);

    if (tau <= 0.0f) return;

    /* Same filter construction as the original filt_set(tau). */
    h[0] = 0.0f;
    for (k = 1; k < g_nsam; k++) {
        const float ratio = (tau + k * g_dt) / tau;
        const float val = ratio * ratio - 1.0f;
        h[k] = sqrtf((val > 0.0f) ? val : 0.0f);
    }

    /* First derivative. */
    for (k = 0; k < g_nsam - 1; k++) {
        h[k] = h[k + 1] - h[k];
    }

    /* Second derivative. */
    for (k = g_nsam - 2; k > 0; k--) {
        h[k] = h[k] - h[k - 1];
    }

    /* Same scaling as the original kirdat_pick(): filt[k] / dt. */
    for (k = 0; k < g_nsam - 1; k++) {
        h[k] /= g_dt;
    }
}

static void build_interpolated_trace(const float *trace, float delta, float *x)
{
    int it;

    memset(x, 0, sizeof(float) * g_nfft);

    /*
     * Same linear interpolation as the original kirdat_pick():
     *   (1-delta)*trace[i] + delta*trace[i+1]
     *
     * For the last sample, trace[nt] is unavailable. We use trace[nt-1]
     * to avoid the out-of-bounds access present in the original pointwise code.
     */
    for (it = 0; it < g_nt - 1; it++) {
        x[it] = (1.0f - delta) * trace[it] + delta * trace[it + 1];
    }
    if (g_nt > 0) x[g_nt - 1] = trace[g_nt - 1];
}

void kirdat_integral_init(int nt, float dt, float length, float fmin, float fmax)
{
    int its;
    float nyq;

    g_nt = nt;
    g_dt = dt;
    g_nsam = (int)(length / g_dt) + 2;

    /* Use zero padding so the frequency summation represents linear convolution. */
    g_nfft = nextpow2_int(g_nt + g_nsam + 2);
    g_nfreq = g_nfft / 2 + 1;

    nyq = 0.5f / g_dt;
    if (fmin < 0.0f) fmin = 0.0f;
    if (fmax <= 0.0f || fmax > nyq) fmax = nyq;
    if (fmin > fmax) ERROR(("Invalid frequency band: fmin=%g, fmax=%g", fmin, fmax));

    g_iw_min = (int)ceilf(fmin * (float)g_nfft * g_dt - 1.0e-6f);
    g_iw_max = (int)floorf(fmax * (float)g_nfft * g_dt + 1.0e-6f);

    if (g_iw_min < 0) g_iw_min = 0;
    if (g_iw_max > g_nfreq - 1) g_iw_max = g_nfreq - 1;
    if (g_iw_min > g_iw_max) ERROR(("The selected frequency band contains no DFT samples."));

#ifdef _OPENMP
    g_mts = omp_get_max_threads();
#else
    g_mts = 1;
#endif

    g_ws = (kirdat_integral_workspace_t *)calloc(g_mts, sizeof(kirdat_integral_workspace_t));
    if (g_ws == NULL) ERROR(("Failed to allocate frequency-integral workspaces."));

    for (its = 0; its < g_mts; its++) {
        kirdat_integral_workspace_t *w = &g_ws[its];

        w->h = (float *)fftwf_malloc(sizeof(float) * g_nfft);
        w->x = (float *)fftwf_malloc(sizeof(float) * g_nfft);
        w->out = (float *)malloc(sizeof(float) * g_nt);

        w->H = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * g_nfreq);
        w->X = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * g_nfreq);

        if (w->h == NULL || w->x == NULL || w->out == NULL || w->H == NULL || w->X == NULL) {
            ERROR(("Failed to allocate FFTW arrays."));
        }

        memset(w->h, 0, sizeof(float) * g_nfft);
        memset(w->x, 0, sizeof(float) * g_nfft);
        memset(w->out, 0, sizeof(float) * g_nt);

        /* Planning is done outside OpenMP parallel regions. */
        w->plan_h = fftwf_plan_dft_r2c_1d(g_nfft, w->h, w->H, FFTW_ESTIMATE);
        w->plan_x = fftwf_plan_dft_r2c_1d(g_nfft, w->x, w->X, FFTW_ESTIMATE);

        if (w->plan_h == NULL || w->plan_x == NULL) {
            ERROR(("Failed to create FFTW plans."));
        }
    }

    INFO(("Frequency-integral filter initialized: nt=%d, nsam=%d, nfft=%d, nfreq=%d, iw=[%d,%d], threads=%d.",
          g_nt, g_nsam, g_nfft, g_nfreq, g_iw_min, g_iw_max, g_mts));
}

void kirdat_integral_close(void)
{
    int its;

    if (g_ws == NULL) return;

    for (its = 0; its < g_mts; its++) {
        kirdat_integral_workspace_t *w = &g_ws[its];

        if (w->plan_h != NULL) fftwf_destroy_plan(w->plan_h);
        if (w->plan_x != NULL) fftwf_destroy_plan(w->plan_x);

        if (w->h != NULL) fftwf_free(w->h);
        if (w->x != NULL) fftwf_free(w->x);
        if (w->out != NULL) free(w->out);

        if (w->H != NULL) fftwf_free(w->H);
        if (w->X != NULL) fftwf_free(w->X);
    }

    free(g_ws);
    g_ws = NULL;

    fftwf_cleanup();
}

const float *kirdat_filter_trace_integral(const float *trace, float tau, float delta)
{
    int iw, n;
    const int its = current_thread_id();
    kirdat_integral_workspace_t *w = &g_ws[its];

    if (tau <= 0.0f) {
        memset(w->out, 0, sizeof(float) * g_nt);
        return w->out;
    }

    build_tau_filter(tau, w->h);
    build_interpolated_trace(trace, delta, w->x);

    fftwf_execute(w->plan_h);
    fftwf_execute(w->plan_x);

    memset(w->out, 0, sizeof(float) * g_nt);

    /*
     * Direct inverse real-DFT summation:
     *
     *   out[n] = 1/N * {P[0] + P[N/2](-1)^n +
     *                    2 sum_{iw=1}^{N/2-1} Re[P[iw] exp(i*2*pi*iw*n/N)] }
     *
     * where P[iw] = X[iw] * H[iw]. No inverse FFT is used here.
     */
    for (iw = g_iw_min; iw <= g_iw_max; iw++) {
        const float xr = w->X[iw][0];
        const float xi = w->X[iw][1];
        const float hr = w->H[iw][0];
        const float hi = w->H[iw][1];

        const float pr0 = xr * hr - xi * hi;
        const float pi0 = xr * hi + xi * hr;

        const int is_dc = (iw == 0);
        const int is_nyq = ((g_nfft % 2 == 0) && (iw == g_nfft / 2));
        const float mirror_weight = (is_dc || is_nyq) ? 1.0f : 2.0f;

        const float theta = 2.0f * (float)M_PI * (float)iw / (float)g_nfft;
        const float wr = cosf(theta);
        const float wi = sinf(theta);

        float er = 1.0f;
        float ei = 0.0f;

        for (n = 0; n < g_nt; n++) {
            w->out[n] += mirror_weight * (pr0 * er - pi0 * ei);

            /* e^{i theta (n+1)} = e^{i theta n} e^{i theta} */
            {
                const float ter = er * wr - ei * wi;
                const float tei = er * wi + ei * wr;
                er = ter;
                ei = tei;
            }
        }
    }

    for (n = 0; n < g_nt; n++) {
        w->out[n] /= (float)g_nfft;
    }

    return w->out;
}

static void reverse_trace(int nt, int nh, int ns, float ***tr)
{
    INFO(("Reversing traces..."));
#ifdef _OPENMP
#pragma omp parallel for
#endif
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

int main(int argc, char *argv[])
{
    se_par_init(argc, argv);

    int verb;
    int it, it0, nt, ih, nh, is, ns, nsg, nrg;
    int left, right, ic, aper, c, cc, hh;
    int ir, nr, jump, sleft, sright, tap;
    int cmp = 1; /* cmp=1: h-axis is relative receiver coordinate; cmp=0: h-axis is absolute receiver coordinate */

    float sdatum = 0.0f;
    float rdatum = 0.0f;
    float length, dt, h0, dh, s0, ds, sg0, dsg, rg0, drg;
    float dist, tau, delta, r, dr, s, h, coef, amp;
    float fmin, fmax;

    float ***tr_in = NULL;
    float ***tr_out = NULL;
    float **stable = NULL;
    float **rtable = NULL;

    sep_t *in = NULL;
    sep_t *out = NULL;
    sep_t *sgreen = NULL;
    sep_t *rgreen = NULL;
    sep_t *interm = NULL;
    sep_t *model = NULL;

    int should_datum = 0;
    char *in_f = NULL;
    char *out_f = NULL;
    char *sgreen_f = NULL;
    char *rgreen_f = NULL;
    char *interm_f = NULL;
    char *model_f = NULL;

    if (!se_have_par("input_file"))  ERROR(("Need input_file="));  else in_f = se_get_par_str("input_file");
    if (!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");
    if (!se_have_par("sgreen_file")) ERROR(("Need sgreen_file=")); else sgreen_f = se_get_par_str("sgreen_file");
    if (!se_have_par("rgreen_file")) ERROR(("Need rgreen_file=")); else rgreen_f = se_get_par_str("rgreen_file");
    if (!se_have_par("model_file"))  ERROR(("Need model_file="));  else model_f = se_get_par_str("model_file");

    in = sep_open(in_f, SEP_READ, 0);
    out = sep_open(out_f, SEP_WRITE, 0);
    model = sep_open(model_f, SEP_READ, 0);

    if (!se_have_par("verb")) verb = 1; else verb = se_get_par_int("verb");
    if (!se_have_par("cmp")) cmp = 1; else cmp = se_get_par_int("cmp");
    if (cmp != 0 && cmp != 1) {
        ERROR(("cmp must be 0 or 1. cmp=1: relative receiver coordinate; cmp=0: absolute receiver coordinate."));
    }
    INFO(("Read cmp=%d (%s receiver coordinates).", cmp, cmp ? "relative" : "absolute"));

    if (sep_have_hdr_int(model, "should_datum")) {
        should_datum = sep_get_hdr_int(model, "should_datum", 0);
    } else {
        ERROR(("Need should_datum in model_file header"));
    }
    INFO(("Read should_datum=%d from %s", should_datum, model_f));

    if (!se_have_par("aperture")) aper = 50; else aper = se_get_par_int("aperture");
    if (!se_have_par("taper")) tap = 10; else tap = se_get_par_int("taper");
    if (!se_have_par("length")) length = 0.025f; else length = se_get_par_float("length");

    /* Read input. */
    if (in->headers->ndim < 3) ERROR(("Input must be 3D."));
    nt = in->headers->n[0];
    nh = in->headers->n[1];
    ns = in->headers->n[2];
    dt = in->headers->d[0];
    h0 = in->headers->o[1];
    dh = in->headers->d[1];
    s0 = in->headers->o[2];
    ds = in->headers->d[2];

    if (!se_have_par("fmin")) fmin = 0.0f; else fmin = se_get_par_float("fmin");
    if (!se_have_par("fmax")) fmax = 0.5f / dt; else fmax = se_get_par_float("fmax");

    sdatum = (float)(model->headers->o[0] + should_datum * model->headers->d[0]);
    rdatum = sdatum;

    tr_in = alloc3float(nt, nh, ns);
    se_fsio_read_float(in->data->io, tr_in[0][0], nt * nh * ns);

    reverse_trace(nt, nh, ns, tr_in);

    tr_out = alloc3float(nt, nh, ns);
    memset(tr_out[0][0], 0, sizeof(float) * nt * nh * ns);

    /* Read source Green's function. */
    sgreen = sep_open(sgreen_f, SEP_READ, 0);
    if (sgreen->headers->ndim < 2) ERROR(("Source Green's function must be 2D."));
    nsg = sgreen->headers->n[0];
    sg0 = sgreen->headers->o[0];
    dsg = sgreen->headers->d[0];

    stable = alloc2float(nsg, nsg);
    se_fsio_read_float(sgreen->data->io, stable[0], nsg * nsg);
    sep_close(sgreen);
    sgreen = NULL;

    /* Read receiver Green's function. */
    rgreen = sep_open(rgreen_f, SEP_READ, 0);
    if (rgreen->headers->ndim < 2) ERROR(("Receiver Green's function must be 2D."));
    nrg = rgreen->headers->n[0];
    rg0 = rgreen->headers->o[0];
    drg = rgreen->headers->d[0];

    rtable = alloc2float(nrg, nrg);
    se_fsio_read_float(rgreen->data->io, rtable[0], nrg * nrg);
    sep_close(rgreen);
    rgreen = NULL;

    if (!se_have_par("interm")) interm_f = NULL; else interm_f = se_get_par_str("interm");
    if (interm_f != NULL) {
        interm = sep_open(interm_f, SEP_WRITE, 0);
    }

    /* Initialize direct frequency-integral summation. */
    kirdat_integral_init(nt, dt, length, fmin, fmax);

    /* ------------------------------------------------------------------ */
    /* Receiver-side continuation: common-shot gather.                     */
    /* ------------------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel for private(ih,c,left,right,ic,cc,coef,tau,dist,it0,delta,amp,it)
#endif
    for (is = 0; is < ns; is++) {
        if (verb) INFO(("Processing common-shot gather %d of %d.", is + 1, ns));

        for (ih = 0; ih < nh; ih++) {
            c = (int)(((cmp ? (s0 + is * ds) : 0.0f) + h0 + ih * dh - rg0) / drg + 0.5f);
            if (c < 0 || c > nrg - 1) ERROR(("Receiver table too small."));

            left  = (ih - aper < 0) ? 0 : ih - aper;
            right = (ih + aper > nh - 1) ? nh - 1 : ih + aper;

            for (ic = left; ic <= right; ic++) {
                cc = (int)(((cmp ? (s0 + is * ds) : 0.0f) + h0 + ic * dh - rg0) / drg + 0.5f);
                if (cc < 0 || cc > nrg - 1) ERROR(("Receiver table too small."));

                coef = taper_weight(ic, left, right, tap, 1);

                tau = rtable[cc][c];
                if (tau <= 0.0f) continue;

                dist = rdatum * rdatum + (ic - ih) * dh * (ic - ih) * dh;
                if (dist <= 0.0f) continue;

                amp = coef / (float)M_PI * dh * rdatum * tau / dist;

                it0 = first_output_index(tau, dt, nt, &delta);
                if (it0 >= nt) continue;

                const float *filtered = kirdat_filter_trace_integral((const float *)tr_in[is][ic], tau, delta);

                for (it = it0; it < nt; it++) {
                    const int shift = it - it0;
                    tr_out[is][ih][it] += amp * filtered[shift];
                }
            }
        }
    }

    if (interm_f != NULL) {
        se_fsio_write_float(interm->data->io, tr_out[0][0], nt * nh * ns);
        sep_close(interm);
        interm = NULL;
    }

    memset(tr_in[0][0], 0, sizeof(float) * nt * nh * ns);

    /* ------------------------------------------------------------------ */
    /* Source-side continuation: common-receiver gather.                   */
    /* ------------------------------------------------------------------ */
    if (cmp == 1) {
        /* Original relative-coordinate geometry. */
        s = fabsf((ns - 1) * ds);
        h = fabsf((nh - 1) * dh);

        if (fabsf(ds) >= fabsf(dh)) {
            dr = fabsf(dh);
            jump = 1;
        } else {
            dr = fabsf(ds);
            jump = (int)(dh / ds + 0.5f);
            if (jump < 1) jump = 1;
        }

        nr = (int)(s / dr + h / dr + 1.5f);

#ifdef _OPENMP
#pragma omp parallel for private(r,sleft,sright,is,c,ih,left,right,ic,cc,hh,coef,tau,dist,it0,delta,amp,it)
#endif
        for (ir = 0; ir < nr; ir++) {
            if (verb) INFO(("Processing common-receiver gather %d of %d.", ir + 1, nr));

            r = ir * dr + ((ds <= 0.0f) ? -1.0f : 0.0f) * s + ((dh <= 0.0f) ? -1.0f : 0.0f) * h;

            sleft  = (int)((ir * dr + ((ds <= 0.0f) ? -1.0f : 0.0f) * s + ((ds <= 0.0f) ?  0.0f : -1.0f) * h) / ds + 0.5f);
            sright = (int)((ir * dr + ((ds <= 0.0f) ? -1.0f : 0.0f) * s + ((ds <= 0.0f) ? -1.0f :  0.0f) * h) / ds + 0.5f);

            if (sleft < 0) sleft = 0;
            if (sright > ns - 1) sright = ns - 1;

            left = (int)((r - sleft * ds) / dh + 0.5f);
            if (left < 0 || left > nh - 1) sleft++;

            right = (int)((r - sright * ds) / dh + 0.5f);
            if (right < 0 || right > nh - 1) sright--;

            for (is = sleft; is <= sright; is += jump) {
                c = (int)((s0 + is * ds - sg0) / dsg + 0.5f);
                if (c < 0 || c > nsg - 1) ERROR(("Source table too small."));

                ih = (int)((r - is * ds) / dh + 0.5f);
                if (ih < 0 || ih > nh - 1) continue;

                left  = (is - jump * aper < sleft)  ? sleft  : is - jump * aper;
                right = (is + jump * aper > sright) ? sright : is + jump * aper;

                for (ic = left; ic <= right; ic += jump) {
                    cc = (int)((s0 + ic * ds - sg0) / dsg + 0.5f);
                    if (cc < 0 || cc > nsg - 1) ERROR(("Source table too small."));

                    hh = (int)((r - ic * ds) / dh + 0.5f);
                    if (hh < 0 || hh > nh - 1) continue;

                    coef = taper_weight(ic, left, right, tap, jump);

                    tau = stable[cc][c];
                    if (tau <= 0.0f) continue;

                    dist = sdatum * sdatum + (ic - is) * ds * (ic - is) * ds;
                    if (dist <= 0.0f) continue;

                    amp = coef / (float)M_PI * ds * sdatum * tau / dist;

                    it0 = first_output_index(tau, dt, nt, &delta);
                    if (it0 >= nt) continue;

                    const float *filtered = kirdat_filter_trace_integral((const float *)tr_out[ic][hh], tau, delta);

                    for (it = it0; it < nt; it++) {
                        const int shift = it - it0;
                        tr_in[is][ih][it] += amp * filtered[shift];
                    }
                }
            }
        }
    } else {
        /* Absolute receiver-coordinate geometry: fixed ih is already a common-receiver gather. */
#ifdef _OPENMP
#pragma omp parallel for private(is,c,left,right,ic,cc,coef,tau,dist,it0,delta,amp,it)
#endif
        for (ih = 0; ih < nh; ih++) {
            if (verb) INFO(("Processing absolute-coordinate receiver gather %d of %d.", ih + 1, nh));

            for (is = 0; is < ns; is++) {
                c = (int)((s0 + is * ds - sg0) / dsg + 0.5f);
                if (c < 0 || c > nsg - 1) ERROR(("Source table too small."));

                left  = (is - aper < 0) ? 0 : is - aper;
                right = (is + aper > ns - 1) ? ns - 1 : is + aper;

                for (ic = left; ic <= right; ic++) {
                    cc = (int)((s0 + ic * ds - sg0) / dsg + 0.5f);
                    if (cc < 0 || cc > nsg - 1) ERROR(("Source table too small."));

                    coef = taper_weight(ic, left, right, tap, 1);

                    tau = stable[cc][c];
                    if (tau <= 0.0f) continue;

                    dist = sdatum * sdatum + (ic - is) * ds * (ic - is) * ds;
                    if (dist <= 0.0f) continue;

                    amp = coef / (float)M_PI * ds * sdatum * tau / dist;

                    it0 = first_output_index(tau, dt, nt, &delta);
                    if (it0 >= nt) continue;

                    const float *filtered = kirdat_filter_trace_integral((const float *)tr_out[ic][ih], tau, delta);

                    for (it = it0; it < nt; it++) {
                        const int shift = it - it0;
                        tr_in[is][ih][it] += amp * filtered[shift];
                    }
                }
            }
        }
    }

    out->headers->ndim = 3;
    out->headers->n[0] = nt;
    out->headers->n[1] = nh;
    out->headers->n[2] = ns;
    out->headers->d[0] = dt;
    out->headers->d[1] = dh;
    out->headers->d[2] = ds;
    out->headers->o[0] = 0.0f;
    out->headers->o[1] = h0;
    out->headers->o[2] = s0;

    reverse_trace(nt, nh, ns, tr_in);

    se_fsio_write_float(out->data->io, tr_in[0][0], nt * nh * ns);

    kirdat_integral_close();

    if (in != NULL) sep_close(in);
    if (out != NULL) sep_close(out);
    if (model != NULL) sep_close(model);

    if (tr_in != NULL) free3float(tr_in);
    if (tr_out != NULL) free3float(tr_out);
    if (stable != NULL) free2float(stable);
    if (rtable != NULL) free2float(rtable);

    INFO(("Done."));
    return 0;
}
