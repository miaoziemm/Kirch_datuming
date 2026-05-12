#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <fftw3.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static inline size_t freq_trace_offset(int is, int ih, int nh, int nf)
{
    return ((size_t)is * (size_t)nh + (size_t)ih) * (size_t)nf;
}


static inline fftwf_complex* freq_trace(fftwf_complex *f, int is, int ih, int nh, int nf)
{
    return f + freq_trace_offset(is, ih, nh, nf);
}


static inline const fftwf_complex* freq_trace_const(const fftwf_complex *f, int is, int ih, int nh, int nf)
{
    return f + freq_trace_offset(is, ih, nh, nf);
}


static inline int next_pow2(int n)
{
    int m = 1;
    while (m < n) m <<= 1;
    return m;
}


void reverse_trace(int nt, int nh, int ns, float ***tr)
{
    INFO(("Reversing traces..."));

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int is = 0; is < ns; is++) {
        for (int ih = 0; ih < nh; ih++) {
            for (int it = 0; it < nt/2; it++) {
                float tmp = tr[is][ih][it];
                tr[is][ih][it] = tr[is][ih][nt-1-it];
                tr[is][ih][nt-1-it] = tmp;
            }
        }
    }
}


static void time_to_freq_all(
    float ***tr,
    fftwf_complex *ftr,
    int nt,
    int nh,
    int ns,
    int nfft,
    int nf
)
{
    float *buf_t = (float*)fftwf_malloc(sizeof(float) * (size_t)nfft);
    fftwf_complex *buf_f = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)nf);

    if (buf_t == NULL || buf_f == NULL) {
        ERROR(("FFTW buffer allocation failed in time_to_freq_all."));
    }

    fftwf_plan plan_f = fftwf_plan_dft_r2c_1d(nfft, buf_t, buf_f, FFTW_ESTIMATE);
    if (plan_f == NULL) {
        ERROR(("FFTW forward plan creation failed."));
    }

    for (int is = 0; is < ns; is++) {
        for (int ih = 0; ih < nh; ih++) {

            memset(buf_t, 0, sizeof(float) * (size_t)nfft);

            for (int it = 0; it < nt; it++) {
                buf_t[it] = tr[is][ih][it];
            }

            fftwf_execute(plan_f);

            fftwf_complex *dst = freq_trace(ftr, is, ih, nh, nf);
            memcpy(dst, buf_f, sizeof(fftwf_complex) * (size_t)nf);
        }
    }

    fftwf_destroy_plan(plan_f);
    fftwf_free(buf_t);
    fftwf_free(buf_f);
}


static void freq_to_time_all(
    const fftwf_complex *ftr,
    float ***tr,
    int nt,
    int nh,
    int ns,
    int nfft,
    int nf
)
{
    float *buf_t = (float*)fftwf_malloc(sizeof(float) * (size_t)nfft);
    fftwf_complex *buf_f = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)nf);

    if (buf_t == NULL || buf_f == NULL) {
        ERROR(("FFTW buffer allocation failed in freq_to_time_all."));
    }

    fftwf_plan plan_b = fftwf_plan_dft_c2r_1d(nfft, buf_f, buf_t, FFTW_ESTIMATE);
    if (plan_b == NULL) {
        ERROR(("FFTW backward plan creation failed."));
    }

    for (int is = 0; is < ns; is++) {
        for (int ih = 0; ih < nh; ih++) {

            const fftwf_complex *src = freq_trace_const(ftr, is, ih, nh, nf);
            memcpy(buf_f, src, sizeof(fftwf_complex) * (size_t)nf);

            fftwf_execute(plan_b);

            for (int it = 0; it < nt; it++) {
                tr[is][ih][it] = buf_t[it] / (float)nfft;
            }
        }
    }

    fftwf_destroy_plan(plan_b);
    fftwf_free(buf_t);
    fftwf_free(buf_f);
}


static inline void add_delayed_trace_freq(
    fftwf_complex *out,
    const fftwf_complex *in,
    const std::vector<double> &omega,
    int iw_min,
    int iw_max,
    double amp,
    double tau
)
{
    if (amp == 0.0) return;

    for (int iw = iw_min; iw <= iw_max; iw++) {

        double phase = -omega[iw] * tau;

        double cp = cos(phase);
        double sp = sin(phase);

        double xr = (double)in[iw][0];
        double xi = (double)in[iw][1];

        out[iw][0] = (float)(out[iw][0] + amp * (cp * xr - sp * xi));
        out[iw][1] = (float)(out[iw][1] + amp * (sp * xr + cp * xi));
    }
}


int main(int argc, char* argv[])
{
    se_par_init(argc, argv);

    int verb;
    int nt, nh, ns;
    int nsg, nrg;
    int aper, tap;
    int nr, jump;

    float sdatum = 0.0f;
    float rdatum = 0.0f;
    float length;
    float dt, h0, dh, s0, ds;
    float sg0, dsg, rg0, drg;

    float ***tr_in = NULL;
    float ***tr_out = NULL;
    float **stable = NULL;
    float **rtable = NULL;

    sep_t *in = NULL;
    sep_t *out = NULL;
    sep_t *sgreen = NULL;
    sep_t *rgreen = NULL;
    sep_t *interm = NULL;

    char *in_f = NULL;
    char *out_f = NULL;
    char *sgreen_f = NULL;
    char *rgreen_f = NULL;
    char *interm_f = NULL;

    if (!se_have_par("input_file"))  ERROR(("Need input_file="));
    else in_f = se_get_par_str("input_file");

    if (!se_have_par("output_file")) ERROR(("Need output_file="));
    else out_f = se_get_par_str("output_file");

    if (!se_have_par("sgreen_file")) ERROR(("Need sgreen_file="));
    else sgreen_f = se_get_par_str("sgreen_file");

    if (!se_have_par("rgreen_file")) ERROR(("Need rgreen_file="));
    else rgreen_f = se_get_par_str("rgreen_file");

    in = sep_open(in_f, SEP_READ, 0);
    out = sep_open(out_f, SEP_WRITE, 0);

    if (!se_have_par("verb")) verb = 1;
    else verb = se_get_par_int("verb");

    if (!se_have_par("sdatum")) ERROR(("Need sdatum="));
    else sdatum = se_get_par_float("sdatum");

    if (!se_have_par("rdatum")) ERROR(("Need rdatum="));
    else rdatum = se_get_par_float("rdatum");

    if (!se_have_par("aperture")) aper = 50;
    else aper = se_get_par_int("aperture");

    if (!se_have_par("taper")) tap = 10;
    else tap = se_get_par_int("taper");

    if (!se_have_par("length")) length = 0.025f;
    else length = se_get_par_float("length");

    /*
     * Optional frequency band.
     * If not provided, all positive frequencies are used.
     */
    float fmin = 0.0f;
    float fmax = -1.0f;

    if (se_have_par("fmin")) fmin = se_get_par_float("fmin");
    if (se_have_par("fmax")) fmax = se_get_par_float("fmax");

    /* read input */
    if (in->headers->ndim < 3) ERROR(("Input must be 3D."));

    nt = in->headers->n[0];
    nh = in->headers->n[1];
    ns = in->headers->n[2];

    dt = in->headers->d[0];
    h0 = in->headers->o[1];
    dh = in->headers->d[1];
    s0 = in->headers->o[2];
    ds = in->headers->d[2];

    if (fmax <= 0.0f) fmax = 0.5f / dt;

    if (verb) {
        INFO(("Input dimensions: nt=%d, nh=%d, ns=%d", nt, nh, ns));
        INFO(("dt=%g, dh=%g, ds=%g", dt, dh, ds));
        INFO(("Frequency-domain implementation: fmin=%g Hz, fmax=%g Hz", fmin, fmax));
        INFO(("The original length=%g parameter is kept for compatibility but is not used by this pure phase-delay version.", length));
    }

    tr_in = alloc3float(nt, nh, ns);
    se_fsio_read_float(in->data->io, tr_in[0][0], nt * nh * ns);

    reverse_trace(nt, nh, ns, tr_in);

    /* read Green's function: source table */
    sgreen = sep_open(sgreen_f, SEP_READ, 0);

    if (sgreen->headers->ndim < 2) {
        ERROR(("Source Green's function must be 2D."));
    }

    nsg = sgreen->headers->n[0];
    sg0 = sgreen->headers->o[0];
    dsg = sgreen->headers->d[0];

    stable = alloc2float(nsg, nsg);
    se_fsio_read_float(sgreen->data->io, stable[0], nsg * nsg);
    sep_close(sgreen);

    /* read Green's function: receiver table */
    rgreen = sep_open(rgreen_f, SEP_READ, 0);

    if (rgreen->headers->ndim < 2) {
        ERROR(("Receiver Green's function must be 2D."));
    }

    nrg = rgreen->headers->n[0];
    rg0 = rgreen->headers->o[0];
    drg = rgreen->headers->d[0];

    rtable = alloc2float(nrg, nrg);
    se_fsio_read_float(rgreen->data->io, rtable[0], nrg * nrg);
    sep_close(rgreen);

    /* intermediate output */
    if (!se_have_par("interm")) {
        interm_f = NULL;
    } else {
        interm_f = se_get_par_str("interm");
    }

    if (interm_f != NULL) {
        interm = sep_open(interm_f, SEP_WRITE, 0);

        interm->headers->ndim = 3;
        interm->headers->n[0] = nt;
        interm->headers->n[1] = nh;
        interm->headers->n[2] = ns;
        interm->headers->d[0] = dt;
        interm->headers->d[1] = dh;
        interm->headers->d[2] = ds;
        interm->headers->o[0] = 0.0f;
        interm->headers->o[1] = h0;
        interm->headers->o[2] = s0;
    }

    /*
     * FFT setup.
     * Use zero padding to reduce circular-shift wraparound.
     */
    int nfft = next_pow2(2 * nt);
    int nf = nfft / 2 + 1;

    if (verb) {
        INFO(("Using nfft=%d, nf=%d", nfft, nf));
    }

    std::vector<double> omega(nf);
    for (int iw = 0; iw < nf; iw++) {
        double freq = (double)iw / ((double)nfft * (double)dt);
        omega[iw] = 2.0 * M_PI * freq;
    }

    int iw_min = (int)ceil((double)fmin * (double)nfft * (double)dt);
    int iw_max = (int)floor((double)fmax * (double)nfft * (double)dt);

    if (iw_min < 0) iw_min = 0;
    if (iw_max > nf - 1) iw_max = nf - 1;
    if (iw_min > iw_max) {
        ERROR(("Invalid frequency band: iw_min > iw_max."));
    }

    size_t ntrace = (size_t)ns * (size_t)nh;
    size_t nfreq_total = ntrace * (size_t)nf;

    fftwf_complex *f_in  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * nfreq_total);
    fftwf_complex *f_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * nfreq_total);

    if (f_in == NULL || f_out == NULL) {
        ERROR(("FFTW frequency-domain array allocation failed."));
    }

    memset(f_in,  0, sizeof(fftwf_complex) * nfreq_total);
    memset(f_out, 0, sizeof(fftwf_complex) * nfreq_total);

    /*
     * Convert input traces to frequency domain.
     */
    if (verb) INFO(("Transforming input traces to frequency domain..."));
    time_to_freq_all(tr_in, f_in, nt, nh, ns, nfft, nf);

    /*
     * ------------------------------------------------------------------
     * 1. Common-shot gather receiver-side datuming in frequency domain.
     *
     * Time-domain original:
     * tr_out[is][ih][it] += W_r * kirdat_pick(..., tr_in[is][ic], ...)
     *
     * Frequency-domain version:
     * F_OUT[is][ih](w) += W_r * exp(-i*w*tau_r) * F_IN[is][ic](w)
     * ------------------------------------------------------------------
     */

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int is = 0; is < ns; is++) {

        if (verb) INFO(("Processing common-shot gather %d of %d.", is + 1, ns));

        for (int ih = 0; ih < nh; ih++) {

            int c = (int)((s0 + is * ds + h0 + ih * dh - rg0) / drg + 0.5f);
            if (c < 0 || c > nrg - 1) {
                ERROR(("Receiver table too small."));
            }

            int left  = (ih - aper < 0)      ? 0      : ih - aper;
            int right = (ih + aper > nh - 1) ? nh - 1 : ih + aper;

            fftwf_complex *out_trace = freq_trace(f_out, is, ih, nh, nf);

            for (int ic = left; ic <= right; ic++) {

                int cc = (int)((s0 + is * ds + h0 + ic * dh - rg0) / drg + 0.5f);
                if (cc < 0 || cc > nrg - 1) {
                    ERROR(("Receiver table too small."));
                }

                float coef = 1.0f;

                if (tap > 0) {
                    coef *= (ic - left  >= tap) ? 1.0f : (float)(ic - left) / (float)tap;
                    coef *= (right - ic >= tap) ? 1.0f : (float)(right - ic) / (float)tap;
                }

                float tau = rtable[cc][c];

                float dxr = (float)(ic - ih) * dh;
                float dist = rdatum * rdatum + dxr * dxr;

                if (dist <= 0.0f) continue;
                if (tau < 0.0f) continue;

                double amp = (double)coef / M_PI
                           * (double)dh
                           * (double)rdatum
                           * (double)tau
                           / (double)dist;

                const fftwf_complex *in_trace = freq_trace_const(f_in, is, ic, nh, nf);

                add_delayed_trace_freq(
                    out_trace,
                    in_trace,
                    omega,
                    iw_min,
                    iw_max,
                    amp,
                    (double)tau
                );
            }
        }
    }

    /*
     * If requested, write intermediate receiver-datumed data in time domain.
     */
    if (interm_f != NULL) {
        if (verb) INFO(("Writing intermediate result..."));

        tr_out = alloc3float(nt, nh, ns);

        freq_to_time_all(f_out, tr_out, nt, nh, ns, nfft, nf);
        se_fsio_write_float(interm->data->io, tr_out[0][0], nt * nh * ns);

        free3float(tr_out);
        tr_out = NULL;
    }

    /*
     * Reuse f_in as final frequency-domain output.
     */
    memset(f_in, 0, sizeof(fftwf_complex) * nfreq_total);

    /*
     * Acquisition geometry for common-receiver gather.
     */
    float s = fabsf((ns - 1) * ds);
    float h = fabsf((nh - 1) * dh);

    float dr;

    if (fabsf(ds) >= fabsf(dh)) {
        dr = fabsf(dh);
        jump = 1;
    } else {
        dr = fabsf(ds);
        jump = (int)(dh / ds + 0.5f);
        if (jump < 1) jump = 1;
    }

    nr = (int)((s + h) / dr + 1.5f);

    /*
     * ------------------------------------------------------------------
     * 2. Common-receiver gather source-side datuming in frequency domain.
     *
     * Time-domain original:
     * tr_in[is][ih][it] += W_s * kirdat_pick(..., tr_out[ic][hh], ...)
     *
     * Frequency-domain version:
     * F_IN[is][ih](w) += W_s * exp(-i*w*tau_s) * F_OUT[ic][hh](w)
     * ------------------------------------------------------------------
     */

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int ir = 0; ir < nr; ir++) {

        if (verb) INFO(("Processing common-receiver gather %d of %d.", ir + 1, nr));

        float r = ir * dr
                + ((ds <= 0.0f) ? -1.0f : 0.0f) * s
                + ((dh <= 0.0f) ? -1.0f : 0.0f) * h;

        /*
         * Source-receiver reciprocity.
         */
        int sleft = (int)((ir * dr
                 + ((ds <= 0.0f) ? -1.0f : 0.0f) * s
                 + ((ds <= 0.0f) ?  0.0f : -1.0f) * h) / ds + 0.5f);

        int sright = (int)((ir * dr
                  + ((ds <= 0.0f) ? -1.0f : 0.0f) * s
                  + ((ds <= 0.0f) ? -1.0f : 0.0f) * h) / ds + 0.5f);

        if (sleft < 0) sleft = 0;
        if (sright > ns - 1) sright = ns - 1;

        int left_check = (int)((r - sleft * ds) / dh + 0.5f);
        if (left_check < 0 || left_check > nh - 1) sleft++;

        int right_check = (int)((r - sright * ds) / dh + 0.5f);
        if (right_check < 0 || right_check > nh - 1) sright--;

        for (int is = sleft; is <= sright; is += jump) {

            int c = (int)((s0 + is * ds - sg0) / dsg + 0.5f);
            if (c < 0 || c > nsg - 1) {
                ERROR(("Source table too small."));
            }

            int ih = (int)((r - is * ds) / dh + 0.5f);
            if (ih < 0 || ih > nh - 1) continue;

            int left  = (is - jump * aper < sleft)  ? sleft  : is - jump * aper;
            int right = (is + jump * aper > sright) ? sright : is + jump * aper;

            fftwf_complex *out_trace = freq_trace(f_in, is, ih, nh, nf);

            for (int ic = left; ic <= right; ic += jump) {

                int cc = (int)((s0 + ic * ds - sg0) / dsg + 0.5f);
                if (cc < 0 || cc > nsg - 1) {
                    ERROR(("Source table too small."));
                }

                int hh = (int)((r - ic * ds) / dh + 0.5f);
                if (hh < 0 || hh > nh - 1) continue;

                float coef = 1.0f;

                if (tap > 0) {
                    coef *= (ic - left  >= tap) ? 1.0f : (float)(ic - left) / (float)jump / (float)tap;
                    coef *= (right - ic >= tap) ? 1.0f : (float)(right - ic) / (float)jump / (float)tap;
                }

                float tau = stable[cc][c];

                float dxs = (float)(ic - is) * ds;
                float dist = sdatum * sdatum + dxs * dxs;

                if (dist <= 0.0f) continue;
                if (tau < 0.0f) continue;

                double amp = (double)coef / M_PI
                           * (double)ds
                           * (double)sdatum
                           * (double)tau
                           / (double)dist;

                const fftwf_complex *in_trace = freq_trace_const(f_out, ic, hh, nh, nf);

                add_delayed_trace_freq(
                    out_trace,
                    in_trace,
                    omega,
                    iw_min,
                    iw_max,
                    amp,
                    (double)tau
                );
            }
        }
    }

    /*
     * Convert final frequency-domain result back to time domain.
     */
    if (verb) INFO(("Transforming final result back to time domain..."));
    freq_to_time_all(f_in, tr_in, nt, nh, ns, nfft, nf);

    /*
     * Output header.
     */
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

    /*
     * Write output.
     */
    se_fsio_write_float(out->data->io, tr_in[0][0], nt * nh * ns);

    /*
     * Clean up.
     */
    sep_close(in);
    sep_close(out);

    if (interm != NULL) {
        sep_close(interm);
    }

    free3float(tr_in);
    free2float(stable);
    free2float(rtable);

    fftwf_free(f_in);
    fftwf_free(f_out);

    INFO(("Done."));

    return 0;
}