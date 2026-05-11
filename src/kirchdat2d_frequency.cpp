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
    while (p < n) {
        p <<= 1;
    }
    return p;
}

static float max_table_value(float **table, int n)
{
    float maxv = 0.0f;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float v = table[i][j];
            if (v > maxv) {
                maxv = v;
            }
        }
    }
    return maxv;
}

static void build_kirdat_filter(float tau, float dt, int nsam, std::vector<float> &filt)
{
    filt.assign(nsam, 0.0f);
    if (tau <= 0.0f) {
        return;
    }

    for (int isam = 1; isam < nsam; isam++) {
        float ratio = (tau + isam * dt) / tau;
        filt[isam] = sqrtf(ratio * ratio - 1.0f);
    }

    for (int isam = 0; isam < nsam - 1; isam++) {
        filt[isam] = filt[isam + 1] - filt[isam];
    }

    for (int isam = nsam - 2; isam > 0; isam--) {
        filt[isam] = filt[isam] - filt[isam - 1];
    }

    for (int isam = 0; isam < nsam - 1; isam++) {
        filt[isam] /= dt;
    }
}

static int compute_it0(float tau, float dt, float *delta)
{
    int it0 = (int)ceilf(tau / dt);
    if (it0 < 0) {
        it0 = 0;
    }
    if (delta) {
        *delta = ((float)it0 * dt - tau) / dt;
    }
    return it0;
}

static void fill_fft_input(float *dst, float ***src, int nt, int nh, int ns, int nfft)
{
    int ntraces = ns * nh;
    for (int itr = 0; itr < ntraces; itr++) {
        int is = itr / nh;
        int ih = itr % nh;
        float *trace = src[is][ih];
        float *buf = dst + (size_t)itr * nfft;
        memcpy(buf, trace, (size_t)nt * sizeof(float));
        memset(buf + nt, 0, (size_t)(nfft - nt) * sizeof(float));
    }
}

static void apply_ifft_output(float ***dst, const float *src, const float *corr,
                              int nt, int nh, int ns, int nfft, float scale)
{
    int ntraces = ns * nh;
    for (int itr = 0; itr < ntraces; itr++) {
        int is = itr / nh;
        int ih = itr % nh;
        const float *buf = src + (size_t)itr * nfft;
        const float *corr_trace = corr + (size_t)itr * nt;
        for (int it = 0; it < nt; it++) {
            dst[is][ih][it] = buf[it] * scale + corr_trace[it];
        }
    }
}

int main(int argc, char *argv[])
{
    se_par_init(argc, argv);
    int verb;
    int nt, ih, nh, is, ns, nsg, nrg, left, right, ic, aper, c, cc, hh;
    int ir, nr, jump, sleft, sright, tap;
    float sdatum = 0.0f;
    float rdatum = 0.0f;
    float length, dt, h0, dh, s0, ds, sg0, dsg, rg0, drg, dist, tau;
    float r, dr, s, h, coef;
    float ***tr_in, ***tr_out, **stable, **rtable;
    sep_t *in, *out, *sgreen, *rgreen, *interm = NULL;
    char *in_f = NULL, *out_f = NULL, *sgreen_f = NULL, *rgreen_f = NULL, *interm_f = NULL;

    if (!se_have_par("input_file")) ERROR(("Need input_file=")); else in_f = se_get_par_str("input_file");
    if (!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");
    if (!se_have_par("sgreen_file")) ERROR(("Need sgreen_file=")); else sgreen_f = se_get_par_str("sgreen_file");
    if (!se_have_par("rgreen_file")) ERROR(("Need rgreen_file=")); else rgreen_f = se_get_par_str("rgreen_file");

    in = sep_open(in_f, SEP_READ, 0);
    out = sep_open(out_f, SEP_WRITE, 0);

    if (!se_have_par("verb")) verb = 1; else verb = se_get_par_int("verb");
    if (!se_have_par("sdatum")) ERROR(("Need sdatum=")); else sdatum = se_get_par_float("sdatum");
    if (!se_have_par("rdatum")) ERROR(("Need rdatum=")); else rdatum = se_get_par_float("rdatum");

    if (!se_have_par("aperture")) aper = 50; else aper = se_get_par_int("aperture");
    if (!se_have_par("taper")) tap = 10; else tap = se_get_par_int("taper");
    if (!se_have_par("length")) length = 0.025f; else length = se_get_par_float("length");

    if (in->headers->ndim < 3) ERROR(("Input must be 3D."));
    nt = in->headers->n[0];
    nh = in->headers->n[1];
    ns = in->headers->n[2];
    dt = in->headers->d[0];
    h0 = in->headers->o[1];
    dh = in->headers->d[1];
    s0 = in->headers->o[2];
    ds = in->headers->d[2];

    tr_in = alloc3float(nt, nh, ns);
    se_fsio_read_float(in->data->io, tr_in[0][0], nt * nh * ns);

    reverse_trace(nt, nh, ns, tr_in);

    tr_out = alloc3float(nt, nh, ns);

    sgreen = sep_open(sgreen_f, SEP_READ, 0);
    if (sgreen->headers->ndim < 2) ERROR(("Source Green's function must be 2D."));
    nsg = sgreen->headers->n[0];
    sg0 = sgreen->headers->o[0];
    dsg = sgreen->headers->d[0];

    stable = alloc2float(nsg, nsg);
    se_fsio_read_float(sgreen->data->io, stable[0], nsg * nsg);
    sep_close(sgreen);

    rgreen = sep_open(rgreen_f, SEP_READ, 0);
    if (rgreen->headers->ndim < 2) ERROR(("Receiver Green's function must be 2D."));
    nrg = rgreen->headers->n[0];
    rg0 = rgreen->headers->o[0];
    drg = rgreen->headers->d[0];

    rtable = alloc2float(nrg, nrg);
    se_fsio_read_float(rgreen->data->io, rtable[0], nrg * nrg);
    sep_close(rgreen);

    if (!se_have_par("interm")) interm_f = NULL; else interm_f = se_get_par_str("interm");
    if (interm_f != NULL) {
        interm = sep_open(interm_f, SEP_WRITE, 0);
    }

    float max_tau = max_table_value(stable, nsg);
    float max_tau_r = max_table_value(rtable, nrg);
    if (max_tau_r > max_tau) {
        max_tau = max_tau_r;
    }

    int nsam = (int)(length / dt) + 2;
    int nfilt = nsam - 1;
    int it0_max = compute_it0(max_tau, dt, NULL);
    int nfft_min = nt + nfilt - 1 + it0_max;
    int nfft = next_pow2(nfft_min);
    int nf = nfft / 2 + 1;
    int ntraces = ns * nh;

    float *fft_time = (float *)fftwf_malloc(sizeof(float) * (size_t)ntraces * nfft);
    float *fft_time_out = (float *)fftwf_malloc(sizeof(float) * (size_t)ntraces * nfft);
    fftwf_complex *fft_freq = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (size_t)ntraces * nf);
    fftwf_complex *fft_freq_accum = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (size_t)ntraces * nf);
    float *corr = alloc1float((size_t)ntraces * nt);

    int max_threads = 1;
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#endif

    std::vector<float *> filt_time_thread(max_threads, NULL);
    std::vector<fftwf_complex *> filt_freq_thread(max_threads, NULL);
    std::vector<fftwf_plan> plan_filter_thread(max_threads, NULL);

    for (int t = 0; t < max_threads; t++) {
        filt_time_thread[t] = (float *)fftwf_malloc(sizeof(float) * (size_t)nfft);
        filt_freq_thread[t] = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (size_t)nf);
        if (filt_time_thread[t] && filt_freq_thread[t]) {
            plan_filter_thread[t] = fftwf_plan_dft_r2c_1d(nfft, filt_time_thread[t], filt_freq_thread[t], FFTW_ESTIMATE);
        }
    }

    if (!fft_time || !fft_time_out || !fft_freq || !fft_freq_accum || !corr) {
        ERROR(("FFTW memory allocation failed."));
    }

    for (int t = 0; t < max_threads; t++) {
        if (!filt_time_thread[t] || !filt_freq_thread[t] || !plan_filter_thread[t]) {
            ERROR(("FFTW plan creation failed."));
        }
    }

    int n[1] = { nfft };
    fftwf_plan plan_forward = fftwf_plan_many_dft_r2c(1, n, ntraces,
                                                      fft_time, NULL, 1, nfft,
                                                      fft_freq, NULL, 1, nf,
                                                      FFTW_ESTIMATE);
    fftwf_plan plan_inverse = fftwf_plan_many_dft_c2r(1, n, ntraces,
                                                      fft_freq_accum, NULL, 1, nf,
                                                      fft_time_out, NULL, 1, nfft,
                                                      FFTW_ESTIMATE);
    if (!plan_forward || !plan_inverse) {
        ERROR(("FFTW plan creation failed."));
    }

    std::vector<float> cos_wdt(nf);
    std::vector<float> sin_wdt(nf);
    for (int k = 0; k < nf; k++) {
        float wdt = 2.0f * (float)M_PI * (float)k / (float)nfft;
        cos_wdt[k] = cosf(wdt);
        sin_wdt[k] = sinf(wdt);
    }

    fill_fft_input(fft_time, tr_in, nt, nh, ns, nfft);
    fftwf_execute(plan_forward);

    memset(fft_freq_accum, 0, sizeof(fftwf_complex) * (size_t)ntraces * nf);
    memset(corr, 0, sizeof(float) * (size_t)ntraces * nt);

    std::vector<std::vector<float>> filt_thread(max_threads, std::vector<float>(nsam, 0.0f));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (is = 0; is < ns; is++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        float *filt_time = filt_time_thread[tid];
        fftwf_complex *filt_freq = filt_freq_thread[tid];
        fftwf_plan plan_filter = plan_filter_thread[tid];
        std::vector<float> &filt = filt_thread[tid];

        if (verb) INFO(("Processing common-shot gather %d of %d.", is + 1, ns));
        for (ih = 0; ih < nh; ih++) {
            int out_idx = is * nh + ih;
            c = (int)((s0 + is * ds + h0 + ih * dh - rg0) / drg + 0.5f);
            if (c < 0 || c > nrg - 1) ERROR(("Receiver table too small."));

            left = (ih - aper < 0) ? 0 : ih - aper;
            right = (ih + aper > nh - 1) ? nh - 1 : ih + aper;

            for (ic = left; ic <= right; ic++) {
                int in_idx = is * nh + ic;
                cc = (int)((s0 + is * ds + h0 + ic * dh - rg0) / drg + 0.5f);
                if (cc < 0 || cc > nrg - 1) ERROR(("Receiver table too small."));

                coef = 1.0f;
                coef *= (ic - left >= tap) ? 1.0f : (float)(ic - left) / (float)tap;
                coef *= (right - ic >= tap) ? 1.0f : (float)(right - ic) / (float)tap;

                tau = rtable[cc][c];
                dist = rdatum * rdatum + (ic - ih) * dh * (ic - ih) * dh;
                if (dist <= 0.0f) {
                    continue;
                }

                float delta = 0.0f;
                int it0 = compute_it0(tau, dt, &delta);
                if (it0 >= nt) {
                    continue;
                }

                build_kirdat_filter(tau, dt, nsam, filt);
                memset(filt_time, 0, sizeof(float) * (size_t)nfft);
                memcpy(filt_time, filt.data(), sizeof(float) * (size_t)nfilt);
#ifdef _OPENMP
#pragma omp critical(fftw_filter_exec)
#endif
                fftwf_execute(plan_filter);

                float weight = coef / (float)M_PI * dh * rdatum * tau / dist;

                fftwf_complex *out_spec = fft_freq_accum + (size_t)out_idx * nf;
                fftwf_complex *in_spec = fft_freq + (size_t)in_idx * nf;

                for (int k = 0; k < nf; k++) {
                    float interp_r = (1.0f - delta) + delta * cos_wdt[k];
                    float interp_i = delta * sin_wdt[k];

                    float phase = -2.0f * (float)M_PI * (float)k * (float)it0 / (float)nfft;
                    float cos_p = cosf(phase);
                    float sin_p = sinf(phase);

                    float hr = cos_p * interp_r - sin_p * interp_i;
                    float hi = cos_p * interp_i + sin_p * interp_r;

                    float fr = filt_freq[k][0];
                    float fi = filt_freq[k][1];

                    float kr = fr * hr - fi * hi;
                    float ki = fr * hi + fi * hr;

                    float xr = in_spec[k][0];
                    float xi = in_spec[k][1];

                    out_spec[k][0] += weight * (kr * xr - ki * xi);
                    out_spec[k][1] += weight * (kr * xi + ki * xr);
                }

                if (it0 < nt) {
                    int nmax = nfilt - 2;
                    if (nmax > nt - 1 - it0) {
                        nmax = nt - 1 - it0;
                    }
                    if (nmax >= 0) {
                        float x0 = tr_in[is][ic][0];
                        float corr_scale = -weight * delta * x0;
                        float *corr_trace = corr + (size_t)out_idx * nt + it0;
                        for (int nidx = 0; nidx <= nmax; nidx++) {
                            corr_trace[nidx] += corr_scale * filt[nidx + 1];
                        }
                    }
                }
            }
        }
    }

    fftwf_execute(plan_inverse);
    apply_ifft_output(tr_out, fft_time_out, corr, nt, nh, ns, nfft, 1.0f / (float)nfft);

    if (NULL != interm_f) {
        se_fsio_write_float(interm->data->io, tr_out[0][0], nt * nh * ns);
    }

    fill_fft_input(fft_time, tr_out, nt, nh, ns, nfft);
    fftwf_execute(plan_forward);

    memset(fft_freq_accum, 0, sizeof(fftwf_complex) * (size_t)ntraces * nf);
    memset(corr, 0, sizeof(float) * (size_t)ntraces * nt);

    s = fabsf((ns - 1) * ds);
    h = fabsf((nh - 1) * dh);

    if (fabsf(ds) >= fabsf(dh)) {
        dr = fabsf(dh);
        jump = 1;
    } else {
        dr = fabsf(ds);
        jump = (int)(dh / ds + 0.5f);
    }

    nr = (int)((s + h) / dr + 1.5f);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (ir = 0; ir < nr; ir++) {
    #ifdef _OPENMP
        int tid = omp_get_thread_num();
    #else
        int tid = 0;
    #endif
        float *filt_time = filt_time_thread[tid];
        fftwf_complex *filt_freq = filt_freq_thread[tid];
        fftwf_plan plan_filter = plan_filter_thread[tid];
        std::vector<float> &filt = filt_thread[tid];

        if (verb) INFO(("Processing common-receiver gather %d of %d.", ir + 1, nr));

        r = ir * dr + ((ds <= 0.0f) ? -1.0f : 0.0f) * s + ((dh <= 0.0f) ? -1.0f : 0.0f) * h;

        sleft = (int)((ir * dr + ((ds <= 0.0f) ? -1.0f : 0.0f) * s + ((ds <= 0.0f) ? 0.0f : -1.0f) * h) / ds + 0.5f);
        sright = (int)((ir * dr + ((ds <= 0.0f) ? -1.0f : 0.0f) * s + ((ds <= 0.0f) ? -1.0f : 0.0f) * h) / ds + 0.5f);

        if (sleft < 0) sleft = 0;
        if (sright > ns - 1) sright = ns - 1;

        left = (int)((r - sleft * ds) / dh + 0.5f);
        if (left < 0 || left > nh - 1) sleft++;

        right = (int)((r - sright * ds) / dh + 0.5f);
        if (right < 0 || right > nh - 1) sright--;

        for (is = sleft; is <= sright; is = is + jump) {
            c = (int)((s0 + is * ds - sg0) / dsg + 0.5f);
            if (c < 0 || c > nsg - 1) ERROR(("Source table too small."));

            ih = (int)((r - is * ds) / dh + 0.5f);
            if (ih < 0 || ih > nh - 1) {
                continue;
            }

            int out_idx = is * nh + ih;

            left = (is - jump * aper < sleft) ? sleft : is - jump * aper;
            right = (is + jump * aper > sright) ? sright : is + jump * aper;

            for (ic = left; ic <= right; ic = ic + jump) {
                cc = (int)((s0 + ic * ds - sg0) / dsg + 0.5f);
                if (cc < 0 || cc > nsg - 1) ERROR(("Source table too small."));

                hh = (int)((r - ic * ds) / dh + 0.5f);
                if (hh < 0 || hh > nh - 1) {
                    continue;
                }

                coef = 1.0f;
                coef *= (ic - left >= tap) ? 1.0f : (float)(ic - left) / (float)jump / (float)tap;
                coef *= (right - ic >= tap) ? 1.0f : (float)(right - ic) / (float)jump / (float)tap;

                tau = stable[cc][c];
                dist = sdatum * sdatum + (ic - is) * ds * (ic - is) * ds;
                if (dist <= 0.0f) {
                    continue;
                }

                float delta = 0.0f;
                int it0 = compute_it0(tau, dt, &delta);
                if (it0 >= nt) {
                    continue;
                }

                build_kirdat_filter(tau, dt, nsam, filt);
                memset(filt_time, 0, sizeof(float) * (size_t)nfft);
                memcpy(filt_time, filt.data(), sizeof(float) * (size_t)nfilt);
#ifdef _OPENMP
#pragma omp critical(fftw_filter_exec)
#endif
                fftwf_execute(plan_filter);

                float weight = coef / (float)M_PI * ds * sdatum * tau / dist;

                int in_idx = ic * nh + hh;
                fftwf_complex *out_spec = fft_freq_accum + (size_t)out_idx * nf;
                fftwf_complex *in_spec = fft_freq + (size_t)in_idx * nf;

                for (int k = 0; k < nf; k++) {
                    float interp_r = (1.0f - delta) + delta * cos_wdt[k];
                    float interp_i = delta * sin_wdt[k];

                    float phase = -2.0f * (float)M_PI * (float)k * (float)it0 / (float)nfft;
                    float cos_p = cosf(phase);
                    float sin_p = sinf(phase);

                    float hr = cos_p * interp_r - sin_p * interp_i;
                    float hi = cos_p * interp_i + sin_p * interp_r;

                    float fr = filt_freq[k][0];
                    float fi = filt_freq[k][1];

                    float kr = fr * hr - fi * hi;
                    float ki = fr * hi + fi * hr;

                    float xr = in_spec[k][0];
                    float xi = in_spec[k][1];

                    out_spec[k][0] += weight * (kr * xr - ki * xi);
                    out_spec[k][1] += weight * (kr * xi + ki * xr);
                }

                if (it0 < nt) {
                    int nmax = nfilt - 2;
                    if (nmax > nt - 1 - it0) {
                        nmax = nt - 1 - it0;
                    }
                    if (nmax >= 0) {
                        float x0 = tr_out[ic][hh][0];
                        float corr_scale = -weight * delta * x0;
                        float *corr_trace = corr + (size_t)out_idx * nt + it0;
                        for (int nidx = 0; nidx <= nmax; nidx++) {
                            corr_trace[nidx] += corr_scale * filt[nidx + 1];
                        }
                    }
                }
            }
        }
    }

    fftwf_execute(plan_inverse);
    apply_ifft_output(tr_in, fft_time_out, corr, nt, nh, ns, nfft, 1.0f / (float)nfft);

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

    sep_close(in);
    sep_close(out);
    if (interm_f != NULL) {
        sep_close(interm);
    }

    fftwf_destroy_plan(plan_forward);
    fftwf_destroy_plan(plan_inverse);

    for (int t = 0; t < max_threads; t++) {
        if (plan_filter_thread[t]) {
            fftwf_destroy_plan(plan_filter_thread[t]);
        }
    }

    fftwf_free(fft_time);
    fftwf_free(fft_time_out);
    fftwf_free(fft_freq);
    fftwf_free(fft_freq_accum);

    for (int t = 0; t < max_threads; t++) {
        if (filt_time_thread[t]) {
            fftwf_free(filt_time_thread[t]);
        }
        if (filt_freq_thread[t]) {
            fftwf_free(filt_freq_thread[t]);
        }
    }

    free1float(corr);
    free3float(tr_in);
    free3float(tr_out);
    free2float(stable);
    free2float(rtable);

    INFO(("Done."));
    return 0;
}
