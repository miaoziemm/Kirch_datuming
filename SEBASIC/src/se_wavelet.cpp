#include "se_wavelet.h"
#include <math.h>
#include <stdlib.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

void se_wavelet_ricker(se_ricker_params_t *params, float *wavelet)
{
    if (!params || !wavelet) return;

    float f0 = params->f0;
    float amp = params->amp;
    float dt = params->dt;
    int nt = params->nt;
    float t0 = params->t0;

    float pi2 = PI * PI;
    float f02 = f0 * f0;

    for (int i = 0; i < nt; ++i) {
        float t = i * dt - t0;
        float t2 = t * t;
        float arg = pi2 * f02 * t2;
        
        wavelet[i] = amp * (1.0f - 2.0f * arg) * expf(-arg);
    }
}

void se_wavelet_gaussian(se_gaussian_params_t *params, float *wavelet)
{
    if (!params || !wavelet) return;

    float sigma = params->sigma;
    float amp = params->amp;
    float dt = params->dt;
    int nt = params->nt;
    float t0 = params->t0;

    /* Avoid division by zero */
    if (sigma == 0.0f) return;

    float sigma2 = sigma * sigma;

    for (int i = 0; i < nt; ++i) {
        float t = i * dt - t0;
        float t2 = t * t;
        
        wavelet[i] = amp * expf(-t2 / (2.0f * sigma2));
    }
}

void se_wavelet_ricker_freq(se_ricker_freq_params_t *params, float *real, float *imag)
{
    if (!params || !real || !imag) return;

    float f0 = params->f0;
    float amp = params->amp;
    float df = params->df;
    int nf = params->nf;
    float t0 = params->t0;

    /* Avoid division by zero */
    if (f0 == 0.0f) return;

    float f02 = f0 * f0;
    float f03 = f02 * f0;
    float sqrt_pi = sqrtf((float)PI);
    float factor = 2.0f / (sqrt_pi * f03);

    for (int i = 0; i < nf; ++i) {
        float f = i * df;
        float f2 = f * f;
        
        /* Amplitude spectrum: A(f) = (2 * f^2 / (sqrt(pi) * f0^3)) * exp(-f^2 / f0^2) */
        float A = amp * factor * f2 * expf(-f2 / f02);
        
        /* Phase shift: exp(-i * 2 * pi * f * t0) */
        float phase = -2.0f * (float)PI * f * t0;
        
        real[i] = A * cosf(phase);
        imag[i] = A * sinf(phase);
    }
}

void se_wavelet_ricker_freq(se_ricker_freq_params_t *params, std::complex<float> *wavelet)
{
    if (!params || !wavelet) return;

    float f0 = params->f0;
    float amp = params->amp;
    float df = params->df;
    int nf = params->nf;
    float t0 = params->t0;

    /* Avoid division by zero */
    if (f0 == 0.0f) return;

    float f02 = f0 * f0;
    float f03 = f02 * f0;
    float sqrt_pi = sqrtf((float)PI);
    float factor = 2.0f / (sqrt_pi * f03);

    for (int i = 0; i < nf; ++i) {
        float f = i * df;
        float f2 = f * f;
        
        /* Amplitude spectrum: A(f) = (2 * f^2 / (sqrt(pi) * f0^3)) * exp(-f^2 / f0^2) */
        float A = amp * factor * f2 * expf(-f2 / f02);
        
        /* Phase shift: exp(-i * 2 * pi * f * t0) */
        float phase = -2.0f * (float)PI * f * t0;
        
        wavelet[i] = std::complex<float>(A * cosf(phase), A * sinf(phase));
    }
}

void se_wavelet_gaussian_freq(se_gaussian_freq_params_t *params, float *real, float *imag)
{
    if (!params || !real || !imag) return;

    float sigma = params->sigma;
    float amp = params->amp;
    float df = params->df;
    int nf = params->nf;
    float t0 = params->t0;

    float sigma2 = sigma * sigma;
    float sqrt_2pi = sqrtf(2.0f * (float)PI);
    float pi2 = (float)PI * (float)PI;
    float factor = sigma * sqrt_2pi;

    for (int i = 0; i < nf; ++i) {
        float f = i * df;
        float f2 = f * f;
        
        /* Amplitude spectrum: A(f) = sigma * sqrt(2*pi) * exp(-2 * pi^2 * sigma^2 * f^2) */
        float A = amp * factor * expf(-2.0f * pi2 * sigma2 * f2);
        
        /* Phase shift: exp(-i * 2 * pi * f * t0) */
        float phase = -2.0f * (float)PI * f * t0;
        
        real[i] = A * cosf(phase);
        imag[i] = A * sinf(phase);
    }
}

void se_wavelet_gaussian_freq(se_gaussian_freq_params_t *params, std::complex<float> *wavelet)
{
    if (!params || !wavelet) return;

    float sigma = params->sigma;
    float amp = params->amp;
    float df = params->df;
    int nf = params->nf;
    float t0 = params->t0;

    float sigma2 = sigma * sigma;
    float sqrt_2pi = sqrtf(2.0f * (float)PI);
    float pi2 = (float)PI * (float)PI;
    float factor = sigma * sqrt_2pi;

    for (int i = 0; i < nf; ++i) {
        float f = i * df;
        float f2 = f * f;
        
        /* Amplitude spectrum: A(f) = sigma * sqrt(2*pi) * exp(-2 * pi^2 * sigma^2 * f^2) */
        float A = amp * factor * expf(-2.0f * pi2 * sigma2 * f2);
        
        /* Phase shift: exp(-i * 2 * pi * f * t0) */
        float phase = -2.0f * (float)PI * f * t0;
        
        wavelet[i] = std::complex<float>(A * cosf(phase), A * sinf(phase));
    }
}