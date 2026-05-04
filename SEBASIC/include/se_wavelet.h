#ifndef SE_WAVELET_H
#define SE_WAVELET_H

#ifdef __cplusplus
extern "C" {
#endif

/* Ricker Wavelet Parameters */
typedef struct se_ricker_params_s {
    float f0;       /* Peak frequency (Hz) */
    float amp;      /* Peak amplitude */
    float dt;       /* Sampling interval (s) */
    int nt;         /* Number of samples */
    float t0;       /* Time delay of peak (s) */
} se_ricker_params_t;

/* Gaussian Wavelet Parameters */
typedef struct se_gaussian_params_s {
    float sigma;    /* Standard deviation (width control) */
    float amp;      /* Peak amplitude */
    float dt;       /* Sampling interval (s) */
    int nt;         /* Number of samples */
    float t0;       /* Time delay of peak (s) */
} se_gaussian_params_t;

/* Ricker Wavelet Frequency Domain Parameters */
typedef struct se_ricker_freq_params_s {
    float f0;       /* Peak frequency (Hz) */
    float amp;      /* Peak amplitude in time domain */
    float df;       /* Frequency interval (Hz) */
    int nf;         /* Number of frequency samples */
    float t0;       /* Time delay of peak (s) */
} se_ricker_freq_params_t;

/* Gaussian Wavelet Frequency Domain Parameters */
typedef struct se_gaussian_freq_params_s {
    float sigma;    /* Standard deviation (width control) */
    float amp;      /* Peak amplitude in time domain */
    float df;       /* Frequency interval (Hz) */
    int nf;         /* Number of frequency samples */
    float t0;       /* Time delay of peak (s) */
} se_gaussian_freq_params_t;

/* Function Prototypes */

/**
 * @brief Generates a Ricker wavelet.
 * 
 * @param params Pointer to parameter structure.
 * @param wavelet Output array (must be allocated with size params->nt).
 */
void se_wavelet_ricker(se_ricker_params_t *params, float *wavelet);

/**
 * @brief Generates a Gaussian wavelet.
 * 
 * @param params Pointer to parameter structure.
 * @param wavelet Output array (must be allocated with size params->nt).
 */
void se_wavelet_gaussian(se_gaussian_params_t *params, float *wavelet);

/**
 * @brief Generates a Ricker wavelet in frequency domain.
 * 
 * @param params Pointer to parameter structure.
 * @param real Output real part array (must be allocated with size params->nf).
 * @param imag Output imaginary part array (must be allocated with size params->nf).
 */
void se_wavelet_ricker_freq(se_ricker_freq_params_t *params, float *real, float *imag);

/**
 * @brief Generates a Gaussian wavelet in frequency domain.
 * 
 * @param params Pointer to parameter structure.
 * @param real Output real part array (must be allocated with size params->nf).
 * @param imag Output imaginary part array (must be allocated with size params->nf).
 */
void se_wavelet_gaussian_freq(se_gaussian_freq_params_t *params, float *real, float *imag);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#include <complex>

/**
 * @brief Generates a Ricker wavelet in frequency domain (Complex output).
 * 
 * @param params Pointer to parameter structure.
 * @param wavelet Output complex array (must be allocated with size params->nf).
 */
void se_wavelet_ricker_freq(se_ricker_freq_params_t *params, std::complex<float> *wavelet);

/**
 * @brief Generates a Gaussian wavelet in frequency domain (Complex output).
 * 
 * @param params Pointer to parameter structure.
 * @param wavelet Output complex array (must be allocated with size params->nf).
 */
void se_wavelet_gaussian_freq(se_gaussian_freq_params_t *params, std::complex<float> *wavelet);

#endif

#endif /* SE_WAVELET_H */