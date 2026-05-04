#ifndef SE_MATH_PHYSICS_H
#define SE_MATH_PHYSICS_H

#ifdef __cplusplus
#include <complex>
extern "C" {
#endif

/**
 * @brief Bessel function of the first kind of order 0.
 * Compatible wrapper for j0.
 */
double se_bessel_j0(double x);

/**
 * @brief Bessel function of the second kind of order 0.
 * Compatible wrapper for y0.
 */
double se_bessel_y0(double x);

/**
 * @brief Bessel function of the first kind of order 1.
 * Compatible wrapper for j1.
 */
double se_bessel_j1(double x);

/**
 * @brief Bessel function of the second kind of order 1.
 * Compatible wrapper for y1.
 */
double se_bessel_y1(double x);

/**
 * @brief Bessel function of the first kind of order n.
 * Compatible wrapper for jn.
 */
double se_bessel_jn(int n, double x);

/**
 * @brief Bessel function of the second kind of order n.
 * Compatible wrapper for yn.
 */
double se_bessel_yn(int n, double x);

/**
 * @brief Hankel function of the first kind of order 0.
 * H0(1)(x) = J0(x) + i Y0(x)
 * Returns real and imaginary parts via pointers.
 */
void se_bessel_hankel1_0(double x, double *real, double *imag);

/**
 * @brief Hankel function of the second kind of order 0.
 * H0(2)(x) = J0(x) - i Y0(x)
 * Returns real and imaginary parts via pointers.
 */
void se_bessel_hankel2_0(double x, double *real, double *imag);

/**
 * @brief Hankel function of the first kind of order 1.
 * H1(1)(x) = J1(x) + i Y1(x)
 * Returns real and imaginary parts via pointers.
 */
void se_bessel_hankel1_1(double x, double *real, double *imag);

/**
 * @brief Hankel function of the second kind of order 1.
 * H1(2)(x) = J1(x) - i Y1(x)
 * Returns real and imaginary parts via pointers.
 */
void se_bessel_hankel2_1(double x, double *real, double *imag);

/**
 * @brief Hankel function of the first kind of order n.
 * Hn(1)(x) = Jn(x) + i Yn(x)
 * Returns real and imaginary parts via pointers.
 */
void se_bessel_hankel1_n(int n, double x, double *real, double *imag);

/**
 * @brief Hankel function of the second kind of order n.
 * Hn(2)(x) = Jn(x) - i Yn(x)
 * Returns real and imaginary parts via pointers.
 */
void se_bessel_hankel2_n(int n, double x, double *real, double *imag);

/**
 * @brief Legendre polynomials of the first kind Pn(x).
 * Defined for |x| <= 1.
 */
double se_legendre_P(int n, double x);

/**
 * @brief Spherical Bessel function of the first kind jn(x).
 */
double se_sph_bessel_j(int n, double x);

/**
 * @brief Spherical Bessel function of the second kind yn(x) (Spherical Neumann).
 */
double se_sph_bessel_y(int n, double x);

/**
 * @brief Hermite polynomials Hn(x).
 * Physically relevant for Gaussian beams.
 */
double se_hermite(int n, double x);

/**
 * @brief Error function erf(x).
 */
double se_math_erf(double x);

/**
 * @brief Complementary error function erfc(x).
 */
double se_math_erfc(double x);

/**
 * @brief Compute the Hilbert transform of a real signal.
 * Uses FFT (float precision internally).
 * @param n Number of samples.
 * @param input Input signal (real).
 * @param output Output signal (Hilbert transform, real).
 */
void se_hilbert(int n, const double* input, double* output);

/**
 * @brief Compute the Analytic signal of a real signal.
 * Uses FFT (float precision internally).
 * @param n Number of samples.
 * @param input Input signal (real).
 * @param output Output signal (Analytic signal, complex).
 */
void se_analytic_signal(int n, const double* input, std::complex<double>* output);

#ifdef __cplusplus
}

/* C++ Overloads returning std::complex */

/**
 * @brief Hankel function of the first kind of order 0.
 * Returns std::complex<double>.
 */
std::complex<double> se_bessel_hankel1_0(double x);

/**
 * @brief Hankel function of the second kind of order 0.
 * Returns std::complex<double>.
 */
std::complex<double> se_bessel_hankel2_0(double x);

/**
 * @brief Hankel function of the first kind of order 1.
 * Returns std::complex<double>.
 */
std::complex<double> se_bessel_hankel1_1(double x);

/**
 * @brief Hankel function of the second kind of order 1.
 * Returns std::complex<double>.
 */
std::complex<double> se_bessel_hankel2_1(double x);

/**
 * @brief Hankel function of the first kind of order n.
 * Returns std::complex<double>.
 */
std::complex<double> se_bessel_hankel1_n(int n, double x);

/**
 * @brief Hankel function of the second kind of order n.
 * Returns std::complex<double>.
 */
std::complex<double> se_bessel_hankel2_n(int n, double x);

/* Float Overloads */

inline float se_bessel_j0(float x) { return (float)se_bessel_j0((double)x); }
inline float se_bessel_y0(float x) { return (float)se_bessel_y0((double)x); }
inline float se_bessel_j1(float x) { return (float)se_bessel_j1((double)x); }
inline float se_bessel_y1(float x) { return (float)se_bessel_y1((double)x); }
inline float se_bessel_jn(int n, float x) { return (float)se_bessel_jn(n, (double)x); }
inline float se_bessel_yn(int n, float x) { return (float)se_bessel_yn(n, (double)x); }

inline float se_legendre_P(int n, float x) { return (float)se_legendre_P(n, (double)x); }
inline float se_sph_bessel_j(int n, float x) { return (float)se_sph_bessel_j(n, (double)x); }
inline float se_sph_bessel_y(int n, float x) { return (float)se_sph_bessel_y(n, (double)x); }
inline float se_hermite(int n, float x) { return (float)se_hermite(n, (double)x); }
inline float se_math_erf(float x) { return (float)se_math_erf((double)x); }
inline float se_math_erfc(float x) { return (float)se_math_erfc((double)x); }

inline std::complex<float> se_bessel_hankel1_0(float x) {
    std::complex<double> res = se_bessel_hankel1_0((double)x);
    return std::complex<float>((float)res.real(), (float)res.imag());
}
inline std::complex<float> se_bessel_hankel2_0(float x) {
    std::complex<double> res = se_bessel_hankel2_0((double)x);
    return std::complex<float>((float)res.real(), (float)res.imag());
}
inline std::complex<float> se_bessel_hankel1_1(float x) {
    std::complex<double> res = se_bessel_hankel1_1((double)x);
    return std::complex<float>((float)res.real(), (float)res.imag());
}
inline std::complex<float> se_bessel_hankel2_1(float x) {
    std::complex<double> res = se_bessel_hankel2_1((double)x);
    return std::complex<float>((float)res.real(), (float)res.imag());
}
inline std::complex<float> se_bessel_hankel1_n(int n, float x) {
    std::complex<double> res = se_bessel_hankel1_n(n, (double)x);
    return std::complex<float>((float)res.real(), (float)res.imag());
}
inline std::complex<float> se_bessel_hankel2_n(int n, float x) {
    std::complex<double> res = se_bessel_hankel2_n(n, (double)x);
    return std::complex<float>((float)res.real(), (float)res.imag());
}

/**
 * @brief Compute the Hilbert transform of a real signal (float version).
 */
void se_hilbert(int n, const float* input, float* output);

/**
 * @brief Compute the Analytic signal of a real signal (float version).
 */
void se_analytic_signal(int n, const float* input, std::complex<float>* output);

/**
 * @brief Apply Butterworth Bandpass filter.
 * @param n Number of samples.
 * @param input Input signal.
 * @param output Output signal.
 * @param dt Sampling interval (seconds).
 * @param low_cut Low cut frequency (Hz).
 * @param high_cut High cut frequency (Hz).
 * @param order Filter order (poles).
 * @param zero_phase If 1, apply zero-phase filtering (forward and backward).
 */
void se_filter_bandpass(int n, const double* input, double* output, double dt, double low_cut, double high_cut, int order, int zero_phase);

/**
 * @brief Apply Butterworth Lowpass filter.
 * @param n Number of samples.
 * @param input Input signal.
 * @param output Output signal.
 * @param dt Sampling interval (seconds).
 * @param cutoff Cutoff frequency (Hz).
 * @param order Filter order.
 * @param zero_phase If 1, apply zero-phase filtering.
 */
void se_filter_lowpass(int n, const double* input, double* output, double dt, double cutoff, int order, int zero_phase);

/**
 * @brief Apply Butterworth Highpass filter.
 * @param n Number of samples.
 * @param input Input signal.
 * @param output Output signal.
 * @param dt Sampling interval (seconds).
 * @param cutoff Cutoff frequency (Hz).
 * @param order Filter order.
 * @param zero_phase If 1, apply zero-phase filtering.
 */
void se_filter_highpass(int n, const double* input, double* output, double dt, double cutoff, int order, int zero_phase);

/* Float Overloads for Filters */
void se_filter_bandpass(int n, const float* input, float* output, double dt, double low_cut, double high_cut, int order, int zero_phase);
void se_filter_lowpass(int n, const float* input, float* output, double dt, double cutoff, int order, int zero_phase);
void se_filter_highpass(int n, const float* input, float* output, double dt, double cutoff, int order, int zero_phase);

/* Wavelet Transform */

typedef enum {
    SE_WAVELET_HAAR = 0,
    SE_WAVELET_DB4  = 1
} se_wavelet_type_t;

/**
 * @brief Discrete Wavelet Transform (1D).
 * Decomposes signal into [Approximation_LevelN, Detail_LevelN, Detail_LevelN-1, ..., Detail_Level1].
 * Output array must be same size as input (or slightly larger for padding, but here we assume power of 2 or handle padding internally).
 * Ideally n should be power of 2. If not, zero-padding is applied internally but output size is n.
 * 
 * @param n Input signal length.
 * @param input Input signal array.
 * @param output Output coefficients array (size n).
 * @param level Decomposition level.
 * @param type Wavelet type (SE_WAVELET_HAAR or SE_WAVELET_DB4).
 */
void se_dwt(int n, const double* input, double* output, int level, se_wavelet_type_t type);

/**
 * @brief Inverse Discrete Wavelet Transform (1D).
 * Reconstructs signal from coefficients.
 * 
 * @param n Signal length (must match the DWT input length).
 * @param input Input coefficients array (from se_dwt).
 * @param output Reconstructed signal array.
 * @param level Decomposition level (must match DWT).
 * @param type Wavelet type.
 */
void se_idwt(int n, const double* input, double* output, int level, se_wavelet_type_t type);

/* Float Overloads for Wavelet */
void se_dwt(int n, const float* input, float* output, int level, se_wavelet_type_t type);
void se_idwt(int n, const float* input, float* output, int level, se_wavelet_type_t type);

#endif

#endif // SE_MATH_PHYSICS_H
