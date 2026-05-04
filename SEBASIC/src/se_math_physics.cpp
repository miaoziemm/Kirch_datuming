#include "se_math_physics.h"
#include "se_su_fft.h"
#include "../include/se_bandpass.h"
#include <cmath>
#include <math.h>
#include <vector>
#include <complex>

/*
 * Implementation of physics-related math functions.
 * Provides compatibility wrappers for Bessel functions and others.
 */

double se_bessel_j0(double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    // C++17 standard
    return std::cyl_bessel_j(0.0, x);
#elif defined(_MSC_VER)
    // MSVC specific
    return _j0(x);
#else
    // POSIX / BSD standard (usually available in math.h)
    return j0(x);
#endif
}

double se_bessel_y0(double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    // C++17 standard
    return std::cyl_neumann(0.0, x);
#elif defined(_MSC_VER)
    // MSVC specific
    return _y0(x);
#else
    // POSIX / BSD standard (usually available in math.h)
    return y0(x);
#endif
}

double se_bessel_j1(double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    // C++17 standard
    return std::cyl_bessel_j(1.0, x);
#elif defined(_MSC_VER)
    // MSVC specific
    return _j1(x);
#else
    // POSIX / BSD standard (usually available in math.h)
    return j1(x);
#endif
}

double se_bessel_y1(double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    // C++17 standard
    return std::cyl_neumann(1.0, x);
#elif defined(_MSC_VER)
    // MSVC specific
    return _y1(x);
#else
    // POSIX / BSD standard (usually available in math.h)
    return y1(x);
#endif
}

double se_bessel_jn(int n, double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    // C++17 standard
    return std::cyl_bessel_j((double)n, x);
#elif defined(_MSC_VER)
    // MSVC specific
    return _jn(n, x);
#else
    // POSIX / BSD standard (usually available in math.h)
    return jn(n, x);
#endif
}

double se_bessel_yn(int n, double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    // C++17 standard
    return std::cyl_neumann((double)n, x);
#elif defined(_MSC_VER)
    // MSVC specific
    return _yn(n, x);
#else
    // POSIX / BSD standard (usually available in math.h)
    return yn(n, x);
#endif
}

void se_bessel_hankel1_0(double x, double *real, double *imag) {
    if (real) *real = se_bessel_j0(x);
    if (imag) *imag = se_bessel_y0(x);
}

void se_bessel_hankel2_0(double x, double *real, double *imag) {
    if (real) *real = se_bessel_j0(x);
    if (imag) *imag = -se_bessel_y0(x);
}

void se_bessel_hankel1_1(double x, double *real, double *imag) {
    if (real) *real = se_bessel_j1(x);
    if (imag) *imag = se_bessel_y1(x);
}

void se_bessel_hankel2_1(double x, double *real, double *imag) {
    if (real) *real = se_bessel_j1(x);
    if (imag) *imag = -se_bessel_y1(x);
}

void se_bessel_hankel1_n(int n, double x, double *real, double *imag) {
    if (real) *real = se_bessel_jn(n, x);
    if (imag) *imag = se_bessel_yn(n, x);
}

void se_bessel_hankel2_n(int n, double x, double *real, double *imag) {
    if (real) *real = se_bessel_jn(n, x);
    if (imag) *imag = -se_bessel_yn(n, x);
}

double se_legendre_P(int n, double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    return std::legendre(n, x);
#else
    // Recurrence relation: (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)
    if (n < 0) return 0.0;
    if (n == 0) return 1.0;
    if (n == 1) return x;
    
    double p0 = 1.0;
    double p1 = x;
    double pn = 0.0;
    
    for (int i = 1; i < n; ++i) {
        pn = ((2.0 * i + 1.0) * x * p1 - i * p0) / (i + 1.0);
        p0 = p1;
        p1 = pn;
    }
    return p1;
#endif
}

double se_sph_bessel_j(int n, double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    return std::sph_bessel(n, x);
#else
    // Explicit formulas for low orders or recurrence
    // j_n(x) = sqrt(pi/(2x)) * J_{n+0.5}(x)
    // Recurrence: j_{n+1}(x) = (2n+1)/x * j_n(x) - j_{n-1}(x)
    // Note: Forward recurrence is unstable for j_n. Backward is better.
    // But for simplicity and low n, we use explicit or forward (with caution).
    // Better: Use explicit for n=0,1 and forward.
    
    if (n < 0) return 0.0;
    if (std::abs(x) < 1e-10) return (n == 0) ? 1.0 : 0.0;

    double j0 = std::sin(x) / x;
    if (n == 0) return j0;
    
    double j1 = j0 / x - std::cos(x) / x;
    if (n == 1) return j1;

    double j_prev2 = j0;
    double j_prev1 = j1;
    double j_curr = 0.0;

    for (int i = 1; i < n; ++i) {
        j_curr = (2.0 * i + 1.0) / x * j_prev1 - j_prev2;
        j_prev2 = j_prev1;
        j_prev1 = j_curr;
    }
    return j_prev1;
#endif
}

double se_sph_bessel_y(int n, double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    return std::sph_neumann(n, x);
#else
    // y_0(x) = -cos(x)/x
    // y_1(x) = -cos(x)/x^2 - sin(x)/x
    // Recurrence: y_{n+1}(x) = (2n+1)/x * y_n(x) - y_{n-1}(x)
    
    if (n < 0) return 0.0;
    if (std::abs(x) < 1e-10) return -std::numeric_limits<double>::infinity();

    double y0 = -std::cos(x) / x;
    if (n == 0) return y0;

    double y1 = y0 / x - std::sin(x) / x;
    if (n == 1) return y1;

    double y_prev2 = y0;
    double y_prev1 = y1;
    double y_curr = 0.0;

    for (int i = 1; i < n; ++i) {
        y_curr = (2.0 * i + 1.0) / x * y_prev1 - y_prev2;
        y_prev2 = y_prev1;
        y_prev1 = y_curr;
    }
    return y_prev1;
#endif
}

double se_hermite(int n, double x) {
#if defined(__cplusplus) && __cplusplus >= 201703L && defined(__cpp_lib_math_special_functions)
    return std::hermite(n, x);
#else
    // H_0(x) = 1
    // H_1(x) = 2x
    // H_{n+1}(x) = 2xH_n(x) - 2nH_{n-1}(x)
    if (n < 0) return 0.0;
    if (n == 0) return 1.0;
    if (n == 1) return 2.0 * x;

    double h0 = 1.0;
    double h1 = 2.0 * x;
    double hn = 0.0;

    for (int i = 1; i < n; ++i) {
        hn = 2.0 * x * h1 - 2.0 * i * h0;
        h0 = h1;
        h1 = hn;
    }
    return h1;
#endif
}

double se_math_erf(double x) {
    return erf(x);
}

double se_math_erfc(double x) {
    return erfc(x);
}

/* C++ Overloads implementation */

std::complex<double> se_bessel_hankel1_0(double x) {
    return std::complex<double>(se_bessel_j0(x), se_bessel_y0(x));
}

std::complex<double> se_bessel_hankel2_0(double x) {
    return std::complex<double>(se_bessel_j0(x), -se_bessel_y0(x));
}

std::complex<double> se_bessel_hankel1_1(double x) {
    return std::complex<double>(se_bessel_j1(x), se_bessel_y1(x));
}

std::complex<double> se_bessel_hankel2_1(double x) {
    return std::complex<double>(se_bessel_j1(x), -se_bessel_y1(x));
}

std::complex<double> se_bessel_hankel1_n(int n, double x) {
    return std::complex<double>(se_bessel_jn(n, x), se_bessel_yn(n, x));
}

std::complex<double> se_bessel_hankel2_n(int n, double x) {
    return std::complex<double>(se_bessel_jn(n, x), -se_bessel_yn(n, x));
}

void se_analytic_signal(int n, const double* input, std::complex<double>* output) {
    if (n <= 0 || !input || !output) return;

    // 1. Determine FFT size (Prime Factor Algorithm)
    int nfft = su_fft_npfa(n);
    
    // 2. Prepare complex buffer
    std::vector<std::complex<float>> z(nfft);
    for (int i = 0; i < n; ++i) {
        z[i] = std::complex<float>((float)input[i], 0.0f);
    }
    // Zero padding
    for (int i = n; i < nfft; ++i) {
        z[i] = std::complex<float>(0.0f, 0.0f);
    }

    // 3. Forward FFT
    // isign = -1 for forward transform (standard convention)
    su_fft_pfacc(-1, nfft, z.data());

    // 4. Apply Hilbert mask in frequency domain
    // H(f) = 1 for f=0
    // H(f) = 2 for f > 0
    // H(f) = 0 for f < 0
    // Nyquist (if even) is usually kept as is or doubled? 
    // For analytic signal, we want to suppress negative frequencies and double positive ones.
    // DC component is kept as is.
    // Nyquist component (at nfft/2) is shared, usually kept as is (or 0 for strict analytic?).
    // Standard approach:
    // z[0] *= 1
    // z[1 ... nfft/2 - 1] *= 2
    // z[nfft/2] *= 1 (if nfft even)
    // z[nfft/2 + 1 ... nfft - 1] *= 0

    int half = nfft / 2;
    for (int i = 1; i < half; ++i) {
        z[i] *= 2.0f;
    }
    if (nfft % 2 == 0) {
        // Nyquist frequency
        // z[half] *= 1.0f; // Keep as is
    } else {
        // If odd, half is (nfft-1)/2, which is the last positive freq
        z[half] *= 2.0f;
    }
    
    // Zero out negative frequencies
    for (int i = half + 1; i < nfft; ++i) {
        z[i] = 0.0f;
    }

    // 5. Inverse FFT
    su_fft_pfacc(1, nfft, z.data());

    // 6. Scale and output
    float scale = 1.0f / nfft;
    for (int i = 0; i < n; ++i) {
        output[i] = std::complex<double>((double)(z[i].real() * scale), (double)(z[i].imag() * scale));
    }
}

void se_hilbert(int n, const double* input, double* output) {
    if (n <= 0 || !input || !output) return;
    
    // Reuse se_analytic_signal logic or reimplement for efficiency?
    // Reimplementing to avoid allocating complex<double> buffer if not needed.
    
    int nfft = su_fft_npfa(n);
    std::vector<std::complex<float>> z(nfft);
    
    for (int i = 0; i < n; ++i) {
        z[i] = std::complex<float>((float)input[i], 0.0f);
    }
    for (int i = n; i < nfft; ++i) {
        z[i] = std::complex<float>(0.0f, 0.0f);
    }

    su_fft_pfacc(-1, nfft, z.data());

    int half = nfft / 2;
    for (int i = 1; i < half; ++i) {
        z[i] *= 2.0f;
    }
    if (nfft % 2 != 0) {
        z[half] *= 2.0f;
    }
    for (int i = half + 1; i < nfft; ++i) {
        z[i] = 0.0f;
    }

    su_fft_pfacc(1, nfft, z.data());

    float scale = 1.0f / nfft;
    for (int i = 0; i < n; ++i) {
        // Hilbert transform is the imaginary part of the analytic signal
        output[i] = (double)(z[i].imag() * scale);
    }
}

void se_analytic_signal(int n, const float* input, std::complex<float>* output) {
    if (n <= 0 || !input || !output) return;

    int nfft = su_fft_npfa(n);
    std::vector<std::complex<float>> z(nfft);
    
    for (int i = 0; i < n; ++i) {
        z[i] = std::complex<float>(input[i], 0.0f);
    }
    for (int i = n; i < nfft; ++i) {
        z[i] = std::complex<float>(0.0f, 0.0f);
    }

    su_fft_pfacc(-1, nfft, z.data());

    int half = nfft / 2;
    for (int i = 1; i < half; ++i) {
        z[i] *= 2.0f;
    }
    if (nfft % 2 != 0) {
        z[half] *= 2.0f;
    }
    for (int i = half + 1; i < nfft; ++i) {
        z[i] = 0.0f;
    }

    su_fft_pfacc(1, nfft, z.data());

    float scale = 1.0f / nfft;
    for (int i = 0; i < n; ++i) {
        output[i] = std::complex<float>(z[i].real() * scale, z[i].imag() * scale);
    }
}

void se_hilbert(int n, const float* input, float* output) {
    if (n <= 0 || !input || !output) return;
    
    int nfft = su_fft_npfa(n);
    std::vector<std::complex<float>> z(nfft);
    
    for (int i = 0; i < n; ++i) {
        z[i] = std::complex<float>(input[i], 0.0f);
    }
    for (int i = n; i < nfft; ++i) {
        z[i] = std::complex<float>(0.0f, 0.0f);
    }

    su_fft_pfacc(-1, nfft, z.data());

    int half = nfft / 2;
    for (int i = 1; i < half; ++i) {
        z[i] *= 2.0f;
    }
    if (nfft % 2 != 0) {
        z[half] *= 2.0f;
    }
    for (int i = half + 1; i < nfft; ++i) {
        z[i] = 0.0f;
    }

    su_fft_pfacc(1, nfft, z.data());

    float scale = 1.0f / nfft;
    for (int i = 0; i < n; ++i) {
        output[i] = z[i].imag() * scale;
    }
}

/* Filter Implementations */

static void se_filter_generic(int n, const float* input, float* output, double dt, double flo, double fhi, int nplo, int nphi, int zero_phase) {
    if (n <= 0 || !input || !output) return;

    // Initialize bandpass filter
    se_bandpass_t* bp = se_bandpass_init(dt, flo, fhi, nplo, nphi, zero_phase);
    if (!bp) return;

    // Copy input to output (se_apply_bandpass works in-place)
    if (input != output) {
        for (int i = 0; i < n; ++i) output[i] = input[i];
    }

    // Apply filter
    // se_apply_bandpass takes reala* (float*)
    se_apply_bandpass(bp, output, n, 1);

    // Cleanup
    se_bandpass_destroy(bp);
}

void se_filter_bandpass(int n, const float* input, float* output, double dt, double low_cut, double high_cut, int order, int zero_phase) {
    se_filter_generic(n, input, output, dt, low_cut, high_cut, order, order, zero_phase);
}

void se_filter_lowpass(int n, const float* input, float* output, double dt, double cutoff, int order, int zero_phase) {
    // flo=0 means no low cut (lowpass behavior determined by fhi)
    se_filter_generic(n, input, output, dt, 0.0, cutoff, 0, order, zero_phase);
}

void se_filter_highpass(int n, const float* input, float* output, double dt, double cutoff, int order, int zero_phase) {
    // fhi=0 means no high cut (highpass behavior determined by flo)
    // se_bandpass_init treats fhi <= 0 as Nyquist
    se_filter_generic(n, input, output, dt, cutoff, 0.0, order, 0, zero_phase);
}

/* Double versions */

void se_filter_bandpass(int n, const double* input, double* output, double dt, double low_cut, double high_cut, int order, int zero_phase) {
    if (n <= 0 || !input || !output) return;
    std::vector<float> tmp(n);
    for(int i=0; i<n; ++i) tmp[i] = (float)input[i];
    
    se_filter_bandpass(n, tmp.data(), tmp.data(), dt, low_cut, high_cut, order, zero_phase);
    
    for(int i=0; i<n; ++i) output[i] = (double)tmp[i];
}

void se_filter_lowpass(int n, const double* input, double* output, double dt, double cutoff, int order, int zero_phase) {
    if (n <= 0 || !input || !output) return;
    std::vector<float> tmp(n);
    for(int i=0; i<n; ++i) tmp[i] = (float)input[i];
    
    se_filter_lowpass(n, tmp.data(), tmp.data(), dt, cutoff, order, zero_phase);
    
    for(int i=0; i<n; ++i) output[i] = (double)tmp[i];
}

void se_filter_highpass(int n, const double* input, double* output, double dt, double cutoff, int order, int zero_phase) {
    if (n <= 0 || !input || !output) return;
    std::vector<float> tmp(n);
    for(int i=0; i<n; ++i) tmp[i] = (float)input[i];
    
    se_filter_highpass(n, tmp.data(), tmp.data(), dt, cutoff, order, zero_phase);
    
    for(int i=0; i<n; ++i) output[i] = (double)tmp[i];
}

/* Wavelet Transform Implementation */

// Helper: Single level DWT step
static void se_dwt_step(const double* in, double* out, int n, se_wavelet_type_t type) {
    // Filter coefficients
    // Haar
    static const double h0_haar[] = { 0.7071067811865475, 0.7071067811865475 }; // Low pass
    static const double h1_haar[] = { -0.7071067811865475, 0.7071067811865475 }; // High pass (reversed for convolution?) 
    // Standard DWT definition: 
    // cA[k] = sum(x[i] * h0[2k-i]) ... depends on convention.
    // Usually: cA = x * h0 (downsampled), cD = x * h1 (downsampled)
    // Haar: cA[k] = (x[2k] + x[2k+1])/sqrt(2)
    //       cD[k] = (x[2k] - x[2k+1])/sqrt(2) (or similar)

    // DB4 (4 coefficients)
    static const double h0_db4[] = { 0.4829629131445341, 0.8365163037378079, 0.2241438680420134, -0.1294095225512604 };
    static const double h1_db4[] = { -0.1294095225512604, -0.2241438680420134, 0.8365163037378079, -0.4829629131445341 };

    const double* h0 = (type == SE_WAVELET_DB4) ? h0_db4 : h0_haar;
    const double* h1 = (type == SE_WAVELET_DB4) ? h1_db4 : h1_haar;
    int len = (type == SE_WAVELET_DB4) ? 4 : 2;

    int half = n / 2;
    for (int i = 0; i < half; ++i) {
        double s0 = 0.0;
        double s1 = 0.0;
        for (int j = 0; j < len; ++j) {
            int idx = 2 * i + j;
            // Periodic boundary handling
            if (idx >= n) idx -= n; 
            
            // Standard Mallat algorithm
            // cA = conv(x, h0) downsample
            // cD = conv(x, h1) downsample
            // Note: Implementation details vary (filter order, sign). 
            // Using standard decomposition filters.
            s0 += in[idx] * h0[j];
            s1 += in[idx] * h1[j];
        }
        out[i] = s0;        // Approximation (Low freq)
        out[half + i] = s1; // Detail (High freq)
    }
}

// Helper: Single level IDWT step
static void se_idwt_step(const double* in, double* out, int n, se_wavelet_type_t type) {
    // Reconstruction filters
    // For orthogonal wavelets (Haar, DB4), reconstruction filters g0, g1 are related to h0, h1.
    // g0[n] = h0[n]
    // g1[n] = -h1[n] ... wait, usually:
    // g0[n] = h0[L-1-n]
    // g1[n] = h1[L-1-n] * (-1)^n ... 
    // Actually for orthogonal:
    // Inv Low pass = h0 (time reversed?)
    // Let's use the property: x = Upsample(cA)*g0 + Upsample(cD)*g1
    
    // Haar Reconstruction
    // x[2k] = (cA[k] + cD[k])/sqrt(2)
    // x[2k+1] = (cA[k] - cD[k])/sqrt(2)
    
    // DB4 Reconstruction
    // g0 = h0 reversed
    // g1 = h1 reversed
    // static const double g0_db4[] = { -0.1294095225512604, 0.2241438680420134, 0.8365163037378079, 0.4829629131445341 }; // h0 reversed? No.
    // Correct reconstruction filters for DB4 (from standard literature)
    // static const double rec_lo_db4[] = { 0.4829629131445341, 0.8365163037378079, 0.2241438680420134, -0.1294095225512604 }; // Same as dec_lo for orthogonal? No.
    // Actually, for orthogonal wavelets:
    // Synthesis Low Pass (g0) = Decomposition Low Pass (h0) reversed order?
    // Let's use the explicit coefficients for reconstruction to be safe.
    // DB4 Synthesis Low:  0.48296..., 0.83651..., 0.22414..., -0.12940... (Same as analysis h0? No, usually reversed)
    // Wait, for Daubechies, Analysis Low (Lo_D) and Synthesis Low (Lo_R) are related.
    // Lo_R[k] = Lo_D[N-1-k]
    // Hi_R[k] = Hi_D[N-1-k]
    
    // static const double lo_r_db4[] = { -0.1294095225512604, 0.2241438680420134, 0.8365163037378079, 0.4829629131445341 }; // h0 reversed
    // static const double hi_r_db4[] = { -0.4829629131445341, 0.8365163037378079, -0.2241438680420134, -0.1294095225512604 }; // h1 reversed

    // static const double lo_r_haar[] = { 0.7071067811865475, 0.7071067811865475 };
    // static const double hi_r_haar[] = { 0.7071067811865475, -0.7071067811865475 };

    // const double* g0 = (type == SE_WAVELET_DB4) ? lo_r_db4 : lo_r_haar;
    // const double* g1 = (type == SE_WAVELET_DB4) ? hi_r_db4 : hi_r_haar;
    int len = (type == SE_WAVELET_DB4) ? 4 : 2;

    int half = n / 2;
    // Clear output
    for(int i=0; i<n; ++i) out[i] = 0.0;

    for (int i = 0; i < half; ++i) {
        // Upsampling and convolution
        // Contribution of cA[i] (at out[2*i ...]) and cD[i] (at out[2*i ...])
        // This is tricky with periodic boundary.
        // Standard formula: x[j] = sum_k ( cA[k]*g0[j-2k] + cD[k]*g1[j-2k] )
        
        double ca = in[i];
        double cd = in[half + i];

        for (int j = 0; j < len; ++j) {
            // int out_idx = 2 * i + j - (len/2) + 1; // Alignment adjustment?
            // For Haar (len=2): 2*i + j. 
            // For DB4 (len=4): 2*i + j - 2?
            // Let's stick to simple periodic reconstruction matching the forward step.
            // Forward: s0 += in[2i+j] * h0[j]
            // Adjoint (Inverse): out[2i+j] += s0 * h0[j]
            
            // Using the adjoint property for orthogonal wavelets:
            int idx = 2 * i + j;
            if (type == SE_WAVELET_DB4) idx -= 2; // Shift for DB4 alignment
            
            while (idx < 0) idx += n;
            while (idx >= n) idx -= n;
            
            // Note: Using the same coefficients as forward for adjoint if orthogonal?
            // For orthogonal matrix W, W^-1 = W^T.
            // So we can just use the forward coefficients but accumulate into input positions.
            
            // Let's use the Transpose method which is robust for orthogonal wavelets.
            // Forward: out[i] += in[2i+j] * h0[j]
            // Inverse: out[2i+j] += in[i] * h0[j]
            
            // Re-fetching forward coefficients for Transpose implementation
            static const double h0_haar_f[] = { 0.7071067811865475, 0.7071067811865475 };
            static const double h1_haar_f[] = { -0.7071067811865475, 0.7071067811865475 };
            static const double h0_db4_f[] = { 0.4829629131445341, 0.8365163037378079, 0.2241438680420134, -0.1294095225512604 };
            static const double h1_db4_f[] = { -0.1294095225512604, -0.2241438680420134, 0.8365163037378079, -0.4829629131445341 };
            
            const double* f0 = (type == SE_WAVELET_DB4) ? h0_db4_f : h0_haar_f;
            const double* f1 = (type == SE_WAVELET_DB4) ? h1_db4_f : h1_haar_f;
            
            out[idx] += ca * f0[j];
            out[idx] += cd * f1[j];
        }
    }
}

void se_dwt(int n, const double* input, double* output, int level, se_wavelet_type_t type) {
    if (n <= 0 || !input || !output || level < 1) return;
    
    // Copy input to output initially
    std::vector<double> buf(n);
    for(int i=0; i<n; ++i) buf[i] = input[i];
    
    int current_n = n;
    for (int l = 0; l < level; ++l) {
        if (current_n % 2 != 0) break; // Must be even
        
        std::vector<double> tmp_out(current_n);
        se_dwt_step(buf.data(), tmp_out.data(), current_n, type);
        
        // Copy back to buf (approximation part will be processed next)
        for(int i=0; i<current_n; ++i) buf[i] = tmp_out[i];
        
        // Copy detail coefficients to final output
        // Detail coefficients are at the second half of tmp_out
        // They should be placed at the correct position in output
        // Structure: [A_L, D_L, D_{L-1}, ..., D_1]
        // In this loop, we generate D_{l+1} at each step.
        // Step 1: [A1, D1] -> D1 is final. A1 is processed.
        // Step 2: A1 -> [A2, D2].
        // So buf always contains the current approximation followed by details from previous steps?
        // No, buf is working buffer.
        
        // Let's write directly to output?
        // It's easier to keep the structure in `buf` and copy to `output` at the end.
        // After step 1: buf = [A1, D1]
        // After step 2: buf = [A2, D2, D1] (conceptually)
        // Implementation:
        // se_dwt_step writes to tmp_out [A, D].
        // We copy tmp_out back to buf[0...current_n-1].
        // Next step processes buf[0...current_n/2 - 1].
        
        current_n /= 2;
    }
    
    for(int i=0; i<n; ++i) output[i] = buf[i];
}

void se_idwt(int n, const double* input, double* output, int level, se_wavelet_type_t type) {
    if (n <= 0 || !input || !output || level < 1) return;

    std::vector<double> buf(n);
    for(int i=0; i<n; ++i) buf[i] = input[i];

    // Determine the size of the smallest approximation
    int current_n = n;
    for(int l=0; l<level; ++l) current_n /= 2;
    if (current_n < 2) current_n = 2; // Safety

    // Reconstruct from smallest to largest
    // We need to start from the deepest level.
    // Deepest level approximation is at buf[0...current_n-1] ? No.
    // If level=1, n=10. A1=5, D1=5. current_n=5.
    // We process buf[0...9].
    
    // We need to reverse the loop of DWT.
    // We start with current_n (smallest size * 2) and go up to n.
    
    // int start_n = n >> (level - 1); // Size of the first reconstruction block
    // Example: n=8, level=2.
    // DWT: 8 -> 4 ([A1, D1]) -> 2 ([A2, D2]). Output: [A2, D2, D1].
    // IDWT: Start with [A2, D2] (size 2+2=4). Reconstruct A1.
    // Then [A1, D1] (size 4+4=8). Reconstruct Signal.
    
    // Correct logic:
    // The smallest approximation A_L is at buf[0 ... n/(2^L)-1]
    // The corresponding detail D_L is at buf[n/(2^L) ... n/(2^(L-1))-1]
    
    int step_size = n >> (level); // Size of A_L
    step_size *= 2; // Size of (A_L + D_L) -> A_{L-1}
    
    for (int l = 0; l < level; ++l) {
        std::vector<double> tmp_in(step_size);
        std::vector<double> tmp_out(step_size);
        
        // Copy current level coefficients to temp
        for(int i=0; i<step_size; ++i) tmp_in[i] = buf[i];
        
        se_idwt_step(tmp_in.data(), tmp_out.data(), step_size, type);
        
        // Copy reconstructed approximation back to buf
        for(int i=0; i<step_size; ++i) buf[i] = tmp_out[i];
        
        step_size *= 2;
    }

    for(int i=0; i<n; ++i) output[i] = buf[i];
}

void se_dwt(int n, const float* input, float* output, int level, se_wavelet_type_t type) {
    if (n <= 0 || !input || !output) return;
    std::vector<double> in_d(n), out_d(n);
    for(int i=0; i<n; ++i) in_d[i] = (double)input[i];
    se_dwt(n, in_d.data(), out_d.data(), level, type);
    for(int i=0; i<n; ++i) output[i] = (float)out_d[i];
}

void se_idwt(int n, const float* input, float* output, int level, se_wavelet_type_t type) {
    if (n <= 0 || !input || !output) return;
    std::vector<double> in_d(n), out_d(n);
    for(int i=0; i<n; ++i) in_d[i] = (double)input[i];
    se_idwt(n, in_d.data(), out_d.data(), level, type);
    for(int i=0; i<n; ++i) output[i] = (float)out_d[i];
}
