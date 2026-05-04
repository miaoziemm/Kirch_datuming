#ifndef SE_SU_FFT_H
#define SE_SU_FFT_H

#include <complex>

/* Prime Factor FFTs */

/**
 * @brief Return valid n for complex-to-complex PFA FFT.
 * 
 * The returned n will be composed of mutually prime factors from
 * the set {2,3,4,5,7,8,9,11,13,16}.
 * 
 * @param nmin Lower bound on returned value.
 * @return Valid n for prime factor fft.
 */
int su_fft_npfa (int nmin);

/**
 * @brief Return optimal n for complex-to-complex PFA FFT.
 * 
 * The optimal n is chosen to minimize the estimated cost of performing the fft,
 * while satisfying the constraint, if possible, that n not exceed nmax.
 * 
 * @param nmin Lower bound on returned value.
 * @param nmax Desired (but not guaranteed) upper bound on returned value.
 * @return Valid n for prime factor fft.
 */
int su_fft_npfao (int nmin, int nmax);

/**
 * @brief Return valid n for real-to-complex/complex-to-real PFA FFT.
 * 
 * Requires that the transform length n be even and that n/2 be a valid length
 * for a complex-to-complex prime factor fft.
 * 
 * @param nmin Lower bound on returned value.
 * @return Valid n for real-to-complex/complex-to-real prime factor fft.
 */
int su_fft_npfar (int nmin);

/**
 * @brief Return optimal n for real-to-complex/complex-to-real PFA FFT.
 * 
 * @param nmin Lower bound on returned value.
 * @param nmax Desired (but not guaranteed) upper bound on returned value.
 * @return Valid n for real-to-complex/complex-to-real prime factor fft.
 */
int su_fft_npfaro (int nmin, int nmax);

/**
 * @brief 1D PFA complex to complex FFT.
 * 
 * @param isign Sign of exponent in fourier kernel (-1 for forward, 1 for inverse).
 * @param n Length of transform.
 * @param z Array[n] of complex numbers to be transformed in place.
 */
void su_fft_pfacc (int isign, int n, std::complex<float> z[]);

/**
 * @brief 1D PFA real to complex FFT.
 * 
 * @param isign Sign of exponent in fourier kernel.
 * @param n Length of transform; must be even.
 * @param rz Array[n] of real values (may be equivalenced to cz).
 * @param cz Output array[n/2+1] of complex values (may be equivalenced to rz).
 */
void su_fft_pfarc (int isign, int n, float rz[], std::complex<float> cz[]);

/**
 * @brief 1D PFA complex to real FFT.
 * 
 * @param isign Sign of exponent in fourier kernel.
 * @param n Length of transform.
 * @param cz Array[n/2+1] of complex values (may be equivalenced to rz).
 * @param rz Output array[n] of real values (may be equivalenced to cz).
 */
void su_fft_pfacr (int isign, int n, std::complex<float> cz[], float rz[]);

/**
 * @brief 2D PFA complex to complex FFT.
 * 
 * @param isign Sign of exponent in fourier kernel.
 * @param idim Dimension to transform, either 1 or 2.
 * @param n1 1st (fast) dimension of array.
 * @param n2 2nd (slow) dimension of array.
 * @param z Array[n2][n1] of complex elements to be transformed in place.
 */
void su_fft_pfa2cc (int isign, int idim, int n1, int n2, std::complex<float> z[]);

/**
 * @brief 2D PFA real to complex FFT.
 * 
 * @param isign Sign of exponent in fourier kernel.
 * @param idim Dimension to transform, either 1 or 2.
 * @param n1 1st (fast) dimension of array.
 * @param n2 2nd (slow) dimension of array.
 * @param rz Array of real values (may be equivalenced to cz).
 * @param cz Output array of complex values (may be equivalenced to rz).
 */
void su_fft_pfa2rc (int isign, int idim, int n1, int n2, float rz[], std::complex<float> cz[]);

/**
 * @brief 2D PFA complex to real FFT.
 * 
 * @param isign Sign of exponent in fourier kernel.
 * @param idim Dimension to transform, either 1 or 2.
 * @param n1 1st (fast) dimension of array.
 * @param n2 2nd (slow) dimension of array.
 * @param cz Array of complex values (may be equivalenced to rz).
 * @param rz Output array of real values (may be equivalenced to cz).
 */
void su_fft_pfa2cr (int isign, int idim, int n1, int n2, std::complex<float> cz[], float rz[]);

/**
 * @brief Multiple PFA complex to complex FFT.
 * 
 * @param isign Sign of exponent in fourier kernel.
 * @param n Number of complex elements per transform.
 * @param nt Number of transforms.
 * @param k Stride in complex elements within transforms.
 * @param kt Stride in complex elements between transforms.
 * @param z Array of complex elements to be transformed in place.
 */
void su_fft_pfamcc (int isign, int n, int nt, int k, int kt, std::complex<float> z[]);

#endif
