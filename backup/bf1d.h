#ifndef BF1D_H
#define BF1D_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Adaptive segmented multi-level 1D butterfly factorization with mask-aware exact leaf-stage evaluation for
 *     K(row,col) = Amp(row,col) * exp(-i * omega * tau(row,col)).
 *
 * p must satisfy p < leaf_n.  n must be positive.  The segmented factor
 * divides the matrix into active square panels.  Panels larger than the leaf
 * size use the true butterfly recursion; active leaf panels are evaluated by
 * a deterministic exact leaf-stage kernel to avoid interpolation across hard
 * aperture/anti-alias masks.  This is not an error-triggered fallback.
 */
typedef struct BFStrictFactor BFStrictFactor;
typedef struct BFStrictSegmentedFactor BFStrictSegmentedFactor;

BFStrictFactor *bf1d_strict_factor_create_phase_amp(
    int n,
    float **tau_mat,
    const fftwf_complex *Amp_mat,
    float omega,
    int p,
    int leaf_n);

void bf1d_strict_factor_apply(
    const BFStrictFactor *F,
    const fftwf_complex *Uin,
    fftwf_complex *Uout);

void bf1d_strict_factor_destroy(BFStrictFactor *F);

BFStrictSegmentedFactor *bf1d_strict_segmented_create_phase_amp(
    int n,
    float **tau_mat,
    const fftwf_complex *Amp_mat,
    float omega,
    int p,
    int leaf_n,
    int panel_levels,
    float amp_eps,
    float phase_tol);

void bf1d_strict_segmented_apply(
    const BFStrictSegmentedFactor *F,
    const fftwf_complex *Uin,
    fftwf_complex *Uout);

void bf1d_strict_segmented_destroy(BFStrictSegmentedFactor *F);

/* phase_tol controls when a fully active panel is allowed to use multi-level
 * butterfly recursion.  Smaller phase_tol improves accuracy by subdividing
 * more panels to exact leaves; larger phase_tol uses more classical butterfly
 * panels and is faster.  A good default for the supplied tests is 2.0.
 */
void butterfly_apply_1d_phase_amp(
    int n,
    float **tau_mat,
    const fftwf_complex *Amp_mat,
    const fftwf_complex *Uin,
    fftwf_complex *Uout,
    float omega,
    int p,
    int leaf_n);

#ifdef __cplusplus
}
#endif
                                            

#endif