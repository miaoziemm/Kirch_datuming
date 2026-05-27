#ifndef BF1D_H
#define BF1D_H

#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct BFAcaFactor BFAcaFactor;

BFAcaFactor *bf1d_aca_factor_create_phase_amp(int n,
                                              float **tau_mat,
                                              const fftwf_complex *Amp_mat,
                                              float omega,
                                              int p,
                                              int leaf_n);

void bf1d_aca_factor_apply(const BFAcaFactor *F,
                           const fftwf_complex *Uin_is,
                           fftwf_complex *Uout_is);

void bf1d_aca_factor_destroy(BFAcaFactor *F);
void butterfly_apply_1d_phase_amp(int n_org,
                                  float **tau_mat,
                                  const fftwf_complex *Amp_mat,
                                  const fftwf_complex *Uin_is,
                                  fftwf_complex *Uout_is,
                                  float omega,
                                  int p,
                                  int leaf_n);
void butterfly_apply_1d_phase_amp_main_tail(int n,
                                            float **tau_mat,
                                            const fftwf_complex *Amp_mat,
                                            const fftwf_complex *Uin_is,
                                            fftwf_complex *Uout_is,
                                            float omega,
                                            int p,
                                            int leaf_n);
                                            

#endif