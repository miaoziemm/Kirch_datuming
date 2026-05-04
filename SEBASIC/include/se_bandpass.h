#ifndef SE_BANDPASS_H
#define SE_BANDPASS_H
#include <stdint.h>
#include "se_type.h"
#include "se_alloc.h"
#include "se_log.h"
#include "se_util.h"
#include "se_basic_math.h"



typedef struct se_bandpass_s se_bandpass_t;

se_bandpass_t* se_bandpass_init( const double dt, 
                                   const double flo, const double fhi, 
                                   const int nplo, const int nphi, const int phase );

se_bandpass_t* se_bandpass_clone(const se_bandpass_t* src);

void se_bandpass_destroy( se_bandpass_t* bp );

void se_apply_bandpass( se_bandpass_t* bp, reala* data, int nsamples, int ntraces );

void se_bandpass_get_pars( const se_bandpass_t* bp, double* dt, 
                            double* flo, double* fhi, 
                            int* nplo, int* nphi, int* phase);

#endif