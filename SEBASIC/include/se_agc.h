#ifndef SE_AGC_H
#define SE_AGC_H
#include <stdint.h>
#include "se_type.h"
#include "se_alloc.h"
#include "se_log.h"
#include "se_util.h"
#include "se_basic_math.h"


typedef struct se_agc_s se_agc_t;

se_agc_t* se_agc_init(int window, int dwind, int detect, real threshold);

se_agc_t* se_agc_clone(const se_agc_t* src);

void se_agc_destroy( se_agc_t* agc );

void se_apply_agc( se_agc_t* agc, reala* data, int n );

void se_agc_get_pars( const se_agc_t* agc, int *window, int *dwind, int *detect, real *threshold);


#endif