#ifndef SE_APFILT_H
#define SE_APFILT_H
#include <stdint.h>
#include "se_type.h"
#include "se_alloc.h"
#include "se_log.h"
#include "se_util.h"
#include "se_basic_math.h"
/** the PWD filter object */ 
typedef void* se_apfilt_t;

/** create a PWD filter object */
se_apfilt_t se_apfilt_init(int nw /* filter order */);

/** de-allocated a pwd object */
void se_apfilt_destroy(se_apfilt_t* _ap);

/** find filter coefficients */
void se_passfilter (se_apfilt_t _ap,
                      float p  /* slope */, 
                      float* a /* output filter [n+1] */);

/** find coefficients for filter derivative */
void se_aderfilter (se_apfilt_t _ap,
                      float p  /* slope */, 
                      float* a /* output filter [n+1] */);



#endif