#ifndef SE_VCONV_H
#define SE_VCONV_H
#include "se_interpolation.h"
#include "se_type.h"
#include "se_log.h"
#include "se_basic_math.h"
#include "se_assert.h"
#include "se_grid.h"



/*
 * Note on average velocity routines: vavg(t) = z(t)/t; vavg(z)=z/t(z) 
 * Used for vertical shifting in ztwell_calib.  //
 * Also useful if a velocity is to be sub-sampled without changing event depth. sub sampled 
 */

// time to time
void vconv_time_int_to_rms(const axa_t *a, const float *vint, float *vrms);
void vconv_time_int_to_avg(const axa_t *a, const float *vint, float *vavg);
void vconv_time_rms_to_int(const axa_t *a, const float *vrms, float *vint);
void vconv_time_rms_to_avg(const axa_t *a, const float *vrms, float *vavg);
void vconv_time_avg_to_int(const axa_t *a, const float *vavg, float *vint);
void vconv_time_avg_to_rms(const axa_t *a, const float *vavg, float *vrms);

// depth to depth
void vconv_depth_int_to_rms(const axa_t *a, const float *vint, float *vrms);
void vconv_depth_int_to_avg(const axa_t *a, const float *vint, float *vavg);
void vconv_depth_rms_to_int(const axa_t *a, const float *vrms, float *vint);
void vconv_depth_rms_to_avg(const axa_t *a, const float *vrms, float *vavg);
void vconv_depth_avg_to_int(const axa_t *a, const float *vavg, float *vint);
void vconv_depth_avg_to_rms(const axa_t *a, const float *vavg, float *vrms);

// axis conversion based on a single trace
void vconv_depth_axis_from_time(const axa_t *at, const float *vint_t, axa_t *az);
void vconv_time_axis_from_depth(const axa_t *az, const float *vint_z, axa_t *at);

// time to depth
// if size of output axis exceeds that of converted input axis, extend final value and return 1
// otherwise, return 0
int vconv_time_to_depth_int_to_int(const axa_t *at, const float *vint_t, const axa_t *az, float *vint_z, float *work);
int vconv_time_to_depth_int_to_rms(const axa_t *at, const float *vint_t, const axa_t *az, float *vrms_z, float *work);
int vconv_time_to_depth_int_to_avg(const axa_t *at, const float *vint_t, const axa_t *az, float *vavg_z, float *work);
int vconv_time_to_depth_rms_to_int(const axa_t *at, const float *vrms_t, const axa_t *az, float *vint_z, float *work);
int vconv_time_to_depth_rms_to_rms(const axa_t *at, const float *vrms_t, const axa_t *az, float *vrms_z, float *work);
int vconv_time_to_depth_rms_to_avg(const axa_t *at, const float *vrms_t, const axa_t *az, float *vint_z, float *work);
int vconv_time_to_depth_avg_to_int(const axa_t *at, const float *vavg_t, const axa_t *az, float *vint_z, float *work);
int vconv_time_to_depth_avg_to_rms(const axa_t *at, const float *vavg_t, const axa_t *az, float *vrms_z, float *work);
int vconv_time_to_depth_avg_to_avg(const axa_t *at, const float *vavg_t, const axa_t *az, float *vavg_z, float *work);

// compute z(t)
int vconv_time_to_depth(const axa_t *at, const float *vint_t, float *z);

// given f(t), v(z), compute f(z)
int vconv_data_time_to_depth(const axa_t *at, const float *f_t, const float *vint_z, const axa_t *az, float *f_z, float *work);

// depth to time
// if size of output axis exceeds that of converted input axis, extend final value and return 1
// otherwise, return 0
int vconv_depth_to_time_int_to_int(const axa_t *az, const float *vint_z, const axa_t *at, float *vint_t, float *work);
int vconv_depth_to_time_int_to_rms(const axa_t *az, const float *vint_z, const axa_t *at, float *vrms_t, float *work);
int vconv_depth_to_time_int_to_avg(const axa_t *az, const float *vint_z, const axa_t *at, float *vavg_t, float *work);
int vconv_depth_to_time_rms_to_int(const axa_t *az, const float *vrms_z, const axa_t *at, float *vint_t, float *work);
int vconv_depth_to_time_rms_to_rms(const axa_t *az, const float *vrms_z, const axa_t *at, float *vrms_t, float *work);
int vconv_depth_to_time_rms_to_avg(const axa_t *az, const float *vrms_z, const axa_t *at, float *vavg_t, float *work);
int vconv_depth_to_time_avg_to_int(const axa_t *az, const float *vavg_z, const axa_t *at, float *vint_t, float *work);
int vconv_depth_to_time_avg_to_rms(const axa_t *az, const float *vavg_z, const axa_t *at, float *vrms_t, float *work);
int vconv_depth_to_time_avg_to_avg(const axa_t *az, const float *vavg_z, const axa_t *at, float *vavg_t, float *work);

// compute t(z)
int vconv_depth_to_time(const axa_t *az, const float *vint_z, float *t);

// given f(z), v(t), compute f(t)
int vconv_data_depth_to_time(const axa_t *az, const float *f_z, const float *vint_t, const axa_t *at, float *f_t, float *work);

#endif