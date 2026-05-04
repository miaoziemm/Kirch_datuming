#ifndef SE_ALLPASS_H
#define SE_ALLPASS_H
#include <stdint.h>
#include "se_type.h"
#include "se_alloc.h"
#include "se_log.h"
#include "se_util.h"
#include "se_basic_math.h"
typedef void* se_allpass_t;

/**
 * Creates a PWD object.  创建一个PWD对象
 *
 * \param nw filter size (1,2,3,...)
 * \param nj filter step
 * \param nx,ny,nz data size
 * \param pp data [nz*ny*nx]
 */
se_allpass_t se_allpass_init(int nw,
                                 int nj,
                                 int nx, int ny, int nz,
                                 float *pp);

/**
 * Initialize the lookup table for speeding up dip calculation
 * using some precomputed values and interpolation.
 * 为了加速dip计算采用一些预计算值和插值来初始化查找表
 * \param[in,out] _ap PWD object
 * \param[in] tbl_min minimal dip to precompute
 * \param[in] tbl_max maximal dip to precompute
 * \param[in] tbl_n number of dips to precompute
 *
 */
void se_allpass_init_lookup_table(se_allpass_t _ap,
                                    double tbl_min, double tbl_max, int tbl_n);


/**
 * Apply a dip dependent band pass filter to the PWD coefficients in
 * PWD lookup table for anti-aliasing. Since the PWD filter is a
 * convolution filter, the band pass filter on the coefficients acts
 * as if it it is applied to the traces.
 *
 * \param[in,out] _ap PWD object
 *
 */
void se_allpass_lookup_table_init_antialias(se_allpass_t _ap);

/**
 * Deallocate a PWD object 
 */
void se_allpass_destroy(se_allpass_t _ap);

/**
 * In-line plane-wave destruction.
 *
 * \param _ap PWD object
 * \param left indicates left or right prediction
 * \param der derivative flag
 * \param xx input
 * \param yy output
 */
void se_allpass1(const se_allpass_t _ap, 
                   const int left,
                   const int der, 
                   const float* xx, 
                   float* yy);

/**
 * Cross-line plane-wave destruction 
 *
 * \param _ap PWD object
 * \param left indicates left or right prediction
 * \param der derivative flag
 * \param xx input
 * \param yy output
 */
void se_allpass2(const se_allpass_t _ap, 
                   const int left,
                   const int der, 
                   const float* xx, 
                   float* yy);

/**
 * Lookup the PWD filter coefficients from the table. If the
 * dip is outside of the table, the PWD coefficients will be 
 * computed exactly, otherwise they will be linearly interpolated
 * from the table.
 *
 * \param[in] _ap PWD opject
 * \param[in] p the dip
 * \param[out] flt output array for the coefficients
 * \param[out] min_i smallest index of a non-negligible entry in the filter
 * \param[out] max_i largest index of a non-negligible entry in the filter
 *
 */ 
void se_allpass_filter_lookup(se_allpass_t _ap, double p, float *flt, int *min_i, int *max_i);

/**
 * Shift a part of a trace in time using pwd filter.
 *
 * \param _ap PWD opject
 * \param shift amount to shift the trace by
 * \param tr the seismic trace
 * \param n1 number of samples in the trace
 * \param i1 center of part that will be shifted
 * \param w1 length of the part that will be shifted (in samples)
 * \param tr_shift size=w1 - the result of the shift
 */ 
void se_allpass_shift_trace(se_allpass_t _ap,
                              double shift,
                              const float *tr, 
                              int n1, 
                              int i1, 
                              int w1, 
                              double *tr_shift);
/**
 * \param _ap PWD object
 * \param left indicates left or right prediction
 * \param der derivative flag
 * \param u input
 * \param beam_p
 * \param beam_crd
 * \param beam_halfwid
 * \param beam_total
 * \param beam_yy output
 */
void se_beam_allpass1_L2(const se_allpass_t _ap,
                           const int left,
                           const int der,
                           const float* u,
                           const float* beam_p,
                           const int* beam_crd,
                           const int* beam_halfwid,
                           const int beam_total,
                           float* beam_yy);



#endif