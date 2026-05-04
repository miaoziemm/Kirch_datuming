#ifndef SE_SPARSE_INTERPOLATION_H
#define SE_SPARSE_INTERPOLATION_H
#include "se_type.h"
#include "se_rsf.h"
#include "se_array_bit.h"

typedef void (*sparse_interpolator)( rsf3d_t *rsf, 
                                         float no_value, float no_value_tolerance,
                                         double dpow, int weight_only,
                                         double radius);

#define SE_SPINT_TYPE_FILL_HOLES 1
void sparse_fill_holes_interpolator( rsf3d_t *rsf, 
                                         float no_value, float no_value_tolerance,
                                         double dpow, int weight_only,
                                         double radius);

#define SE_SPINT_TYPE_LINEAR 2
void sparse_linear_interpolator( rsf3d_t *rsf, 
                                     float no_value, float no_value_tolerance,
                                     double dpow, int weight_only,
                                     double radius);

#define SE_SPINT_TYPE_ONEHOLE 3
void sparse_onehole_interpolator( rsf3d_t *rsf, 
                                      float no_value, float no_value_tolerance,
                                      double dpow, int weight_only,
                                      double radius);

#define SE_SPINT_TYPE_NN 4
void sparse_nn_interpolator( rsf3d_t *rsf, 
                                 float no_value, float no_value_tolerance,
                                 double dpow, int weight_only,
                                 double radius);

void sparse_trace_linear_interpolator( float *tr, const axa_t *ax, float no_value, float no_value_tolerance);


#endif