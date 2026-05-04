#ifndef SE_RAY_H
#define SE_RAY_H
#include "se_basic_math.h"
#include "se_type.h"
#include "se_log.h"
#include "se_ode_cfg.h"

/**
 * Structure storing a ray.
 */
typedef struct ray_s {
    
    double t; /**< ray time */  //时间

    v3d_t xyz; /**< ray point */ //射线点  vector 3D
    v3d_t pqr; /**< slowness vector */ //慢度矢量
    
} ray_t;


/**
 * Accumulate two ray structures with a constant.
 *
 * \param[in,out] ray_a the accumulating ray structure
 * \param[in] c the multiplicative constant
 * \param[in] ray the second ray structure
 */
inline void ray_accum_scaled(ray_t *ray_a, double c, const ray_t* ray)
{
    ray_a->t += c*ray->t;
    v3d_accum_scaled(&ray_a->xyz, c, &ray->xyz);
    v3d_accum_scaled(&ray_a->pqr, c, &ray->pqr);

}

/**
 * Add two ray structures with a constant.
 *
 * ray_out = a*ray_a + b*ray_b;
 *
 * \param[out] ray_out the accumulating ray structure
 * \param[in] a multiplicative constant
 * \param[in] ray_a the first ray structure
 * \param[in] b multiplicative constant
 * \param[in] ray_b the second ray structure
 */
inline void ray_add_scaled(ray_t *ray_out, 
                                    double a, const ray_t *ray_a, 
                                    double b, const ray_t *ray_b)
{
    ray_out->t = a*ray_a->t + b*ray_b->t;

    v3d_add_scaled(&ray_out->xyz, a, &ray_a->xyz, b, &ray_b->xyz);
    v3d_add_scaled(&ray_out->pqr, a, &ray_a->pqr, b, &ray_b->pqr);

}


/**
 * Ray norm on ray components relative to another ray.
 *
 * \param[in] ray the ray
 * \param[in] ray_rel the norm is rescaled relative to this ray
 *
 * \returns the relative norm
 */
inline double ray_norm_rel(const ray_t *ray, const ray_t *ray_rel)
{
    return sqrt(v3d_dot(&ray->xyz,&ray->xyz) / (v3d_dot(&ray_rel->xyz,&ray_rel->xyz)+1e-20) );
}

double ray_tbl_get_average_slowness(const ray_t* ray_tbl, 
                                         size_t tbl_min_i, 
                                         size_t tbl_max_i);
double ray_tbl_get_time(const ray_t* ray_tbl, 
                            size_t tbl_min_i, 
                            size_t tbl_max_i);
double ray_tbl_get_length(const ray_t* ray_tbl, 
                              size_t tbl_min_i, 
                              size_t tbl_max_i);
double ray_tbl_get_max_vel(const ray_t* ray_tbl, 
                               size_t tbl_min_i, 
                               size_t tbl_max_i);
double ray_get_offset(const ray_t* ray_A,
                          const ray_t* ray_B);
double ray_tbl_get_slowness(const ray_t* ray_tbl, 
                                size_t tbl_i);
double ray_tbl_get_z(const ray_t* ray_tbl, 
                         size_t tbl_i);
void ray_tbl_get_min_max_velocity(const ray_t* ray_tbl, 
                                      size_t tbl_min_i, 
                                      size_t tbl_max_i,
                                      double *min_vel, double *max_vel);
double ray_tbl_get_cos_to_v3d(const ray_t* ray_tbl, 
                                  size_t tbl_i,
                                  const v3d_t *v);

/**
 * Print ray components.
 *
 * \param[in] info_level INFOV level at which to display components
 * \param[in] ray the ray
 * \param[in] desc a description of the ray
 *
 */
inline void ray_info(int verb, const ray_t *ray, char *desc)
{
    INFOV((verb,
           "%s: t = %g xyz->(%g,%g,%g) pqr->(%g,%g,%g)\n",
           desc,
           ray->t, 
           ray->xyz.v[0], ray->xyz.v[1], ray->xyz.v[2],
           ray->pqr.v[0], ray->pqr.v[1], ray->pqr.v[2]));
}



#endif