#ifndef SE_TRACE_RAY_H
#define SE_TRACE_RAY_H
#include "se_type.h"
#include "se_vel.h"
#include "se_rsf.h"
#include "se_ray.h"
#include "se_interpolation.h"
#include "se_ode_cfg.h"


/**
 * Ray tracing domain bounds.
 */
typedef struct ray_bnds_s {

    v3d_t min_xyz;
    v3d_t max_xyz;

    rsf2d_t* topo;

} ray_bnds_t;

void ray_bnds_init(ray_bnds_t *bnds, 
                        double min_x, double max_x, 
                        double min_y, double max_y, 
                        double min_z, double max_z,
                        rsf2d_t* topo); 

int ray_is_in_bnds(const ray_bnds_t *bnds, 
                        const v3d_t *xyz);


/**
 * Fill a ray table. The table is assumed to be allocated by the
 * caller. This routine will fill the table so that the Gaussian beam
 * times for each entry are between ray0->t and to_t, in increasing
 * order.
 *
 * \param[in] ode_cfg the ODE tracing configuration
 * \param[in] v_app the velocity approximation object
 * \param[in] bnds domain bounds for terminating the tracing
 * \param[in] ray the ray parameters to start tracing from
 * \param[in] to_t the time to trace the ray to
 * \param[in] ray_tbl the preallocated table to hold the traced ray
 * \param[out] tbl_len the length of the ray table
 * \param[out] tbl_min_i the first valid index in the table
 * \param[out] tbl_max_i the last valid index in the table
 *
 * \returns Returns an exit code that describes how the tracing ended.
 * 
 */
ode_return_t ray_trace_tbl(const ode_cfg_t *ode_cfg, 
                                   se_vel_approx_t* v_app, 
                                   const ray_bnds_t *bnds, 
                                   const ray_t* ray, double to_t, 
                                   ray_t* ray_tbl, size_t tbl_len, 
                                   size_t *tbl_min_i, size_t *tbl_max_i);


/**
 * Set r so that |p| = 1/c;
 *
 * \param[in,out] ray the beam structure
 * \param[in] v_app the velocity approximation object
 *
 */
inline void ray_set_r_from_vel(ray_t *ray, se_vel_approx_t *v_app)
{
    double c, pq, over_c;
    int sgn;

    se_vel_get_vel_at_pt(v_app, &ray->xyz, &ray->pqr, &c);
    over_c = 1.0/(c);

    sgn = sign(ray->pqr.v[2]);

    ray->pqr.v[2]=0;
    
    pq = v3d_norm(&ray->pqr);
    if (pq<over_c) {
        ray->pqr.v[2] = sgn*se_sqrt(over_c*over_c - pq*pq);
    } else {
        v3d_scale(&ray->pqr,pq*over_c);
    }

}


/**
 * Ensure that |p| = 1/c(x).
 *
 * \param[in,out] ray the beam structure
 * \param[in] v_app the velocity approximation object
 *
 */
inline void ray_enforce_eikonal(ray_t *ray, se_vel_approx_t *v_app)
{
    double c;

    se_vel_get_vel_at_pt(v_app, &ray->xyz, &ray->pqr, &c);
    v3d_set_norm(&ray->pqr,1/c);

}

/**
 * Evaluate the right hand side for the ray tracing ODE system.
 *
 \f{eqnarray*}{
 \dot{\vec{x}}(t) &=& c(\vec{x}) \frac{\vec{p}}{|\vec{p}|} \\
 \dot{\vec{p}}(t) &=& -|\vec{p}| \nabla c(\vec{x}) \\
 *
 * \param[in] ode_cfg the ode configuration
 * \param[in] v_app velocity approximation object
 * \param[in] ray Gaussian beam parameter at which to evaluate the rhs
 * \param[out] rhs the right hand side of the ODE system
 *
 */
void ray_eval_ode_rhs(const ode_cfg_t *ode_cfg,
                           se_vel_approx_t *v_app,
                           const ray_t *ray,
                           ray_t *rhs);

/**
 * Evaluate the right hand side for the ray tracing ODE system from
 * already computed velocity gradients.
 *
 */
void ray_eval_ode_rhs_and_vel(const ode_cfg_t *ode_cfg,
                                   se_vel_approx_t *v_app,
                                   const ray_t *ray,
                                   ray_t *rhs,
                                   double *v, 
                                   v3d_t *gx_v,
                                   v3d_t *gp_v,
                                   m3d_t *hv);



#endif