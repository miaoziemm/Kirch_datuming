#ifndef SE_TRACE_GB_H
#define SE_TRACE_GB_H
#include "se_gb.h"
#include "se_vel.h"
#include "se_ode_cfg.h"
#include "se_module_par_desc.h"
#include "se_ray.h"
#include "se_trace_ray.h"

/**
 * Fill a Gaussian beam propagation table. The table is assumed to be 
 * allocated by the caller. This routine will fill the table so that the
 * Gaussian beam times for each entry are between gb0->t and to_t, in 
 * increasing order. 
 *
 * \param[in] ode_cfg the ODE tracing configuration
 * \param[in] v_app the velocity approximation object
 * \param[in] bnds domain bounds for terminating the tracing
 * \param[in] gb the Gaussian beam to start tracing from
 * \param[in] to_t the time to trace the Gaussian beam to
 * \param[in] gb_tbl the preallocated table to hold the traced Gaussian beam
 * \param[out] tbl_len the length of the Gaussian beam table
 * \param[out] tbl_min_i the first valid index in the table
 * \param[out] tbl_max_i the last valid index in the table
 *
 * \returns Returns an exit code that describes how the tracing ended.
 * 
 */
ode_return_t gb_trace_tbl(const ode_cfg_t *ode_cfg, 
                                   se_vel_approx_t* v_app, 
                                   const ray_bnds_t *bnds, 
                                   const gb_t* gb, double to_t, 
                                   gb_t* gb_tbl, size_t tbl_len, 
                                   size_t *tbl_min_i, size_t *tbl_max_i);

/**
 * Set r so that |p| = 1/c;
 *
 * \param[in,out] gb the beam structure
 * \param[in] v_app the velocity approximation object
 *
 */
inline void gb_set_r_from_vel(gb_t *gb, se_vel_approx_t *v_app)
{
    ray_set_r_from_vel(&gb->ray, v_app);
}

/**
 * Ensure that |p| = 1/c(x).
 *
 * \param[in,out] gb the beam structure
 * \param[in] v_app the velocity approximation object
 *
 */
inline void gb_enforce_eikonal(gb_t *gb, se_vel_approx_t *v_app)
{
    ray_enforce_eikonal(&gb->ray, v_app);
}


/**
 * Evaluate the right hand side for the Gaussian beam ODE system.
 *
 \f{eqnarray*}{
 \dot{\vec{x}}(t) &=& c(\vec{x}) \frac{\vec{p}}{|\vec{p}|} \\
 \dot{\vec{p}}(t) &=& -|\vec{p}| \nabla c(\vec{x}) \\
 \dot{H}_{re}(t) &=&  A + H_{re} B + B^T H_{re} + H_{re}C H_{re} - H_{im} C H_{im} \\
 \dot{H}_{im}(t) &=&  H_{im} B + B^T H_{im} + H_{im}C H_{re} + H_{re} C H_{im} \\
 \dot{a}(t) &=&  a(t) \left( \frac{\vec{p}\cdot\nabla c(\vec{x}) + \frac{\vec{p}\cdot H \vec{p}}{|\vec{p}|^2} - c(\vec{x}) \rm{Tr}[H]}{2|\vec{p}|} \right)
 \f}
 * Note that \f$a(t)\f$ is complex and
 \f{eqnarray*}{
 A &=& -|\vec{p}| \nabla^2 c(\vec{x}) \\
 B &=& -\frac{\vec{p}}{|\vec{p}|} \otimes \nabla c(\vec{x}) \\
 C &=& -c(\vec{x})^2 \left( Id - \frac{\vec{p}\otimes\vec{p}}{|\vec{p}|^2}\right)
 \f}
 *
 * \param[in] ode_cfg the ode configuration
 * \param[in] v_app velocity approximation object
 * \param[in] gb Gaussian beam parameter at which to evaluate the rhs
 * \param[out] rhs the right hand side of the ODE system
 *
 */
void gb_eval_ode_rhs(const ode_cfg_t *ode_cfg,
                          se_vel_approx_t *v_app,
                          const gb_t *gb,
                          gb_t *rhs);



#endif