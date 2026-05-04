#ifndef SE_ODE_CFG_H
#define SE_ODE_CFG_H
#include "se_module_par_desc.h"



/**
 * ODE tracers.
 */
typedef enum ode_tracer_type {
    ODE_TYPE_RK23,   /**< Runge-Kutta 23 ODE tracer. */
    ODE_TYPE_EULER,  /**< Euler ODE tracer. */
    ODE_TYPE_UNKNOWN  /**< Unknown ODE tracer. */
} ode_tracer_type;

/**
 * ODE tracing exit codes.
 */
typedef enum ode_return_t {
    ODE_SUCCESS, /**< Tracing finished successfully */
    ODE_OUTSIDE, /**< Tracing reached edge of velocity before finishing */
    ODE_ATBNDRY, /**< Tracing reached a boundary */
    ODE_UNKNOWN  /**< All other errors */
} ode_return_t;


/**
 * ODE tracing configuration structure.
 */
typedef struct ode_cfg_s {

    int solver_type_id; /**< Solver type ID used to specify the solver type */
    ode_tracer_type solver_type; /**< Solver type */  //求解方法

    // for rk23
    double min_dt; /**< Smallest tracing time step for adaptive tracing */
    double max_dt; /**< Largest tracing time step for adaptive tracing */
    double min_err_tol; /**< Threshhold for increasing the time step for adaptive tracing */ 
    double max_err_tol; /**< Threshhold for decreasing the time step for adaptive tracing */ 

    // for euler
    double fixed_dt; /**< tracing time step for fixed tracing */

    // for wf near salt
    int use_wf; /**< Flag for wave front approximation for ray propagation */ //是否使用波前近似
    double wf_threshold; /**< Wavefront velocity threshold */  //波前速度的阈值
    double wf_size; /**< Wavefront patch size */  //波前片的大小
    int wf_sets; /**< Number of concentric wave front patches */  //中心波前片的个数
    double wf_weights[2]; /**< Coefficients for avraging the edge points and middle point of the wave front */ 

    se_module_par_desc_t *par_desc; 

} ode_cfg_t;


inline int ode_cfg_solver_type_id(ode_tracer_type solver_type)
{
    switch(solver_type) {
    case ODE_TYPE_RK23:
        return 0;
        break;
    case ODE_TYPE_EULER:
        return 1;
        break;
    default:
    case ODE_TYPE_UNKNOWN:
        return -1;
        break;        
    }
}
inline ode_tracer_type ode_cfg_solver_type_from_id(int type_id)
{
    switch(type_id) {
    case 0:
        return ODE_TYPE_RK23;
        break;
    case 1: 
        return ODE_TYPE_EULER;
        break;
    default:
        return ODE_TYPE_UNKNOWN;
        break;        
    }
}

void ode_cfg_set_rk23_opts(ode_cfg_t *cfg,
                               double min_err_tol, double max_err_tol, 
                               double min_dt, double max_dt);
void ode_cfg_set_euler_opts(ode_cfg_t *cfg, 
                                double dt);
void ode_cfg_set_wf_opts(ode_cfg_t *cfg, 
                             int use_wf, double wf_threshold, double wf_size, int wf_sets,
                             double center_pt_weight);

ode_cfg_t *ode_cfg_create(void);

void ode_cfg_get_def(ode_cfg_t* p);

int ode_cfgs_are_equiv(const ode_cfg_t *A, const ode_cfg_t *B);

/**
 * Validate the ODE parameters. 
 *
 * \param cfg structure containing the ode configuration
 */
void ode_cfg_validate(const ode_cfg_t *p);

/**
 * Report the ODE tracing parameters to the user. 
 *
 * \param cfg structure containing the ode configuration
 */
void ode_cfg_report(const ode_cfg_t *cfg);

/**
 * Free memory used by the ODE tracing parameters.
 *
 * \param[in,out] p ODE configuration structure to de-allocate
 */
// Initialize a stack-allocated ode_cfg_t. Allocates and sets up par_desc.
// Use this when ode_cfg_t is embedded in another struct.
void ode_cfg_init(ode_cfg_t* p);

inline void ode_cfg_cleanup(ode_cfg_t* p)
{
    if (p) {
        se_module_par_desc_destroy(p->par_desc);
        p->par_desc = NULL;
    }
}

inline void ode_cfg_destroy(ode_cfg_t* p)
{
    if (p) {
        ode_cfg_cleanup(p);
        free(p);
    }
}




#endif