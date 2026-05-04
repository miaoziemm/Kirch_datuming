#ifndef SE_VEL_H
#define SE_VEL_H
#include <stdint.h>
#include "se_type.h"
#include "se_assert.h"
#include "se_util.h"
#include "se_hash.h"
#include "se_array.h"
#include "se_basic_math.h"
#include "se_grid.h"

/**
 * Structure storing the velocity.
 */
typedef struct se_vel_s {
    
    float  *vel_data; /**< The velocity data. Uses packing order (z,x,y) */
    float  *del_data; /**< The delta data. Uses packing order (z,x,y) */
    float  *eta_data; /**< The eta=epsilon-delta data. Uses packing order (z,x,y) */
    float  *ttix_data; /**< The X dip of tilt (dx/dz) of the axis of symmetry for TTI Uses packing order (z,x,y) */
    float  *ttiy_data; /**< The X dip of tilt (dx/dz) of the axis of symmetry for TTI Uses packing order (z,x,y) */

    int aniso; /**< Flag for anisotropic model */
    int tti; /**< Flag for TTI */

    int is3d; /**< Flag for 3D (even if input is only 2D) */

    size_t n_xyz[3]; /**< number of samples */
    double o_xyz[3]; /**< origin of the velocity cube */
    double d_xyz[3]; /**< sampling of the velocity cube */
    double m_xyz[3]; /**< maximum of the velocity cube */
    double over_d_xyz[3]; /**< 1/dx, 1/dy, 1/dz */

} se_vel_t;

/**
 * Velocity approximation cell.
 * 
 */
typedef struct se_vel_approx_cell_s se_vel_approx_cell_t;

size_t se_vel_approx_cell_get_mem_usage(se_vel_approx_cell_t * vc);

/**
 * Velocity approximation object.
 */
typedef struct se_vel_approx_s {
    
    se_vel_t *vel; /**< the velocity model */
    int aniso; /**< Flag for anisotropic model */
    int tti; /**< Flag for TTI */
    
    se_vel_approx_cell_t *vc_curr; /**< Current local velocity cell, with derivaties, etc */
    int64_t vc_curr_lidx;
    se_vel_approx_cell_t *vc_prev; /**< The previous local velocity cell, with derivaties, etc */
    int64_t vc_prev_lidx;
    int64_t re_used;
    int64_t re_used_cache;
    ssize_t lidx_n_xyz[3];

    int is2d; /**< 2D flag */

    int st_hl; /**< half length of the interpolation stencil */
    int use_hess; /**< Flag for computing hessian, if 0, velocity hessian will be 0 */ 
    int detect_bdry; /**< Flag for detecting large velocity variations */ 
    double bdry_frac; /**< Fraction for detecting large velocity variations (max_vel-min_vel)/max_vel */ 

    void (*eval_grad)(const se_vel_approx_cell_t*, const float*, int, v3d_t*); /**< function that evalutes the first derivatives */
    void (*eval_hess)(const se_vel_approx_cell_t*, const float*, int, m3d_t*); /**< function that evalutes the second derivatives */

} se_vel_approx_t;

/**
 * Initialize the velocity, by copying the given pointers into the velocity structure.
 * For isotropic models, del_data and eta_data should be NULL.
 * 
 * \param[out] vel the velocity to initilize
 * \param[in] vel_data the velocity data
 * \param[in] del_data the delta data
 * \param[in] eta_data the eta data
 * \param[in] ttix_data the tilt from verical of the symmetric axis in X in radiants
 * \param[in] ttiy_data the tilt from verical of the symmetric axis in Y in radiants
 * \param[in] velgrid the velocity grid specification
 * \param[in] force_3d if velocity data is only 2d, extend it to 3d.
 * 
 */
void se_vel_init(se_vel_t *vel, float *vel_data, 
                  float *del_data, float *eta_data,
                  float *ttix_data, float *ttiy_data,
                  const grid3d_t *velgrid,
                  int force_3d);

size_t se_vel_get_mem_usage(se_vel_t *vel);

/**
 * De-allocate the velocity data.
 * 
 * \param[in] vel the velocity to deallocate
 * 
 */
void se_vel_destroy(se_vel_t *vel);

/**
 * Initialize the velocity approximation object.
 * 
 * \param[out] v_app the velocity approximation object to initialize
 * \param[in] vel the velocity
 * \param[in] stencil the stencil length. 
 * \param[in] use_hess flag if hessian is to be computed.
 * \param[in] dtct_bdry flag controlling whether to detect large variations of velocity
 * \param[in] bdry_frac cutoff fraction for detecting large velocity variations 
 * 
 */
void se_vel_approx_init(se_vel_approx_t *v_app, 
                         se_vel_t *vel,
                         int stencil, int use_hess,
                         int dtct_bdry, double bdry_frac);

size_t se_vel_approx_get_mem_usage(const se_vel_approx_t *v_app);

/**
 * De-allocate the velocity approximation object.
 * 
 * \param[in] v_app the velocity approximation object to deallocate
 * 
 */
void se_vel_approx_destroy(se_vel_approx_t *v_app);

/**
 * Test if a point is inside the velocity grid.
 * 
 * \param[in] v_app the velocity approximation object
 * \param[in] pt the point
 * 
 * \returns 0 if point is outside, 1 if inside
 *
 */
inline int se_vel_pt_is_inside(const se_vel_approx_t *v_app, const v3d_t *pt)
{

    if ( pt->v[0] - v_app->vel->m_xyz[0] > v_app->vel->d_xyz[0] ) return 0;
    if ( v_app->vel->o_xyz[0] - pt->v[0] > v_app->vel->d_xyz[0] ) return 0;
    if (pt->v[0] < v_app->vel->o_xyz[0]) return 0;

    if ( pt->v[2] - v_app->vel->m_xyz[2] > v_app->vel->d_xyz[2] ) return 0;
    if ( v_app->vel->o_xyz[2] - pt->v[2] > v_app->vel->d_xyz[2] ) return 0;
    
    if (v_app->is2d) return 1;

    if ( pt->v[1] - v_app->vel->m_xyz[1] > v_app->vel->d_xyz[1] ) return 0;
    if ( v_app->vel->o_xyz[1] - pt->v[1] > v_app->vel->d_xyz[1] ) return 0;

    return 1;

}

/**
 * Test if a point is inside the approximation domain. For these points there is 
 * enough data for the approximation stencil that is being used.
 * 
 * \param[in] v_app the velocity approximation object
 * \param[in] pt the point
 * 
 * \returns 0 if point is outside, 1 if inside
 *
 */
int se_vel_pt_is_inside_approx_domain(const se_vel_approx_t *v_app, const v3d_t *pt);

/**
 * Get the velocity, its gradient and hessian at a point. The
 * velocity, gradient and hessian can be specified as NULL.
 * 
 * \param[in] v_app the velocity approximation object
 * \param[in] pt the point
 * \param[in] dir direction for anisotropic velocity
 * \param[out] v the velocity at the point 
 * \param[out] gx_v the X gradient of the velocity
 * \param[out] gp_v the P gradient of the velocity (anisotropy only)
 * \param[out] hv the hessian matrix of the velocity
 * 
 * \returns 0 if point is outside, 1 if inside. Gradient vector is not update if 
 * the point is outside.
 *
 */
int se_vel_get_vel_grad_hess_at_pt(se_vel_approx_t *v_app, const v3d_t *pt, const v3d_t *dir,
                                    double *v, v3d_t *gx_v, v3d_t *gp_v, m3d_t *hv);

/**
 * Get the velocity at a point.
 * 
 * \param[in] v_app the velocity approximation object
 * \param[in] pt the point
 * \param[in] dir direction for anisotropic velocity
 * \param[out] v the value of the velocity
 * 
 * \returns 0 if point is outside, 1 if inside. Velocity is not updated if point is outside.
 *
 */
inline int se_vel_get_vel_at_pt(se_vel_approx_t *v_app, const v3d_t *pt, const v3d_t *dir, double *v)
{
    return se_vel_get_vel_grad_hess_at_pt(v_app, pt, dir, v, NULL, NULL, NULL);
}

/**
 * Get the velocity gradient at a point.
 * 
 * \param[in] v_app the velocity approximation object
 * \param[in] pt the point
 * \param[in] dir direction for anisotropic velocity
 * \param[out] gv the gradient of the velocity
 * 
 * \returns 0 if point is outside, 1 if inside. Gradient vector is not update if 
 * the point is outside.
 *
 */
inline int se_vel_get_grad_x_at_pt(se_vel_approx_t *v_app, const v3d_t *pt, const v3d_t *dir, v3d_t *gv)
{
    return se_vel_get_vel_grad_hess_at_pt(v_app, pt, dir, NULL, gv, NULL, NULL);
}

/**
 * Get the velocity hessian at a point.
 * 
 * \param[in] v_app the velocity approximation object
 * \param[in] pt the point
 * \param[in] dir direction for anisotropic velocity
 * \param[out] hv the hessian matrix of the velocity
 * 
 * \returns 0 if point is outside, 1 if inside. Hessian matrix is not update if 
 * the point is outside.
 *
 */
inline int se_vel_get_hess_at_pt(se_vel_approx_t *v_app, const v3d_t *pt, const v3d_t *dir, m3d_t *hv)
{
    return se_vel_get_vel_grad_hess_at_pt(v_app, pt, dir, NULL, NULL, NULL, hv);
}

/**
 * Get the velocity change fraction from the gradient. Formula is 
 * \f$\nabla \vec{v} \cdot \vec{d} / c \f$
 * 
 * \param[in] v_app the velocity approximation object
 * \param[in] pt the point
 * \param[in] dir direction for anisotropic velocity
 *
 * \return the value of the ratio
 */
double se_vel_get_vel_change_frac_from_grad(se_vel_approx_t *v_app, const v3d_t *pt, const v3d_t *dir);

int se_vel_has_rapid_change(se_vel_approx_t *v_app, const v3d_t *pt);

int se_vel_dist_to_rapid_change(se_vel_approx_t *v_app, const v3d_t *pt, double *dist, v3d_t *nv);

inline void se_vel_get_grid(const se_vel_t *vel, grid3d_t *grid)
{
    set_grid3d(grid,
                   vel->n_xyz[0], vel->d_xyz[0], vel->o_xyz[0],
                   vel->n_xyz[1], vel->d_xyz[1], vel->o_xyz[1],
                   vel->n_xyz[2], vel->d_xyz[2], vel->o_xyz[2]);
}



#endif