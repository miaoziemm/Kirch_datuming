#ifndef SE_GB_H
#define SE_GB_H
#include "se_type.h"
#include "se_vel.h"
#include "se_ode_cfg.h"
#include "se_module_par_desc.h"
#include "se_ray.h"



/**
 * Structure storing a gaussian beam. Includes the complex quadratic
 * travel time apporoximation and a complex amplitude coefficent.
 * Travel time is approximated as 
 *
 * \f$ T(y) = t + p.(y-x) + .5* (y-x).H_{re}(y-x) \f$
 *  
 */
typedef struct gb_s {
    
    /* double t; /\**< travel time to the base point x *\/ */

    /* v3d_t xyz; /\**< the base point (ray point) *\/ */
    /* v3d_t pqr; /\**< linear term (slowness vector) *\/ */

    ray_t ray; //射线  走时 射线点位置 射线方向 慢度表示
    
    m3d_t H_re; /**< quadratic term of travel time */ //走时的二阶项
    m3d_t H_im; /**< Gaussian taper, should be positive definite */ //高斯衰减 应该是正定的

    double a_re; /**< real part of the amplitude coefficient */  //振幅的实部
    double a_im; /**< imaginary part of the amplitude coefficient */ //振幅的虚部

} gb_t;


/**
 * Accumulate two gaussian beam structures with a constant.
 *
 * \param[in,out] gb_a the accumulating gaussian beam structure
 * \param[in] c the multiplicative constant
 * \param[in] gb the second gaussian beam structure
 */
inline void gb_accum_scaled(gb_t *gb_a, double c, const gb_t* gb)
{
    ray_accum_scaled(&gb_a->ray, c, &gb->ray);

    /* gb_a->t += c*gb->t; */
    /* v3d_accum_scaled(&gb_a->xyz, c, &gb->xyz); */
    /* v3d_accum_scaled(&gb_a->pqr, c, &gb->pqr); */
    m3d_sym_accum_scaled(&gb_a->H_re, c, &gb->H_re);
    m3d_sym_accum_scaled(&gb_a->H_im, c, &gb->H_im);
    gb_a->a_re += c*gb->a_re;
    gb_a->a_im += c*gb->a_im;

}

/**
 * Add two gaussian beam structures with a constant.
 *
 * gb_out = a*gb_a + b*gb_b;
 *
 * \param[out] gb_out the accumulating gaussian beam structure
 * \param[in] a multiplicative constant
 * \param[in] gb_a the first gaussian beam structure
 * \param[in] b multiplicative constant
 * \param[in] gb_b the second gaussian beam structure
 */
inline void gb_add_scaled(gb_t *gb_out, 
                                   double a, const gb_t *gb_a, 
                                   double b, const gb_t *gb_b)
{

    ray_add_scaled(&gb_out->ray, a, &gb_a->ray, b, &gb_b->ray);

    /* gb_out->t = a*gb_a->t + b*gb_b->t; */

    /* v3d_add_scaled(&gb_out->xyz, a, &gb_a->xyz, b, &gb_b->xyz); */
    /* v3d_add_scaled(&gb_out->pqr, a, &gb_a->pqr, b, &gb_b->pqr); */

    m3d_sym_add_scaled(&gb_out->H_re, a, &gb_a->H_re, b, &gb_b->H_re);
    m3d_sym_add_scaled(&gb_out->H_im, a, &gb_a->H_im, b, &gb_b->H_im);

    gb_out->a_re = a*gb_a->a_re + b*gb_b->a_re;
    gb_out->a_im = a*gb_a->a_im + b*gb_b->a_im;

}


/**
 * Gaussian beam norm on ray components relative to another Gaussian beam.
 *
 * \param[in] gb the Gaussian beam
 * \param[in] gb_rel the norm is rescaled relative to this Gaussian beam
 *
 * \returns the relative norm
 */
inline double gb_ray_norm_rel(const gb_t *gb, const gb_t *gb_rel)
{

    // check position only for now ... speeds up the computation
    /* return sqrt(v3d_dot(&gb->xyz,&gb->xyz) / (v3d_dot(&gb_rel->xyz,&gb_rel->xyz)+1e-20) + */
    /*             v3d_dot(&gb->pqr,&gb->pqr) / (v3d_dot(&gb_rel->pqr,&gb_rel->pqr)+1e-20) ); */

    return ray_norm_rel(&gb->ray, &gb_rel->ray);

    /* return sqrt(v3d_dot(&gb->xyz,&gb->xyz) / (v3d_dot(&gb_rel->xyz,&gb_rel->xyz)+1e-20) ); */
}

/**
 * Print Gaussian beam components.
 *
 * \param[in] info_level INFOV level at which to display components
 * \param[in] gb the Gaussian beam
 * \param[in] desc a description of the gaussian beam
 *
 */
inline void gb_info(int verb, const gb_t *gb, char *desc)
{
    ray_info(verb, &gb->ray, desc);

    INFOV((verb, 
          "   [%1.4e, %1.4e, %1.4e]    [%1.4e, %1.4e, %1.4e]\n"
          "Re=[%1.4e, %1.4e, %1.4e] Im=[%1.4e, %1.4e, %1.4e]\n"
          "   [%1.4e, %1.4e, %1.4e]    [%1.4e, %1.4e, %1.4e]\n",
          desc,
          gb->H_re.m[0], gb->H_re.m[1], gb->H_re.m[2],
          gb->H_im.m[0], gb->H_im.m[1], gb->H_im.m[2],
          gb->H_re.m[3], gb->H_re.m[4], gb->H_re.m[5],
          gb->H_im.m[3], gb->H_im.m[4], gb->H_im.m[5],
          gb->H_re.m[6], gb->H_re.m[7], gb->H_re.m[8],
          gb->H_im.m[6], gb->H_im.m[7], gb->H_im.m[8]));
}


#endif