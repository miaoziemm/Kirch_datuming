#ifndef SE_GEOMETRY_H
#define SE_GEOMETRY_H
#include "se_basic_math.h"


/**
 * For two given planes defined by a vector and a base point, find the
 * point that is on both planes and is closest to the two base points.
 *
 * \param[out] b the point to be determined
 * \param[in] p1 the base point for the first plane
 * \param[in] n1 the orthogonal vector for the first plane
 * \param[in] p2 the base point for the second plane
 * \param[in] n2 the orthogonal vector for the second plane
 *
 * \return 0 if point was determined, 1 if point could not be determined
 *
 */
int se_find_closest_pt_on_planes(v3d_t *b, 
                                   const v3d_t *p1, const v3d_t *n1,
                                   const v3d_t *p2, const v3d_t *n2);


/**
 * For two given paraboloids, find the point that is on both of them
 * and closest to their base points.
 *

对于两个给定的抛物面, 找到在它们两个上的且距离他们的基点最近的点.

paraboloid 抛物面 

 * \param[out] b the point to be determined
 * \param[out] F1 the value of the first paraboloid on b
 * \param[out] F2 the value of the second paraboloid on b
 * \param[in] p1 the base point for the first plane
 * \param[in] n1 the orthogonal vector to the first surface at p1
 * \param[in] M1 the curvature matrix for the first surface
 * \param[in] p2 the base point for the second plane
 * \param[in] n2 the orthogonal vector to the second surface at p2
 * \param[in] M2 the curvature matrix for the second surface
 * \param[in] tol tolerance for determining if a point is on the surface
 *
 * \return 0 if point was determined, 1 if point could not be determined
 *
 */
int se_find_closest_pt_on_paraboloids(v3d_t *b, double *F1, double *F2, 
                                        const v3d_t *p1, const v3d_t *n1, const m3d_t *M1,
                                        const v3d_t *p2, const v3d_t *n2, const m3d_t *M2,
                                        double tol);

/**
 * Evaluate a quadratic form at a given point. 
 *
 * Q(x) = a_0 + p_0.(x-x_0) + 1/2 (x-x_0).H_0.(x-x_0)
 *
 * If x0 is NULL, the base point is assumed to be (0,0,0)
 *
 * \param[in] x0 base point (see formula above)
 * \param[in] a0 scalar (see formula above)
 * \param[in] p0 linear (see formula above)
 * \param[in] H0 quadratic (see formula above)
 * \param[in] x point at which to evaluate the quadratic form
 *
 * \return value of the quadratic form
 *
 */
inline double se_eval_quadratic_form(const v3d_t *x0,
                                           double a0, const v3d_t *p0, const m3d_t *H0,
                                           const v3d_t *x)
{
    if (x0) {
        v3d_t x_m_x0;
        v3d_subtract(&x_m_x0,x,x0);
        return a0 + v3d_dot(p0,&x_m_x0)+.5*v3d_m3d_sym_inner(H0, &x_m_x0);
    } else {
        return a0 + v3d_dot(p0,x)+.5*v3d_m3d_sym_inner(H0, x);
    }

}

/**
 * Evaluate the gradient of a quadratic form at a given point. 计算导数 求解角度
 *
 * nabla Q(x) = p_0 + H_0.(x-x_0)
 *
 * If x0 is NULL, the base point is assumed to be (0,0,0)
 *
 * \param[in] x0 base point (see formula above)
 * \param[in] p0 linear (see formula above)
 * \param[in] H0 quadratic (see formula above)
 * \param[in] x point at which to evaluate the quadratic form
 * \param[out] grad the gradient of the quadratic from
 *
 *
 */
inline void se_eval_quadratic_form_grad(v3d_t *grad,
                                              const v3d_t *x0,
                                              const v3d_t *p0, const m3d_t *H0,
                                              const v3d_t *x)
{
    if (x0) {
        v3d_t x_m_x0;
        v3d_subtract(&x_m_x0,x,x0);
        m3d_v3d_mult(grad, H0, &x_m_x0);
    } else {
        m3d_v3d_mult(grad, H0, x);
    }
    v3d_accum(grad,p0);
}


#endif