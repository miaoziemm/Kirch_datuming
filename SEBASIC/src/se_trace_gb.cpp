#include "../include/se_trace_gb.h"


inline double get_gb_t(const gb_t *gb) { return gb->ray.t; }
inline const v3d_t* get_gb_xyz(const gb_t *gb) { return &gb->ray.xyz; }

#define ode_trace_tbl_fct_name gb_trace_tbl
#define ode_U_t gb_t
#define ode_U_add_scaled gb_add_scaled
#define ode_U_accum_scaled gb_accum_scaled
#define ode_U_enforce_constraints gb_enforce_eikonal
#define ode_eval_rhs gb_eval_ode_rhs
#define ode_U_norm_rel gb_ray_norm_rel
#define ode_U_bnds_t ray_bnds_t
#define ode_U_is_in_bnds ray_is_in_bnds
#define ode_U_get_t get_gb_t
#define ode_U_get_xyz get_gb_xyz

#include "se_ode_impl.cpp"


void gb_eval_ode_rhs(const ode_cfg_t *ode_cfg,
                          se_vel_approx_t *v_app,
                          const gb_t *gb,
                          gb_t *rhs)
{

    double v;
    v3d_t gx_v, gp_v, pqr_hat;
    m3d_t hv, B, C, tmp_M1, tmp_M2;
    double pqr_norm;
    double re, im;


    ray_eval_ode_rhs_and_vel(ode_cfg, v_app, &gb->ray, 
                                  &rhs->ray, &v, &gx_v, &gp_v, &hv);

    pqr_norm = v3d_norm(&gb->ray.pqr);
    v3d_assign_scaled(&pqr_hat,1/pqr_norm,&gb->ray.pqr);

    /* B */
    v3d_outer(&B, &pqr_hat, &gx_v);
    m3d_scale(&B,-1);


    /* ======= Mat 3d implementation ======== */

    /* /\* C *\/ */
    /* m3d_t Id; */
    /* m3d_diag_init(&Id, 1, 1, 1); */
    /* v3d_outer(&C, &pqr_hat, &pqr_hat); */
    /* m3d_sym_subtract(&C, &Id, &C); */
    /* m3d_sym_scale(&C, -v*v); */

    /* m3d_t B_T; */
    /* m3d_transp(&B_T, &B); */
    
    /* /\* H_re_dot = H_re * B + B^T * H_re*\/ */
    /* m3d_mult(&rhs->H_re, &gb->H_re, &B); */
    /* m3d_mult(&tmp_M1, &B_T, &gb->H_re); */
    /* m3d_sym_accum(&rhs->H_re, &tmp_M1); */
    /* /\* H_im_dot = H_im * B  +  B^T * H_im *\/ */
    /* m3d_mult(&tmp_M1, &gb->H_im, &B); */
    /* m3d_mult(&tmp_M2, &B_T, &gb->H_im); */
    /* m3d_sym_add(&rhs->H_im, &tmp_M1, &tmp_M2); */

    /* if (v_app->use_hess) { */
    /*     m3d_t A; */
    /*     /\* H_re_dot += A *\/ */
    /*     m3d_assign_scaled(&A, -pqr_norm, &hv); */
    /*     m3d_sym_accum(&rhs->H_re, &A); */
    /* } */

    /* m3d_mult(&tmp_M1, &C, &gb->H_re); */
    /* /\* H_re_dot += H_re * C * H_re *\/ */
    /* m3d_mult(&tmp_M2, &gb->H_re, &tmp_M1); */
    /* m3d_accum(&rhs->H_re, &tmp_M2); */
    /* /\* H_im_dot += 2*H_im * C * H_re *\/ */
    /* m3d_mult(&tmp_M2, &gb->H_im, &tmp_M1); */
    /* m3d_scale(&tmp_M2,2); */
    /* m3d_accum(&rhs->H_im, &tmp_M2); */

    /* m3d_mult(&tmp_M1, &C, &gb->H_im); */
    /* /\* H_re_dot -= H_im * C * H_im *\/ */
    /* m3d_mult(&tmp_M2, &gb->H_im, &tmp_M1); */
    /* m3d_subtract(&rhs->H_re, &rhs->H_re, &tmp_M2); */

    /* ======= Mat 3d end ======== */


    
    /* ======= Direct implementation ======== */

    /* C */
    v3d_outer(&C, &pqr_hat, &pqr_hat);
    m3d_sym_scale(&C, v*v);
    C.m[0]-=v*v;
    C.m[4]-=v*v;
    C.m[8]-=v*v;

    /* H_re_dot = H_re * B + B^T * H_re */

    rhs->H_re.m[0] = 2*( gb->H_re.m[0]*B.m[0] + gb->H_re.m[1]*B.m[3] + gb->H_re.m[2]*B.m[6] );
    rhs->H_re.m[1] = ( gb->H_re.m[0]*B.m[1] + gb->H_re.m[1]*B.m[4] + gb->H_re.m[2]*B.m[7] +
                       B.m[0]*gb->H_re.m[1] + B.m[3]*gb->H_re.m[4] + B.m[6]*gb->H_re.m[7] );
    rhs->H_re.m[2] = ( gb->H_re.m[0]*B.m[2] + gb->H_re.m[1]*B.m[5] + gb->H_re.m[2]*B.m[8] +
                       B.m[0]*gb->H_re.m[2] + B.m[3]*gb->H_re.m[5] + B.m[6]*gb->H_re.m[8] );
    rhs->H_re.m[4] = 2*( gb->H_re.m[1]*B.m[1] + gb->H_re.m[4]*B.m[4] + gb->H_re.m[5]*B.m[7] );
    rhs->H_re.m[5] = ( gb->H_re.m[3]*B.m[2] + gb->H_re.m[4]*B.m[5] + gb->H_re.m[5]*B.m[8] +
                       B.m[1]*gb->H_re.m[2] + B.m[4]*gb->H_re.m[5] + B.m[7]*gb->H_re.m[8] );
    rhs->H_re.m[8] = 2*( gb->H_re.m[2]*B.m[2] + gb->H_re.m[5]*B.m[5] + gb->H_re.m[8]*B.m[8] );

    m3d_sym_fill_lower(&rhs->H_re);

    if (v_app->use_hess) {
        /* H_re_dot += A */
        rhs->H_re.m[0] -= pqr_norm*hv.m[0];
        rhs->H_re.m[1] -= pqr_norm*hv.m[1];
        rhs->H_re.m[2] -= pqr_norm*hv.m[2];
        rhs->H_re.m[4] -= pqr_norm*hv.m[4];
        rhs->H_re.m[5] -= pqr_norm*hv.m[5];
        rhs->H_re.m[8] -= pqr_norm*hv.m[8];
    }


    /* H_im_dot = H_im * B + B^T * H_im */
    rhs->H_im.m[0] = 2*( gb->H_im.m[0]*B.m[0] + gb->H_im.m[1]*B.m[3] + gb->H_im.m[2]*B.m[6] );
    rhs->H_im.m[1] = ( gb->H_im.m[0]*B.m[1] + gb->H_im.m[1]*B.m[4] + gb->H_im.m[2]*B.m[7] +
                       B.m[0]*gb->H_im.m[1] + B.m[3]*gb->H_im.m[4] + B.m[6]*gb->H_im.m[7] );
    rhs->H_im.m[2] = ( gb->H_im.m[0]*B.m[2] + gb->H_im.m[1]*B.m[5] + gb->H_im.m[2]*B.m[8] +
                       B.m[0]*gb->H_im.m[2] + B.m[3]*gb->H_im.m[5] + B.m[6]*gb->H_im.m[8] );
    rhs->H_im.m[4] = 2*( gb->H_im.m[1]*B.m[1] + gb->H_im.m[4]*B.m[4] + gb->H_im.m[5]*B.m[7] );
    rhs->H_im.m[5] = ( gb->H_im.m[3]*B.m[2] + gb->H_im.m[4]*B.m[5] + gb->H_im.m[5]*B.m[8] +
                       B.m[1]*gb->H_im.m[2] + B.m[4]*gb->H_im.m[5] + B.m[7]*gb->H_im.m[8] );
    rhs->H_im.m[8] = 2*( gb->H_im.m[2]*B.m[2] + gb->H_im.m[5]*B.m[5] + gb->H_im.m[8]*B.m[8] );

    m3d_sym_fill_lower(&rhs->H_im);



    /* m3d_sym_mult(&tmp_M1, &C, &gb->H_re); */
    /* /\* H_re_dot += H_re * C * H_re *\/ */
    /* m3d_mult_sym_result(&tmp_M2, &gb->H_re, &tmp_M1); */
    /* m3d_sym_accum(&rhs->H_re, &tmp_M2); */
    /* /\* H_im_dot += 2*H_im * C * H_re *\/ */
    /* m3d_mult_sym_result(&tmp_M2, &gb->H_im, &tmp_M1); */
    /* /\* m3d_sym_scale(&tmp_M2,2); *\/ */
    /* m3d_sym_accum_scaled(&rhs->H_im, 2, &tmp_M2); */

    /* /\* H_re_dot -= H_im * C * H_im *\/ */
    /* m3d_sym_mult(&tmp_M1, &C, &gb->H_im); */
    /* m3d_mult_sym_result(&tmp_M2, &gb->H_im, &tmp_M1); */
    /* m3d_sym_subtract(&rhs->H_re, &rhs->H_re, &tmp_M2); */




    /* /\* H_re_dot += H_re * C * H_re *\/ */
    /* m3d_mult_HCH_sym(&tmp_M1, &gb->H_re, &C); */
    /* m3d_sym_accum(&rhs->H_re, &tmp_M1); */
    /* /\* H_re_dot -= H_im * C * H_im *\/ */
    /* m3d_mult_HCH_sym(&tmp_M2, &gb->H_im, &C); */
    /* m3d_sym_subtract(&rhs->H_re, &rhs->H_re, &tmp_M2); */


    
    /* /\* H_im_dot += 2*H_im * C * H_re *\/ */
    /* m3d_sym_mult(&tmp_M1, &C, &gb->H_re); */
    /* m3d_mult_sym_result(&tmp_M2, &gb->H_re, &tmp_M1); */
    /* m3d_sym_accum_scaled(&rhs->H_im, 2, &tmp_M2); */




    m3d_sym_mult(&tmp_M1, &C, &gb->H_re);
    /* H_re_dot += H_re * C * H_re */
    m3d_mult_sym_result(&tmp_M2, &gb->H_re, &tmp_M1);
    m3d_sym_accum(&rhs->H_re, &tmp_M2);
    /* H_im_dot += 2*H_im * C * H_re */
    m3d_mult_sym_result(&tmp_M2, &gb->H_im, &tmp_M1);
    /* m3d_sym_scale(&tmp_M2,2); */
    m3d_sym_accum_scaled(&rhs->H_im, 2, &tmp_M2);

    /* H_re_dot -= H_im * C * H_im */
    m3d_mult_HCH_sym(&tmp_M2, &gb->H_im, &C);
    m3d_sym_subtract(&rhs->H_re, &rhs->H_re, &tmp_M2);



    /* ======= Direct end ======== */

    re = .5 * v * ( v3d_dot(&gb->ray.pqr, &gx_v)
                    + v3d_m3d_sym_inner(&gb->H_re, &pqr_hat)
                    - v * m3d_trace(&gb->H_re) );
    im = .5 * v * ( v3d_m3d_sym_inner(&gb->H_im, &pqr_hat)
                    - v * m3d_trace(&gb->H_im) );
    rhs->a_re = gb->a_re * re - gb->a_im * im;
    rhs->a_im = gb->a_re * im + gb->a_im * re;

 
}
