#include "../include/se_trace_ray.h"


void ray_bnds_init(ray_bnds_t *bnds, 
                        double min_x, double max_x, 
                        double min_y, double max_y, 
                        double min_z, double max_z,
                        rsf2d_t* topo)
{
    v3d_init(&bnds->min_xyz, min_x, min_y, min_z);
    v3d_init(&bnds->max_xyz, max_x, max_y, max_z);
    bnds->topo = topo;
}

int ray_is_in_bnds(const ray_bnds_t *bnds, 
                        const v3d_t *xyz)
{
    
    if ( ( bnds->min_xyz.v[0] > xyz->v[0] ) || ( bnds->max_xyz.v[0] < xyz->v[0] ) ) return 0;
    if ( ( bnds->min_xyz.v[1] > xyz->v[1] ) || ( bnds->max_xyz.v[1] < xyz->v[1] ) ) return 0;
    if ( ( bnds->min_xyz.v[2] > xyz->v[2] ) || ( bnds->max_xyz.v[2] < xyz->v[2] ) ) return 0;

    if ( (bnds->topo) && 
         ( xyz->v[2] < se_rsf2d_value_binterp(xyz->v[0], xyz->v[1], bnds->topo) ) ) return 0;

    return 1; 
}

inline double get_ray_t(const ray_t *ray) { return ray->t; }
inline const v3d_t* get_ray_xyz(const ray_t *ray) { return &ray->xyz; }

#define ode_trace_tbl_fct_name ray_trace_tbl
#define ode_U_t ray_t
#define ode_U_add_scaled ray_add_scaled
#define ode_U_accum_scaled ray_accum_scaled
#define ode_U_enforce_constraints ray_enforce_eikonal
#define ode_eval_rhs ray_eval_ode_rhs
#define ode_U_norm_rel ray_norm_rel
#define ode_U_bnds_t ray_bnds_t
#define ode_U_is_in_bnds ray_is_in_bnds
#define ode_U_get_t get_ray_t
#define ode_U_get_xyz get_ray_xyz

#include "se_ode_impl.cpp"


static void get_x_on_wf(const ray_t *ray, int is2d, double wf_size,
                             v3d_t *x0, const v3d_t *p, v3d_t *p_n1, v3d_t *p_n2) {

    int i;

    *(x0+0) = ray->xyz;
    v3d_add_scaled(x0+1, 1, x0,  wf_size, p);
    v3d_add_scaled(x0+2, 1, x0, -wf_size, p);

    if (is2d) {
        v3d_get_2d_orthogonal(p_n1,p);

        v3d_add_scaled(x0+3, 1, x0,  wf_size, p_n1);
        v3d_add_scaled(x0+4, 1, x0, -wf_size, p_n1);

    } else {
        v3d_get_orthogonal_pair(p_n1,p_n2,p);
        v3d_normalize(p_n1);
        v3d_normalize(p_n2);

        v3d_add_scaled(x0+3,  wf_size, p_n1,  wf_size, p_n2);
        v3d_add_scaled(x0+4, -wf_size, p_n1,  wf_size, p_n2);
        v3d_add_scaled(x0+5,  wf_size, p_n1, -wf_size, p_n2);
        v3d_add_scaled(x0+6, -wf_size, p_n1, -wf_size, p_n2);

        for (i=3; i<7; i++) {
            v3d_add_scaled(x0+i, 1, x0, 1/M_SQRT2, x0+i);
        }

    }
}

static void get_vel_on_wf_pts(const v3d_t *x0, const v3d_t *p0, se_vel_approx_t *v_app,
                                   double ctr_pt_wght, double out_pt_wght, double *c0 ) {

    int i;
    double c;

    se_vel_get_vel_at_pt(v_app, x0, p0, c0);

    for (i=1;i<(v_app->is2d?5:7);i++) {
        se_vel_get_vel_at_pt(v_app, x0+i, p0, &c);
        c0[i] = ctr_pt_wght*c0[0] + out_pt_wght*c;
    }

}


void ray_eval_ode_rhs(const ode_cfg_t *ode_cfg,
                           se_vel_approx_t *v_app,
                           const ray_t *ray,
                           ray_t *rhs)
{
    double v;
    v3d_t gx_v, gp_v;

    ray_eval_ode_rhs_and_vel(ode_cfg, v_app, ray,
                                  rhs, &v, &gx_v, &gp_v, NULL);

}

void ray_eval_ode_rhs_and_vel(const ode_cfg_t *ode_cfg,
                                   se_vel_approx_t *v_app,
                                   const ray_t *ray,
                                   ray_t *rhs,
                                   double *v, 
                                   v3d_t *gx_v,
                                   v3d_t *gp_v,
                                   m3d_t *hv)
{

    double pqr_norm;

    rhs->t = 1;

    se_vel_get_vel_grad_hess_at_pt(v_app, &ray->xyz, &ray->pqr, 
                                    v, gx_v, gp_v, hv);

    pqr_norm = v3d_norm(&ray->pqr);
    
    if (ode_cfg->use_wf) {

        int k;
        v3d_t lb_1, lb_2;
        double c0[7];
        v3d_t x0[7];
        v3d_t pqr_hat;

        v3d_assign_scaled(&pqr_hat,1/pqr_norm,&ray->pqr);

        v3d_init_zero(&rhs->pqr);

        for (k=1; k<=ode_cfg->wf_sets; k++) {

            double wf_size = (k*ode_cfg->wf_size)/ode_cfg->wf_sets;
            double d_0, d_1;

            get_x_on_wf(ray, v_app->is2d, wf_size, x0, &pqr_hat, &lb_1, &lb_2);
            get_vel_on_wf_pts(x0, &pqr_hat, v_app, ode_cfg->wf_weights[0], ode_cfg->wf_weights[1], c0);

            d_0 = ( c0[2] - c0[1] ) / ( 2 * c0[0] * wf_size);

            if (v_app->is2d) {
                /*                                        _____
                                                          p - bq        p       ( 1- ,/1+b^2 ) p - bq
                                                          dot p = ------------ - ----- =  ----------------------
                                                          dt c |p-bq|   dt c               ______
                                                          dt c ,/1+b^2
                                   
                                                          a dt 
                                                          b = ----------------
                                                          ________________
                                                          ,/ 4w^2 - dt^2 a^2
 
                                                          a = c_{+} - c_{-}

                                                          Limit as dt->0
   
                                                          - a q
                                                          dot p = -----------
                                                          2 c w 
                                
                */

                d_1 = - ( c0[3] - c0[4] ) / (2 * c0[0] * wf_size);
            } else {
                double d_2;

                d_1 = - ( c0[3] - c0[4] + c0[5] - c0[6] ) * M_SQRT2 / (4 * c0[0] * wf_size);

                d_2 = - ( c0[3] - c0[5] + c0[4] - c0[6] ) * M_SQRT2 / (4 * c0[0] * wf_size);

                v3d_accum_scaled(&rhs->pqr, d_2, &lb_2);

            }

            v3d_accum_scaled(&rhs->pqr, d_0, &pqr_hat);
            v3d_accum_scaled(&rhs->pqr, d_1, &lb_1);

        }

        v3d_assign_scaled(&rhs->xyz, *v, &pqr_hat);
        v3d_scale(&rhs->pqr, 1.0/ode_cfg->wf_sets);

    } else {

        v3d_assign_scaled(&rhs->xyz, *v/pqr_norm, &ray->pqr);
        v3d_accum_scaled(&rhs->xyz, pqr_norm, gp_v);
        v3d_assign_scaled(&rhs->pqr, -pqr_norm, gx_v);

    }

 
}
