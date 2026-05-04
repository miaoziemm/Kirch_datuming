#include "../include/se_vel.h"


const int cube_order[8][3]={{0,0,0}, 
                            {0,0,1},
                            {1,0,0},
                            {1,0,1},
                            {0,1,0},
                            {0,1,1},
                            {1,1,0},
                            {1,1,1}};

/**
 * Velocity approximation cell.
 * 
 */
struct se_vel_approx_cell_s {

    float *vel_loc; /**< the local velocity data needed to compute the
                     * derivatives allocating enough space for the
                     * largest stencil, so the data has all of the
                     * values needed to compute at the 8 cell
                     * corners */
    float *del_loc; /**< the local delta data needed to compute the
                     * derivatives */
    float *eta_loc; /**< the local eta data needed to compute the
                     * derivatives */
    v3d_t tti_loc[8]; /**< Axis of symmetry at the corners. Uses packing order (z,x,y). */

    size_t n_xyz[3]; /**< dimensions of the data_local array */
    size_t n_elem; /**< total size of data_local array */
    ssize_t bp_xyz[3]; /**< index of the base point location in data_local */
    double o_xyz[3]; /**< origin of the data_local */
    double d_xyz[3]; /**< sampling of the data_local */
    double over_d_xyz[3]; /**< 1/sampling of the data_local */

    int valid_loc; /**< True vel_data has been filled correctly for the base indexes. */

    double vel[8]; /**< Velocity at the corners. Uses packing order (z,x,y). */
    double del[8]; /**< Delta at the corners. Uses packing order (z,x,y). */
    double eta[8]; /**< Eta at the corners. Uses packing order (z,x,y). */
    int valid_v; /**< True if the v has been filled correctly for the base indexes. */

    v3d_t gvel[8]; /**< Gradient of the velocity at the corners. Uses packing order (z,x,y). */
    v3d_t gdel[8]; /**< Gradient of delta at the corners. Uses packing order (z,x,y). */
    v3d_t geta[8]; /**< Gradient of eta at the corners. Uses packing order (z,x,y). */
    int valid_g; /**< True if the gradient has been filled correctly for the base indexes. */

    m3d_t hvel[8]; /**< Hessian of the velocity at the corners. Uses packing order (z,x,y). */
    int valid_h; /**< True if the hessian has been filled correctly for the base indexes. */

    int rc_vel; /**< Flag, true if velocity is changing rapidly */
    int rc_del; /**< Flag, true if delta is changing rapidly */
    int rc_eta; /**< Flag, true if eta is changing rapidly */

    int n_rc_vel; /**< number of sub cells with rapidly changing velocity */
    int n_rc_del; /**< number of sub cells with rapidly changing delta */
    int n_rc_eta; /**< number of sub cells with rapidly changing eta */

    v3d_t dbp_vel[8]; /**< Base point on the velocity discontinuity */
    v3d_t dbp_del[8]; /**< Base point on the delta discontinuity */
    v3d_t dbp_eta[8]; /**< Base point on the eta discontinuity */

    v3d_t dnv_vel[8]; /**< Normal vector to the velocity discontinuity */
    v3d_t dnv_del[8]; /**< Normal vector to the delta discontinuity */
    v3d_t dnv_eta[8]; /**< Normal vector to the eta discontinuity */

    double min_vel[8]; /**< Minimum velocity in local sub cell */
    double min_del[8]; /**< Minimum delta in local sub cell */
    double min_eta[8]; /**< Minimum eta in local sub cell */

    double max_vel[8]; /**< Maximum velocity in local sub cell */
    double max_del[8]; /**< Maximum delta in local sub cell */
    double max_eta[8]; /**< Maximum eta in local sub cell */

    int side_vel[8]; /**< Average velocity in local sub cell */
    int side_del[8]; /**< Average delta in local sub cell */
    int side_eta[8]; /**< Average eta in local sub cell */

};

/**
 * Implementation of the detection of discontinuities. If detected,
 * this funcition computes the discontinuity base point and normal
 * vector to the discontinuity for each of the corners of the
 * approximation cell.
 * 
 */
static void detect_and_fill_cell_rapid_change(const se_vel_approx_cell_t *vc, 
                                                   const se_vel_approx_t *v_app,
                                                   const float *d_loc,
                                                   int *rc_flag, int* n_rc_cells,
                                                   double *min_d, double *max_d, int *side,
                                                   v3d_t *dbp, v3d_t *dnv);
/**
 * Implementation of 3pt gradient computation. It computes the
 * gradient at the corners of the approximation cell. The data cell
 * need not be the size of the stencil, but must contain enougth
 * information to compute the necessary finite differences.
 * 
 * \param[in,out] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 */
static void fill_cell_grad_3pt_gen(const se_vel_approx_cell_t *vc, const float *data, int is2d, v3d_t *g);

/**
 * Implementation of 3pt hessian matrix computation. It computes the 
 * hessian at the corners of the approximation cell. The data cell
 * need not be the size of the stencil, but must contain enougth
 * information to compute the necessary finite differences.
 * 
 * \param[in,out] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 */
static void fill_cell_hess_3pt_gen(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h);

/**
 * Implementation of 3pt gradient computation. It computes the 
 * gradient at the corners of the approximation cell. 
 * WARN: ASSUMES THAT THE CELL IS PRECISELY THE 3PT STENCIL CELL. 
 * 
 * \param[in,out] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 */
static void fill_cell_grad_3pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, v3d_t *g);

/**
 * Implementation of 3pt hessian matrix computation. It computes the 
 * hessian at the corners of the approximation cell.
 * WARN: ASSUMES THAT THE CELL IS PRECISELY THE 3PT STENCIL CELL. 
 * 
 * \param[in,out] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 */
static void fill_cell_hess_3pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h);

/**
 * Implementation of 5pt gradient computation. It computes the 
 * gradient at the corners of the approximation cell. 
 * WARN: ASSUMES THAT THE CELL IS PRECISELY THE 5PT STENCIL CELL. 
 * 
 * \param[in,out] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 */
static void fill_cell_grad_5pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, v3d_t *g);

/**
 * Implementation of 5pt hessian matrix computation. It computes the 
 * hessian at the corners of the approximation cell.
 * WARN: ASSUMES THAT THE CELL IS PRECISELY THE 5PT STENCIL CELL. 
 * 
 * \param[in] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 * \param[out] h array of hessian matrices at the corners of the aprroximation cell.
 */
static void fill_cell_hess_5pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h);

/**
 * Implementation of 7pt gradient computation. It computes the 
 * gradient at the corners of the approximation cell. 
 * WARN: ASSUMES THAT THE CELL IS PRECISELY THE 7PT STENCIL CELL. 
 * 
 * \param[in,out] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 */
static void fill_cell_grad_7pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, v3d_t *g);

/**
 * Implementation of 7pt hessian matrix computation. It computes the 
 * hessian at the corners of the approximation cell.
 * WARN: ASSUMES THAT THE CELL IS PRECISELY THE 7PT STENCIL CELL. 
 * 
 * \param[in] vc the velocity approximation cell
 * \param[in] is2d 2D flag
 * \param[out] h array of hessian matrices at the corners of the aprroximation cell.
 */
static void fill_cell_hess_7pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h);

/**
 * Initialize a velocity approximation cell.
 * 
 * \param[out] vc the velocity approximation cell to initialize
 * \param[in] st_hl the approximation stencil half length
 * \param[in] d_xyz sampling of the velocity
 * \param[in] over_d_xyz 1 / sampling of the velocity
 * 
 */
static void vel_approx_cell_init(se_vel_approx_cell_t *vc, int aniso, int st_hl, 
                                      const double *d_xyz, const double *over_d_xyz);

/**
 * De-allocate the velocity approximation cell.
 * 
 * \param[in] vc the velocity approximation cell to deallocate
 * 
 */
static void destroy_vel_approx_cell(se_vel_approx_cell_t *vc);


void se_vel_init(se_vel_t *vel, float *vel_data, 
                  float *del_data, float *eta_data,
                  float *ttix_data, float *ttiy_data,
                  const grid3d_t *velgrid,
                  int force_3d)
{

    int64_t i;
    
    if (velgrid) {
        vel->n_xyz[0] = velgrid->x.n;
        vel->o_xyz[0] = velgrid->x.o;
        vel->d_xyz[0] = velgrid->x.d;

        vel->n_xyz[1] = velgrid->y.n;
        vel->o_xyz[1] = velgrid->y.o;
        vel->d_xyz[1] = velgrid->y.d;

        vel->n_xyz[2] = velgrid->z.n;
        vel->o_xyz[2] = velgrid->z.o;
        vel->d_xyz[2] = velgrid->z.d;

        if (vel->n_xyz[1]==1) { 
            vel->d_xyz[1]=1;
        }
        for (i=0; i<3; i++) {
            vel->m_xyz[i] = vel->o_xyz[i] + (vel->n_xyz[i]-1)*vel->d_xyz[i];
            vel->over_d_xyz[i] = 1./vel->d_xyz[i];
        }
    } else {
        for (i=0; i<3; i++) {
            vel->n_xyz[i] = 1;
            vel->o_xyz[i] = 0;
            vel->d_xyz[i] = 1;
            vel->m_xyz[i] = 0;
            vel->over_d_xyz[i] = 1;
        }        
    }

    
    vel->is3d = (vel->n_xyz[1]==1) ? 0 : 1;
    if (force_3d) vel->is3d=1;

    vel->vel_data=vel_data;
    vel->del_data=del_data;
    vel->eta_data=eta_data;
    vel->ttix_data=ttix_data;
    vel->ttiy_data=ttiy_data;

    if ( (vel->del_data) && (vel->eta_data) ) {
        vel->aniso=1;
        if (ttix_data || ttiy_data) {
            vel->tti=1;
        } else {
            vel->tti=0;
        }
    } else {
        vel->aniso=0;
        vel->tti=0;
    }

    if (vel->vel_data) {
        for (i=0; i<(int64_t)vel->n_xyz[0]*(int64_t)vel->n_xyz[1]*(int64_t)vel->n_xyz[2]; i++) {
            if ( ( vel->vel_data[i]<0 ) || ( FISZERO(vel->vel_data[i]) ) ) {
                ERROR(("Velocity constains 0 or negative values. Please correct."));
            }
        }
    }

}

size_t se_vel_get_mem_usage(se_vel_t *vel)
{
    size_t mem=0;
    int64_t nxyz;
    
    nxyz = (int64_t)(vel->n_xyz[0]*vel->n_xyz[1]*vel->n_xyz[2]);

    if (vel->vel_data)   mem += sizeof(*(vel->vel_data))*nxyz;
    if (vel->del_data)   mem += sizeof(*(vel->del_data))*nxyz;
    if (vel->eta_data)   mem += sizeof(*(vel->eta_data))*nxyz;
    if (vel->ttix_data)  mem += sizeof(*(vel->ttix_data))*nxyz;
    if (vel->ttiy_data)  mem += sizeof(*(vel->ttiy_data))*nxyz;

    mem += sizeof(*vel);

    return mem;
}

void se_vel_destroy(se_vel_t *vel)
{
    free(vel->vel_data);
    free(vel->del_data);
    free(vel->eta_data);
    free(vel->ttix_data);
    free(vel->ttiy_data);
}

static void vel_approx_cell_init(se_vel_approx_cell_t *vc, int aniso, int st_hl, 
                                      const double *d_xyz, const double *over_d_xyz) 
{
    int i;

    vc->rc_vel=0;
    vc->rc_del=0;
    vc->rc_eta=0;

    vc->valid_loc = 0;

    vc->valid_v = 0;
    vc->valid_g = 0;
    vc->valid_h = 0;

    vc->n_xyz[0] = 2*st_hl+1+1;
    vc->n_xyz[1] = 2*st_hl+1+1;
    vc->n_xyz[2] = 2*st_hl+1+1;
    vc->n_elem = vc->n_xyz[0]*vc->n_xyz[1]*vc->n_xyz[2];

    vc->bp_xyz[0] = st_hl;
    vc->bp_xyz[1] = st_hl;
    vc->bp_xyz[2] = st_hl;

    vc->vel_loc = alloc1float_zero(vc->n_elem);
    if (aniso) {
        vc->del_loc = alloc1float_zero(vc->n_elem);
        vc->eta_loc = alloc1float_zero(vc->n_elem);

    } else {
        vc->del_loc = NULL;
        vc->eta_loc = NULL;
    }

    for (i=0; i<8; i++) {
        v3d_init(vc->tti_loc+i,0,0,1);
    }

    memcpy(vc->over_d_xyz, over_d_xyz, 3*sizeof(double));
    memcpy(vc->d_xyz, d_xyz, 3*sizeof(double));

    
}

size_t se_vel_approx_cell_get_mem_usage(se_vel_approx_cell_t * vc)
{
    size_t mem;

    if (!vc) return 0;

    mem = sizeof(*vc);
    if (vc->vel_loc) mem += vc->n_elem*sizeof(*vc->vel_loc);
    if (vc->del_loc) mem += vc->n_elem*sizeof(*vc->del_loc);
    if (vc->eta_loc) mem += vc->n_elem*sizeof(*vc->eta_loc);

    return mem;
    
}

static void destroy_vel_approx_cell(se_vel_approx_cell_t *vc) 
{
    free(vc->vel_loc);
    free(vc->del_loc);
    free(vc->eta_loc);
}

void se_vel_approx_init(se_vel_approx_t *v_app, 
                         se_vel_t *vel,
                         int stencil, int use_hess,
                         int dtct_bdry, double bdry_frac)
{

    v_app->vel = vel;
    v_app->aniso = vel->aniso;
    v_app->tti = vel->tti;

    v_app->is2d = !(vel->is3d);

    v_app->st_hl = (stencil/2);
    v_app->use_hess = use_hess;
    v_app->detect_bdry = dtct_bdry;
    v_app->bdry_frac = bdry_frac;

    v_app->vc_curr =(se_vel_approx_cell_t*)malloc(sizeof(se_vel_approx_cell_t));
    memset(v_app->vc_curr, 0, sizeof(se_vel_approx_cell_t));
    vel_approx_cell_init(v_app->vc_curr, vel->del_data?1:0, 
                         v_app->st_hl, v_app->vel->d_xyz, v_app->vel->over_d_xyz);
    v_app->vc_curr_lidx = -1;

    v_app->vc_prev = (se_vel_approx_cell_t*)malloc(sizeof(se_vel_approx_cell_t));
    memset(v_app->vc_prev, 0, sizeof(se_vel_approx_cell_t));

    vel_approx_cell_init(v_app->vc_prev, vel->del_data?1:0, 
                         v_app->st_hl, v_app->vel->d_xyz, v_app->vel->over_d_xyz);
    v_app->vc_prev_lidx = -1;

    v_app->re_used = 0;
    v_app->re_used_cache = 0;

    v_app->lidx_n_xyz[0] = v_app->vel->n_xyz[0] + stencil;
    v_app->lidx_n_xyz[1] = v_app->is2d ? 1 : v_app->vel->n_xyz[1] + stencil;
    v_app->lidx_n_xyz[2] = v_app->vel->n_xyz[2] + stencil;

    switch (v_app->st_hl) {
    case 1:
        v_app->eval_grad=&fill_cell_grad_3pt;
        v_app->eval_hess=&fill_cell_hess_3pt;
        break;
    case 2:
        v_app->eval_grad=&fill_cell_grad_5pt;
        v_app->eval_hess=&fill_cell_hess_5pt;
        break;
    case 3:
        v_app->eval_grad=&fill_cell_grad_7pt;
        v_app->eval_hess=&fill_cell_hess_7pt;
        break;
    default:
        v_app->eval_grad=&fill_cell_grad_3pt_gen;
        v_app->eval_hess=&fill_cell_hess_3pt_gen;
        break;
        /* default: */
        /*     ERROR(("%d is an unsupported stencil half length",v_app->st_hl)); */
    }

}

size_t se_vel_approx_get_mem_usage(const se_vel_approx_t *v_app)
{
    size_t mem;
    if (!v_app) return 0;
    
    mem = sizeof(v_app);
    mem += se_vel_approx_cell_get_mem_usage(v_app->vc_curr);
    mem += se_vel_approx_cell_get_mem_usage(v_app->vc_prev);
    
    return mem;

}

void se_vel_approx_destroy(se_vel_approx_t *v_app)
{

    INFOV((5,"Vel caching: re-used %Ld  --  re-used from cache %Ld",
           (long long int)v_app->re_used, (long long int)v_app->re_used_cache));

    destroy_vel_approx_cell(v_app->vc_curr);
    free(v_app->vc_curr);
    destroy_vel_approx_cell(v_app->vc_prev);
    free(v_app->vc_prev);
}

inline void pt_to_vel_idx(const se_vel_approx_t *v_app, const v3d_t *pt, ssize_t *idx)
{

    idx[0] =                   (ssize_t) ((pt->v[0] - v_app->vel->o_xyz[0]) * v_app->vel->over_d_xyz[0]);
    idx[1] = v_app->is2d ? 0 : (ssize_t) ((pt->v[1] - v_app->vel->o_xyz[1]) * v_app->vel->over_d_xyz[1]);
    idx[2] =                   (ssize_t) ((pt->v[2] - v_app->vel->o_xyz[2]) * v_app->vel->over_d_xyz[2]);

}

inline void get_trilinear_weights(const se_vel_approx_t *v_app, ssize_t *idx, 
                                      const v3d_t *pt, double *wts)
{

    wts[0] =                   1 - ( (pt->v[0] - v_app->vel->o_xyz[0])*v_app->vel->over_d_xyz[0] - idx[0]);
    wts[1] = v_app->is2d ? 1 : 1 - ( (pt->v[1] - v_app->vel->o_xyz[1])*v_app->vel->over_d_xyz[1] - idx[1]);
    wts[2] =                   1 - ( (pt->v[2] - v_app->vel->o_xyz[2])*v_app->vel->over_d_xyz[2] - idx[2]);

}

static void trilinear_approx_double(const double *wts, const double *cube_d, 
                                         int is2d, double *d )
{

    (*d)  = (  wts[2]) *  (  wts[0]) * (  wts[1]) * cube_d[0] ;
    (*d) += (1-wts[2]) *  (  wts[0]) * (  wts[1]) * cube_d[1] ;
    (*d) += (  wts[2]) *  (1-wts[0]) * (  wts[1]) * cube_d[2] ;
    (*d) += (1-wts[2]) *  (1-wts[0]) * (  wts[1]) * cube_d[3] ;

    if (is2d) return;

    (*d) += (  wts[2]) *  (  wts[0]) * (1-wts[1]) * cube_d[4] ;
    (*d) += (1-wts[2]) *  (  wts[0]) * (1-wts[1]) * cube_d[5] ;
    (*d) += (  wts[2]) *  (1-wts[0]) * (1-wts[1]) * cube_d[6] ;
    (*d) += (1-wts[2]) *  (1-wts[0]) * (1-wts[1]) * cube_d[7] ;

}

static void trilinear_approx_vec3d(const double *wts, const v3d_t *grad_cube, 
                                        int is2d, v3d_t *gv )
{
    
    v3d_assign_scaled(gv, (  wts[2]) *  (  wts[0]) * (  wts[1]), grad_cube + 0);
    v3d_accum_scaled (gv, (1-wts[2]) *  (  wts[0]) * (  wts[1]), grad_cube + 1);
    v3d_accum_scaled (gv, (  wts[2]) *  (1-wts[0]) * (  wts[1]), grad_cube + 2);
    v3d_accum_scaled (gv, (1-wts[2]) *  (1-wts[0]) * (  wts[1]), grad_cube + 3);

    if (is2d) return;

    v3d_accum_scaled (gv, (  wts[2]) *  (  wts[0]) * (1-wts[1]), grad_cube + 4);
    v3d_accum_scaled (gv, (1-wts[2]) *  (  wts[0]) * (1-wts[1]), grad_cube + 5);
    v3d_accum_scaled (gv, (  wts[2]) *  (1-wts[0]) * (1-wts[1]), grad_cube + 6);
    v3d_accum_scaled (gv, (1-wts[2]) *  (1-wts[0]) * (1-wts[1]), grad_cube + 7);

}

static void trilinear_approx_double_disc(const double *wts, int max_side, int is2d, 
                                              const double *cube_d, const int *side_d,  
                                              const double *min_d, const double *max_d, 
                                              double *d)
{

    if (max_side) {

        (*d)  = (  wts[2]) *  (  wts[0]) * (  wts[1]) * ( (side_d[0]<.5)? max_d[0] : cube_d[0] ) ;
        (*d) += (1-wts[2]) *  (  wts[0]) * (  wts[1]) * ( (side_d[1]<.5)? max_d[1] : cube_d[1] ) ;
        (*d) += (  wts[2]) *  (1-wts[0]) * (  wts[1]) * ( (side_d[2]<.5)? max_d[2] : cube_d[2] ) ;
        (*d) += (1-wts[2]) *  (1-wts[0]) * (  wts[1]) * ( (side_d[3]<.5)? max_d[3] : cube_d[3] ) ;

        if (is2d) return;

        (*d) += (  wts[2]) *  (  wts[0]) * (1-wts[1]) * ( (side_d[4]<.5)? max_d[4] : cube_d[4] ) ;
        (*d) += (1-wts[2]) *  (  wts[0]) * (1-wts[1]) * ( (side_d[5]<.5)? max_d[5] : cube_d[5] ) ;
        (*d) += (  wts[2]) *  (1-wts[0]) * (1-wts[1]) * ( (side_d[6]<.5)? max_d[6] : cube_d[6] ) ;
        (*d) += (1-wts[2]) *  (1-wts[0]) * (1-wts[1]) * ( (side_d[7]<.5)? max_d[7] : cube_d[7] ) ;

    } else {

        (*d)  = (  wts[2]) *  (  wts[0]) * (  wts[1]) * ( (side_d[0]>-.5)? min_d[0] : cube_d[0] ) ;
        (*d) += (1-wts[2]) *  (  wts[0]) * (  wts[1]) * ( (side_d[1]>-.5)? min_d[1] : cube_d[1] ) ;
        (*d) += (  wts[2]) *  (1-wts[0]) * (  wts[1]) * ( (side_d[2]>-.5)? min_d[2] : cube_d[2] ) ;
        (*d) += (1-wts[2]) *  (1-wts[0]) * (  wts[1]) * ( (side_d[3]>-.5)? min_d[3] : cube_d[3] ) ;

        if (is2d) return;

        (*d) += (  wts[2]) *  (  wts[0]) * (1-wts[1]) * ( (side_d[4]>-.5)? min_d[4] : cube_d[4] ) ;
        (*d) += (1-wts[2]) *  (  wts[0]) * (1-wts[1]) * ( (side_d[5]>-.5)? min_d[5] : cube_d[5] ) ;
        (*d) += (  wts[2]) *  (1-wts[0]) * (1-wts[1]) * ( (side_d[6]>-.5)? min_d[6] : cube_d[6] ) ;
        (*d) += (1-wts[2]) *  (1-wts[0]) * (1-wts[1]) * ( (side_d[7]>-.5)? min_d[7] : cube_d[7] ) ;
    }
}

static void trilinear_approx_double_rc(const double *tri_wts,
                                            int is2d,
                                            const double *cube_d, const int *side_d,
                                            const double *min_d, const double *max_d, 
                                            const v3d_t *dbp, const v3d_t *dnv,
                                            double *d)
{

    double side=0;
    v3d_t bp, nv;

    trilinear_approx_vec3d(tri_wts, dbp, is2d, &bp);
    trilinear_approx_vec3d(tri_wts, dnv, is2d, &nv);

    side = ( ( (1-tri_wts[0]) - bp.v[0] ) * nv.v[0] +
             ( (1-tri_wts[1]) - bp.v[1] ) * nv.v[1] +
             ( (1-tri_wts[2]) - bp.v[2] ) * nv.v[2] );
        
    trilinear_approx_double_disc(tri_wts, (side>0), is2d, cube_d, side_d, min_d, max_d, d);

}

static void trilinear_approx_vec3d_rc(const double *tri_wts,
                                           int is2d, const double *d_xyz,
                                           const v3d_t *grad_cube, 
                                           const double *min_d, const double *max_d,
                                           const v3d_t *dbp, const v3d_t *dnv,
                                           v3d_t *gv)
{
    double dist;
    v3d_t bp, nv;
    double min, max;
    double len;
    double range;

    trilinear_approx_vec3d(tri_wts, dbp, is2d, &bp);
    trilinear_approx_vec3d(tri_wts, dnv, is2d, &nv);
    trilinear_approx_double(tri_wts, min_d, is2d, &min);
    trilinear_approx_double(tri_wts, max_d, is2d, &max);

    trilinear_approx_vec3d(tri_wts, grad_cube, is2d, gv);

    len = sqrt((            nv.v[0]*d_xyz[0] * nv.v[0]*d_xyz[0] )+
               ( is2d ? 0 : nv.v[1]*d_xyz[1] * nv.v[1]*d_xyz[1] )+
               (            nv.v[2]*d_xyz[2] * nv.v[2]*d_xyz[2] ) );

    dist = se_fabs( ( (1-tri_wts[0]) - bp.v[0] ) * nv.v[0]*d_xyz[0] +
                     ( (1-tri_wts[1]) - bp.v[1] ) * nv.v[1]*d_xyz[1] +
                     ( (1-tri_wts[2]) - bp.v[2] ) * nv.v[2]*d_xyz[2] );
        
    range=1*len;
    if (dist < range) {
        v3d_init(gv, 
                     (            nv.v[0]/d_xyz[0] ),
                     ( is2d ? 0 : nv.v[1]/d_xyz[1] ),
                     (            nv.v[2]/d_xyz[2] ) );
        v3d_set_norm(gv, ( max-min ) * ( (1-dist/range) / range) );

    } else {
        v3d_init(gv, 0, 0, 0);
    }
}

static void trilinear_approx_hess(const double *wts, const m3d_t *hess_cube, 
                                       int is2d, m3d_t *hv )
{
    
    m3d_assign_scaled(hv, (  wts[2]) *  (  wts[0]) * (  wts[1]), hess_cube + 0);
    m3d_accum_scaled (hv, (1-wts[2]) *  (  wts[0]) * (  wts[1]), hess_cube + 1);
    m3d_accum_scaled (hv, (  wts[2]) *  (1-wts[0]) * (  wts[1]), hess_cube + 2);
    m3d_accum_scaled (hv, (1-wts[2]) *  (1-wts[0]) * (  wts[1]), hess_cube + 3);

    if (is2d) return;

    m3d_accum_scaled (hv, (  wts[2]) *  (  wts[0]) * (1-wts[1]), hess_cube + 4);
    m3d_accum_scaled (hv, (1-wts[2]) *  (  wts[0]) * (1-wts[1]), hess_cube + 5);
    m3d_accum_scaled (hv, (  wts[2]) *  (1-wts[0]) * (1-wts[1]), hess_cube + 6);
    m3d_accum_scaled (hv, (1-wts[2]) *  (1-wts[0]) * (1-wts[1]), hess_cube + 7);

}


static void invalidate_approx_cell(se_vel_approx_cell_t *vc) 
{
    vc->valid_loc = 0;
    vc->valid_v = 0;
    vc->valid_g = 0;
    vc->valid_h = 0;
}

static int pt_idx_is_at_bdry(const se_vel_approx_t *v_app, const ssize_t *pt_idx) 
{

    if (pt_idx[0]<v_app->st_hl) return 1;
    if (pt_idx[2]<v_app->st_hl) return 1;
    if (pt_idx[0]>=((ssize_t)v_app->vel->n_xyz[0]-v_app->st_hl-1) ) return 1;
    if (pt_idx[2]>=((ssize_t)v_app->vel->n_xyz[2]-v_app->st_hl-1) ) return 1;

    if (v_app->is2d) return 0;

    if (pt_idx[1]<v_app->st_hl) return 1;
    if (pt_idx[1]>=((ssize_t)v_app->vel->n_xyz[1]-v_app->st_hl-1) ) return 1;

    return 0;

}

static void fill_cell_local_data(se_vel_approx_cell_t *vc, 
                                      const se_vel_approx_t *v_app, 
                                      const ssize_t *idx)
{

    const size_t *vc_n_xyz = vc->n_xyz;
    const size_t *vel_n_xyz = v_app->vel->n_xyz;
    const float *const vel_data = v_app->vel->vel_data;
    const float *const del_data = v_app->vel->del_data;
    const float *const eta_data = v_app->vel->eta_data;

    float *const vel_loc = vc->vel_loc;
    float *const del_loc = vc->del_loc;
    float *const eta_loc = vc->eta_loc;

    const int is2d = v_app->is2d;
    const int st_hl = v_app->st_hl;

    vc->valid_loc = 1;

    if (!pt_idx_is_at_bdry(v_app, idx)) {
     
        ssize_t l_j;
        for (l_j=0; l_j<(ssize_t)vc_n_xyz[1]; l_j++) {
            ssize_t f_j, l_i;
            ssize_t f_y, l_y;

            if (is2d) {
                if (l_j!=st_hl) continue;
                f_j = 0;
            } else {
                f_j = idx[1] + l_j - st_hl;
            }

            l_y =  vc_n_xyz[0] * vc_n_xyz[2]  * l_j;
            f_y = vel_n_xyz[0] * vel_n_xyz[2] * f_j;
            

            for (l_i=0; l_i<(ssize_t)vc_n_xyz[0]; l_i++) {
                ssize_t f_i;
                ssize_t f_x, l_x;

                f_i = idx[0] + l_i - st_hl;

                l_x = l_y +  vc_n_xyz[2] * l_i;
                f_x = f_y + vel_n_xyz[2] * f_i + idx[2] - st_hl;

                memcpy(vel_loc + l_x, vel_data + f_x, sizeof(float)*(vc_n_xyz[2]));
                if (v_app->aniso) {
                    memcpy(del_loc + l_x, del_data + f_x, sizeof(float)*(vc_n_xyz[2]));
                    memcpy(eta_loc + l_x, eta_data + f_x, sizeof(float)*(vc_n_xyz[2]));
                }

            }
        }

    } else {
    
        /* If at the boundary, just extend the velocity by using the neareset value */

        ssize_t l_j;
        for (l_j=0; l_j <(ssize_t)vc_n_xyz[1]; l_j++) {
            ssize_t f_j, l_i;
            ssize_t f_y, l_y;

            if (is2d) {
                if ((l_j!=st_hl)) continue;
                f_j = 0; 
            } else {
                f_j = idx[1] + l_j - st_hl;
                if ( f_j >= (ssize_t)vel_n_xyz[1] ) f_j = vel_n_xyz[1]-1;
                if ( f_j <                      0 ) f_j = 0;
            }

            l_y =  vc_n_xyz[0] * vc_n_xyz[2]  * l_j;
            f_y = vel_n_xyz[0] * vel_n_xyz[2] * f_j;


            for (l_i=0; l_i<(ssize_t)vc_n_xyz[0]; l_i++) {
                ssize_t f_i, l_k;
                ssize_t f_x, l_x;

                f_i = idx[0] + l_i - st_hl;
                if (f_i >= (ssize_t)vel_n_xyz[0] ) f_i = vel_n_xyz[0]-1;
                if (f_i < 0 ) f_i = 0;

                l_x = l_y + vc_n_xyz[2]*l_i;
                f_x = f_y + vel_n_xyz[2]*f_i;


                for (l_k=0; l_k<(ssize_t)vc_n_xyz[2]; l_k++) {
                    ssize_t f_k;
                    f_k = idx[2] + l_k - st_hl;
           
                    if (f_k >= (ssize_t)vel_n_xyz[2] ) f_k = vel_n_xyz[2]-1;
                    if (f_k < 0 ) f_k = 0;
                    
                    *(vel_loc + l_x + l_k) = *(vel_data + f_x + f_k);
                    if (v_app->aniso) {
                        *(del_loc + l_x + l_k) = *(del_data + f_x + f_k);
                        *(eta_loc + l_x + l_k) = *(eta_data + f_x + f_k);
                    }
                    
                }
            }
        }
    }

    if ((v_app->aniso) && v_app->tti) {
        int i;
        for (i=0; i<8; i++) {
            int64_t f_i,f_j,f_k,j;
            double dxdz,dydz;

            if ( (is2d) && (i>=4) ) break;
            
            f_i = idx[0] + cube_order[i][0];
            f_j = idx[1] + cube_order[i][1];
            f_k = idx[2] + cube_order[i][2];
            if (f_i>=(int64_t)vel_n_xyz[0]) f_i=vel_n_xyz[0]-1;
            if (f_i<0)             f_i=0;
            if (f_j>=(int64_t)vel_n_xyz[1]) f_j=vel_n_xyz[1]-1;
            if (f_j<0)             f_j=0;
            if (f_k>=(int64_t)vel_n_xyz[2]) f_k=vel_n_xyz[2]-1;
            if (f_k<0)             f_k=0;

            j = (vel_n_xyz[0]*f_j + f_i)*vel_n_xyz[2] + f_k;

            dxdz = (v_app->vel->ttix_data)? v_app->vel->ttix_data[j] : 0;
            dydz = (v_app->vel->ttiy_data)? v_app->vel->ttiy_data[j] : 0;
            v3d_init(vc->tti_loc+i,dxdz,dydz,1);
            v3d_normalize(vc->tti_loc+i);            
        }            
    }

    /* Assuming that the local velocity data array is initialized to 0 before this routine
       otherwise (need the 0 value in 2D), have to 
       memset(vc->data_local, 0, sizeof(float)*vc->n_elem); 
       at the start */
      
    if (v_app->detect_bdry) {

        detect_and_fill_cell_rapid_change(vc, v_app,
                                          vc->vel_loc, 
                                          &vc->rc_vel, &vc->n_rc_vel,
                                          vc->min_vel, vc->max_vel, vc->side_vel,
                                          vc->dbp_vel, vc->dnv_vel);

        if (v_app->aniso) {

            detect_and_fill_cell_rapid_change(vc, v_app,
                                              vc->del_loc, 
                                              &vc->rc_del, &vc->n_rc_del,
                                              vc->min_del, vc->max_del, vc->side_del,
                                              vc->dbp_del, vc->dnv_del);

            detect_and_fill_cell_rapid_change(vc, v_app,
                                              vc->eta_loc, 
                                              &vc->rc_eta, &vc->n_rc_eta,
                                              vc->min_eta, vc->max_eta, vc->side_eta,
                                              vc->dbp_eta, vc->dnv_eta);

        }
        
    }

} 


static void fill_cell_v(se_vel_approx_cell_t *vc,
                             const se_vel_approx_t *v_app,
                             const ssize_t *idx,
                             const size_t *dim)
{

    const float *const vel_loc = vc->vel_loc;
    const float *const del_loc = vc->del_loc;
    const float *const eta_loc = vc->eta_loc;

    double *const vc_vel = vc->vel;
    double *const vc_del = vc->del;
    double *const vc_eta = vc->eta;

    size_t off;

    vc->valid_v = 1;

    off = dim[2]*(dim[0]*idx[1] + idx[0]) + idx[2];

    vc_vel[0] = vel_loc[off              ];
    vc_vel[1] = vel_loc[off          + 1 ];
    vc_vel[2] = vel_loc[off + dim[2]     ];
    vc_vel[3] = vel_loc[off + dim[2] + 1 ];

    if (v_app->aniso) {
        vc_del[0] = del_loc[off              ];
        vc_del[1] = del_loc[off          + 1 ];
        vc_del[2] = del_loc[off + dim[2]     ];
        vc_del[3] = del_loc[off + dim[2] + 1 ];

        vc_eta[0] = eta_loc[off              ];
        vc_eta[1] = eta_loc[off          + 1 ];
        vc_eta[2] = eta_loc[off + dim[2]     ];
        vc_eta[3] = eta_loc[off + dim[2] + 1 ];
    }

    if (v_app->is2d) return;

    off += dim[0]*dim[2];

    vc_vel[4] = vel_loc[off              ];
    vc_vel[5] = vel_loc[off          + 1 ];
    vc_vel[6] = vel_loc[off + dim[2]     ];
    vc_vel[7] = vel_loc[off + dim[2] + 1 ];

    if (v_app->aniso) {
        vc_del[4] = del_loc[off              ];
        vc_del[5] = del_loc[off          + 1 ];
        vc_del[6] = del_loc[off + dim[2]     ];
        vc_del[7] = del_loc[off + dim[2] + 1 ];

        vc_eta[4] = eta_loc[off              ];
        vc_eta[5] = eta_loc[off          + 1 ];
        vc_eta[6] = eta_loc[off + dim[2]     ];
        vc_eta[7] = eta_loc[off + dim[2] + 1 ];
    }

}

double se_vel_get_vel_change_frac_from_grad(se_vel_approx_t *v_app, const v3d_t *pt, const v3d_t *dir)
{
    double v;
    v3d_t gv;

    se_vel_get_vel_at_pt(v_app, pt, dir, &v); 
    se_vel_get_grad_x_at_pt(v_app, pt, dir, &gv); 

    return se_fabs(gv.v[0]*v_app->vel->d_xyz[0] + 
                    gv.v[1]*v_app->vel->d_xyz[1] + 
                    gv.v[2]*v_app->vel->d_xyz[2])/v;

} 

inline int64_t idx_to_lidx(const se_vel_approx_t *v_app, const ssize_t *idx)
{
    
    int64_t key;

    key  = v_app->is2d? 0 : mini64( v_app->lidx_n_xyz[1] - 1, maxi64( idx[1] + v_app->st_hl + 1, 0 ) )*v_app->lidx_n_xyz[0];
    key += mini64( v_app->lidx_n_xyz[0] - 1, maxi64( idx[0] + v_app->st_hl + 1, 0 ) );
    key *= v_app->lidx_n_xyz[2];
    key += mini64( v_app->lidx_n_xyz[2] - 1, maxi64( idx[2] + v_app->st_hl + 1, 0 ) );

    return key;
}

static se_vel_approx_cell_t* get_vc_to_use(se_vel_approx_t *v_app, const ssize_t *pt_idx)
{

    se_vel_approx_cell_t *vc;
    int64_t pt_lidx;

    pt_lidx = idx_to_lidx(v_app, pt_idx);

    if (pt_lidx == v_app->vc_curr_lidx) {
        vc = v_app->vc_curr;
        v_app->re_used++;

    } else {

        if (pt_lidx != v_app->vc_prev_lidx) {
            invalidate_approx_cell(v_app->vc_prev);
        } else {
            v_app->re_used_cache++;
        }

        vc = v_app->vc_prev;
        v_app->vc_prev      = v_app->vc_curr;
        v_app->vc_prev_lidx = v_app->vc_curr_lidx;
        v_app->vc_curr      = vc;
        v_app->vc_curr_lidx = pt_lidx;

    }

    return vc;

}

int se_vel_has_rapid_change(se_vel_approx_t *v_app, const v3d_t *pt)
{
    ssize_t pt_idx[3];
    se_vel_approx_cell_t *vc;

    pt_to_vel_idx(v_app, pt, pt_idx);

    vc = get_vc_to_use(v_app, pt_idx);

    if (!(vc->valid_loc)) {
        fill_cell_local_data(vc, v_app, pt_idx);
    }

    if (v_app->detect_bdry) {
        if ( (vc->rc_vel) ||
             ( ( v_app->aniso) && ( (vc->rc_del) || ( vc->rc_eta) ) ) ) {
            return 1;
        }
    }

    return 0;

}


int se_vel_dist_to_rapid_change(se_vel_approx_t *v_app, const v3d_t *pt, double *dist, v3d_t *nv)
{
    ssize_t pt_idx[3];
    se_vel_approx_cell_t *vc;
    
    pt_to_vel_idx(v_app, pt, pt_idx);

    vc = get_vc_to_use(v_app, pt_idx);

    if (!(vc->valid_loc)) {
        fill_cell_local_data(vc, v_app, pt_idx);
    }

    if (v_app->detect_bdry) {
        if ( (vc->rc_vel) ||
             ( ( v_app->aniso) && ( (vc->rc_del) || ( vc->rc_eta) ) ) ) {

            v3d_t bp;
            double tri_wts[3];

            get_trilinear_weights(v_app, pt_idx, pt, tri_wts);

            trilinear_approx_vec3d(tri_wts, vc->dbp_vel, v_app->is2d, &bp);
            trilinear_approx_vec3d(tri_wts, vc->dnv_vel, v_app->is2d, nv);
            nv->v[0] *= vc->d_xyz[0];
            nv->v[1] *= vc->d_xyz[1];
            nv->v[2] *= vc->d_xyz[2];

            (*dist) = ( ( (1-tri_wts[0]) - bp.v[0] ) * nv->v[0] +
                        ( (1-tri_wts[1]) - bp.v[1] ) * nv->v[1] +
                        ( (1-tri_wts[2]) - bp.v[2] ) * nv->v[2]  );

            v3d_set_norm(nv,1);

            return 1;
        }
    }

    (*dist)=0;
    v3d_init_zero(nv);
    return 0;

}

int se_vel_get_vel_grad_hess_at_pt(se_vel_approx_t *v_app, const v3d_t *pt, const v3d_t *dir,
                                    double *v, v3d_t *gx_v, v3d_t *gp_v, m3d_t *hv)
{

    ssize_t pt_idx[3];
    double tri_wts[3];
    double v_iso;
    se_vel_approx_cell_t *vc;



    if ( (!v) && (!gx_v) && (!gp_v) && (!hv) ) return se_vel_pt_is_inside(v_app, pt);
    pt_to_vel_idx(v_app, pt, pt_idx);

    vc = get_vc_to_use(v_app, pt_idx);

    get_trilinear_weights(v_app, pt_idx, pt, tri_wts);


    if (v || v_app->aniso) {
       
        if ( !(vc->valid_v) ) {
        
            if ( !(vc->valid_loc) ) 
            {
            
                 fill_cell_local_data(vc, v_app, pt_idx); 
               
            }
                      
            fill_cell_v(vc, v_app, vc->bp_xyz, vc->n_xyz);
           
        }
        if (vc->rc_vel) {
            trilinear_approx_double_rc(tri_wts, v_app->is2d, vc->vel, 
                                       vc->side_vel, vc->min_vel, vc->max_vel,  
                                       vc->dbp_vel, vc->dnv_vel, &v_iso);
        } else {
           
            trilinear_approx_double(tri_wts, vc->vel, v_app->is2d, &v_iso);
           
        }
    }

    if (gx_v) {
        if ( !(vc->valid_g) ) {
            if ( !(vc->valid_loc) ) 
                fill_cell_local_data(vc, v_app, pt_idx);
            (*v_app->eval_grad)(vc, vc->vel_loc, v_app->is2d, vc->gvel);
            if (v_app->aniso) {
                (*v_app->eval_grad)(vc, vc->del_loc, v_app->is2d, vc->gdel);
                (*v_app->eval_grad)(vc, vc->eta_loc, v_app->is2d, vc->geta);
            }
        }
        if (vc->rc_vel) {
            trilinear_approx_vec3d_rc(tri_wts, v_app->is2d, v_app->vel->d_xyz, 
                                      vc->gvel, vc->min_vel, vc->max_vel,
                                      vc->dbp_vel, vc->dnv_vel, gx_v);
        } else {
            trilinear_approx_vec3d(tri_wts, vc->gvel, v_app->is2d, gx_v);
        }
    }

    if (hv) {
        if (v_app->use_hess) {
        
            if ( !(vc->valid_h) ) {
                if ( !(vc->valid_loc) ) 
                    fill_cell_local_data(vc, v_app, pt_idx);        
            }
            if (vc->rc_vel) {
                m3d_init_zero(hv);
            } else {
                (*v_app->eval_hess)(vc, vc->vel_loc, v_app->is2d, vc->hvel);
                trilinear_approx_hess(tri_wts, vc->hvel, v_app->is2d, hv);
            }
        } else {
            m3d_init_zero(hv);
        }
    }

    if ( v_app->aniso && ( v || gx_v || gp_v ) ) {
    
        v3d_t asa, wpd;
        double del, eta, sin_phi_sq, sqrt_fact;
        double asa_dot_wpd;

        trilinear_approx_vec3d(tri_wts, vc->tti_loc, v_app->is2d, &asa);

        wpd=*dir;
        v3d_normalize(&wpd);
        asa_dot_wpd = v3d_dot(&asa,&wpd);
        sin_phi_sq = asa_dot_wpd;  
        sin_phi_sq = 1-sin_phi_sq*sin_phi_sq;

        if (vc->rc_del) {
            trilinear_approx_double_rc(tri_wts, v_app->is2d, vc->del, 
                                       vc->side_del, vc->min_del, vc->max_del, 
                                       vc->dbp_del, vc->dnv_del, &del);

        } else {
            trilinear_approx_double(tri_wts, vc->del, v_app->is2d, &del);
        }

        if (vc->rc_eta) {
            trilinear_approx_double_rc(tri_wts, v_app->is2d, vc->eta, 
                                       vc->side_eta, vc->min_eta, vc->max_eta, 
                                       vc->dbp_eta, vc->dnv_eta, &eta);
        } else {
            trilinear_approx_double(tri_wts, vc->eta, v_app->is2d, &eta);
        }

        sqrt_fact = sqrt(1 + 2*del*sin_phi_sq + 2*eta*sin_phi_sq*sin_phi_sq);


        if (gx_v) {
            v3d_t gdel, geta;
            if (vc->rc_del) {
                trilinear_approx_vec3d_rc(tri_wts, v_app->is2d, v_app->vel->d_xyz,
                                          vc->gdel, vc->min_del, vc->max_del,
                                          vc->dbp_del, vc->dnv_del, &gdel);
            } else {
                trilinear_approx_vec3d(tri_wts, vc->gdel, v_app->is2d, &gdel);
            }

            if (vc->rc_eta) {
                trilinear_approx_vec3d_rc(tri_wts, v_app->is2d, v_app->vel->d_xyz,
                                          vc->geta, vc->min_eta, vc->max_eta,
                                          vc->dbp_eta, vc->dnv_eta, &geta);
            } else {
                trilinear_approx_vec3d(tri_wts, vc->geta, v_app->is2d, &geta);
            }
            
            v3d_scale(gx_v, sqrt_fact);
            v3d_accum_scaled(gx_v,            sin_phi_sq*v_iso/sqrt_fact, &gdel);
            v3d_accum_scaled(gx_v, sin_phi_sq*sin_phi_sq*v_iso/sqrt_fact, &geta);        
        }

        if (gp_v) {
            v3d_assign_scaled(gp_v, asa_dot_wpd, &wpd);
            v3d_subtract(gp_v, gp_v, &asa);
            v3d_scale(gp_v, ( v_iso*v_iso *
                                  (2*del+4*eta*sin_phi_sq) *
                                  asa_dot_wpd ) );
        }

        if (v) {
            (*v)=sqrt_fact*v_iso;
        }

    } else {
        if (v) (*v)=v_iso;
        if (gp_v) v3d_init(gp_v,0,0,0);
    }

    return se_vel_pt_is_inside(v_app, pt);

}

int se_vel_pt_is_inside_approx_domain(const se_vel_approx_t *v_app, const v3d_t *pt) 
{
    ssize_t pt_idx[3];
    pt_to_vel_idx(v_app, pt, pt_idx);

    return !pt_idx_is_at_bdry(v_app, pt_idx); 
 
}

static void detect_and_fill_cell_rapid_change(const se_vel_approx_cell_t *vc,
                                                   const se_vel_approx_t *v_app, 
                                                   const float *d_loc,
                                                   int *rc_flag, int* n_rc_cells,
                                                   double *min_d, double *max_d, int *side,
                                                   v3d_t *dbp, v3d_t *dnv)
{

    int i;
    double mag, ave, min, max;
    int l_j;

    const size_t * const vc_n_xyz = vc->n_xyz;
    const ssize_t * const vc_bp_xyz = vc->bp_xyz;

    const int is2d = v_app->is2d;
    const int st_hl = v_app->st_hl;
    const double frac = v_app->bdry_frac;

    (*rc_flag)=0;
    (*n_rc_cells)=0;

    min=1e10;
    max=-1e10;

    for (l_j=0; l_j<(int)vc_n_xyz[1]; l_j++) {
        int l_i;
        if ( (is2d) && (l_j!=st_hl) ) continue;
        for (l_i=0; l_i<(int)vc_n_xyz[0]; l_i++) {
            int l_k;
            for (l_k=0; l_k<(int)vc_n_xyz[2]; l_k++) {
                int j;
                j = vc_n_xyz[2]*(l_i + l_j*vc_n_xyz[0]) + l_k;
                if (min > d_loc[j]) { 
                    min = d_loc[j];
                }
                if (max < d_loc[j]) { 
                    max = d_loc[j];
                }
            }
        }
    }

    ave = (min+max)/2;
    mag  = .5*(se_fabs(max)+se_fabs(min));

    if ( ( !FISZERO(mag) ) && ( (se_fabs(max-min)/mag) > frac) ) {

        v3d_t bp, nv;

        v3d_init_zero(&bp);
        v3d_init_zero(&nv);

        for (i=0; i<8; i++) {
        
            double w_bg, w_sm;
            v3d_t pt_bg, pt_sm;

            if ( (is2d) && (i>=4) ) break;
                
            w_bg = 0;
            w_sm = 0;

            v3d_init_zero(&pt_sm);
            v3d_init_zero(&pt_bg);

            min_d[i]=0;
            max_d[i]=0;

            for (l_j=cube_order[i][1]; l_j<(int)(vc_n_xyz[1]-1+cube_order[i][1]); l_j++) {
                int l_i;
                if ( (is2d) && (l_j!=st_hl) ) continue;

                for (l_i=cube_order[i][0]; l_i<(int)(vc_n_xyz[0]-1+cube_order[i][0]); l_i++) {
                    int l_k;
                    for (l_k=cube_order[i][2]; l_k<(int)(vc_n_xyz[2]-1+cube_order[i][2]); l_k++) {
                        int j;
                        j = vc_n_xyz[2]*(l_i + l_j*vc_n_xyz[0]) + l_k;

                        if (d_loc[j]>ave) {
                            v3d_accum_doubles(&pt_bg,
                                                  (double)(l_i),
                                                  (double)(l_j),
                                                  (double)(l_k) );
                            max_d[i] += d_loc[j];
                            w_bg += 1;
                        }

                        if (d_loc[j]<ave) {
                            v3d_accum_doubles(&pt_sm,
                                                  (double)(l_i),
                                                  (double)(l_j),
                                                  (double)(l_k) );

                            min_d[i] += d_loc[j];
                            w_sm += 1;
                        }

                    }
                }
            }

            if (FISZERO(max_d[i])) max_d[i]=ave;
            if (FISZERO(min_d[i])) min_d[i]=ave;

            if ( (!FISZERO(w_bg)) && (!FISZERO(w_sm)) ) {

                (*rc_flag)=1;

                v3d_accum_doubles(&pt_sm,
                                      (double)(-w_sm*vc_bp_xyz[0]),
                                      (double)(-w_sm*vc_bp_xyz[1]),
                                      (double)(-w_sm*vc_bp_xyz[2]) );

                v3d_accum_doubles(&pt_bg,
                                      (double)(-w_bg*vc_bp_xyz[0]),
                                      (double)(-w_bg*vc_bp_xyz[1]),
                                      (double)(-w_bg*vc_bp_xyz[2]) );

                v3d_scale(&pt_sm, 1/w_sm);
                min_d[i] /= w_sm;

                v3d_scale(&pt_bg, 1/w_bg);
                max_d[i] /= w_bg;

                v3d_subtract(dnv+i, &pt_bg, &pt_sm);
                v3d_normalize(dnv+i);

                v3d_add_scaled(dbp+i,   w_bg/(float)(w_bg+w_sm), &pt_sm, w_sm/(float)(w_bg+w_sm), &pt_bg);

                v3d_accum(&bp, dbp+i);
                v3d_accum(&nv, dnv+i);

                (*n_rc_cells) += 1;

            } else {
                max_d[i]=ave;
                min_d[i]=ave;
                v3d_init_zero(dbp+i);
                v3d_init_zero(dnv+i);
            }

        }

        if (is2d) v3d_set_elem(&bp, 0.0, 2);

        if ((*n_rc_cells)>0 ) {
            v3d_scale(&nv, 1.0/(*n_rc_cells));
            v3d_scale(&bp, 1.0/(*n_rc_cells));

            memset(side,0,8*sizeof(int));

            for (i=0; i<8; i++)  {
                int64_t j;

                if ( (is2d) && (i>=4) ) break;

                if (FISZERO(v3d_norm(dnv+i))) {
                    dnv[i] = nv;
                    dbp[i] = bp;
                }
                
                j = vc_n_xyz[2]*( (cube_order[i][0]+vc_bp_xyz[0]) + (cube_order[i][1]+vc_bp_xyz[1])*vc_n_xyz[0] ) 
                    + cube_order[i][2]+vc_bp_xyz[2];

                if (d_loc[j] > ave ) {
                    side[i] =  1;
                }
                if (d_loc[j] < ave ) {
                    side[i] = -1;
                }
            }
        }
    }

}

static void fill_cell_grad_3pt_gen(const se_vel_approx_cell_t *vc, const float* data, int is2d, v3d_t *g) {

    int i;

    const size_t * const vc_n_xyz = vc->n_xyz;

    for (i=0; i<8; i++) {
        const int x_i=vc->bp_xyz[0]+cube_order[i][0];
        const int y_i=vc->bp_xyz[1]+cube_order[i][1];
        const int z_i=vc->bp_xyz[2]+cube_order[i][2];

        if ( (is2d) && (i>=4)) return;

        v3d_init(g+i,
                     .5*vc->over_d_xyz[0] *
                     (  - data[( (x_i-1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i]
                        + data[( (x_i+1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] ),
                     is2d ? 0 :
                     .5*vc->over_d_xyz[1] *
                     (  - data[( x_i + (y_i-1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i]
                        + data[( x_i + (y_i+1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] ),
                     .5*vc->over_d_xyz[2] *
                     (  - data[( x_i + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i-1)]
                        + data[( x_i + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i+1)] ) );

    }

}


static void fill_cell_hess_3pt_gen(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h) {

    int i;

    const size_t * const vc_n_xyz = vc->n_xyz;

    for (i=0; i<8; i++) {
        const int x_i=vc->bp_xyz[0]+cube_order[i][0];
        const int y_i=vc->bp_xyz[1]+cube_order[i][1];
        const int z_i=vc->bp_xyz[2]+cube_order[i][2];

        if ( (is2d) && (i>=4)) return;

        m3d_sym_init(h+i, 
                         vc->over_d_xyz[0] * vc->over_d_xyz[0] * /* xx */
                         (  +  data[((x_i-1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] 
                            -2*data[((x_i  ) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] 
                            +  data[((x_i+1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] ), 
                         is2d ? 0 :
                         vc->over_d_xyz[0] * vc->over_d_xyz[1] * /* xy */
                         (  + data[( (x_i-1) + (y_i-1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] 
                            - data[( (x_i+1) + (y_i-1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] 
                            - data[( (x_i-1) + (y_i+1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i]
                            + data[( (x_i+1) + (y_i+1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] ), 
                         vc->over_d_xyz[0] * vc->over_d_xyz[2] * /* xz */
                         (  + data[( (x_i-1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i-1)] 
                            - data[( (x_i-1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i+1)] 
                            - data[( (x_i+1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i-1)] 
                            + data[( (x_i+1) + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i+1)] ), 
                         is2d ? 0 :
                         vc->over_d_xyz[1] * vc->over_d_xyz[1] * /* yy */
                         (  +  data[( x_i + (y_i-1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] 
                            -2*data[( x_i + (y_i  ) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] 
                            +  data[( x_i + (y_i+1) * vc_n_xyz[0] ) * vc_n_xyz[2] + z_i] ), 
                         is2d ? 0 :
                         vc->over_d_xyz[1] * vc->over_d_xyz[2] * /* yz */
                         (  + data[( x_i + (y_i-1) * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i-1)] 
                            - data[( x_i + (y_i-1) * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i+1)] 
                            - data[( x_i + (y_i+1) * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i-1)] 
                            + data[( x_i + (y_i+1) * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i+1)] ), 
                         vc->over_d_xyz[2] * vc->over_d_xyz[2] * /* zz */
                         (  +  data[( x_i + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i-1)] 
                            -2*data[( x_i + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i  )] 
                            +  data[( x_i + y_i * vc_n_xyz[0] ) * vc_n_xyz[2] + (z_i+1)] ) ); 
    }

}

static void fill_cell_grad_3pt(const se_vel_approx_cell_t *vc, const float* data, int is2d, v3d_t *g) {

    int i;

    for (i=0; i<8; i++) {

        /* WARN: ASSUMES THAT THE CELL IS PRECISELY THE 3PT STENCIL CELL. 
           can be made generic by replacing: 
           --> 1 with vc->bp_xyz[0, 1 or 2] as appropriate 
           --> 4 with vc->n_xyz[0, 1 or 2] as appropriate 
        */

        const int x_i = 1 + cube_order[i][0];
        const int y_i = 1 + cube_order[i][1];
        const int z_i = 1 + cube_order[i][2];

        if ( (is2d) && (i>=4)) return;

        v3d_init(g+i,
                     .5*vc->over_d_xyz[0] *
                     (  - data[( (x_i-1) + y_i * 4 ) * 4 + z_i ]
                        + data[( (x_i+1) + y_i * 4 ) * 4 + z_i ] ),
                     is2d ? 0 :
                     .5*vc->over_d_xyz[1] *
                     (  - data[( x_i + (y_i-1) * 4 ) * 4 + z_i ]
                        + data[( x_i + (y_i+1) * 4 ) * 4 + z_i ] ),
                     .5*vc->over_d_xyz[2] *
                     (  - data[( x_i + y_i * 4 ) * 4 + (z_i-1) ]
                        + data[( x_i + y_i * 4 ) * 4 + (z_i+1) ] ) );

    }

}


static void fill_cell_hess_3pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h) {

    int i;

    for (i=0; i<8; i++) {

        /* WARN: ASSUMES THAT THE CELL IS PRECISELY THE 3PT STENCIL CELL. 
           can be made generic by replacing: 
           --> 1 with vc->bp_xyz[0, 1 or 2] as appropriate 
           --> 4 with vc->n_xyz[0, 1 or 2] as appropriate 
        */

        const int x_i = 1 + cube_order[i][0];
        const int y_i = 1 + cube_order[i][1];
        const int z_i = 1 + cube_order[i][2];

        if ( (is2d) && (i>=4)) return;

        m3d_sym_init(h+i, 
                         vc->over_d_xyz[0] * vc->over_d_xyz[0] * /* xx */
                         (  +  data[( (x_i-1) + y_i * 4 ) * 4 + z_i] 
                            -2*data[( (x_i  ) + y_i * 4 ) * 4 + z_i] 
                            +  data[( (x_i+1) + y_i * 4 ) * 4 + z_i] ), 
                         is2d ? 0 :
                         vc->over_d_xyz[0] * vc->over_d_xyz[1] * /* xy */
                         (  + data[( (x_i-1) + (y_i-1) * 4 ) * 4 + z_i] 
                            - data[( (x_i+1) + (y_i-1) * 4 ) * 4 + z_i] 
                            - data[( (x_i-1) + (y_i+1) * 4 ) * 4 + z_i]
                            + data[( (x_i+1) + (y_i+1) * 4 ) * 4 + z_i] ), 
                         vc->over_d_xyz[0] * vc->over_d_xyz[2] * /* xz */
                         (  + data[( (x_i-1) + y_i * 4 ) * 4 + (z_i-1)] 
                            - data[( (x_i-1) + y_i * 4 ) * 4 + (z_i+1)] 
                            - data[( (x_i+1) + y_i * 4 ) * 4 + (z_i-1)] 
                            + data[( (x_i+1) + y_i * 4 ) * 4 + (z_i+1)] ), 
                         is2d ? 0 :
                         vc->over_d_xyz[1] * vc->over_d_xyz[1] * /* yy */
                         (  +  data[( x_i + (y_i-1) * 4 ) * 4 + z_i] 
                            -2*data[( x_i + (y_i  ) * 4 ) * 4 + z_i] 
                            +  data[( x_i + (y_i+1) * 4 ) * 4 + z_i] ), 
                         is2d ? 0 :
                         vc->over_d_xyz[1] * vc->over_d_xyz[2] * /* yz */
                         (  + data[( x_i + (y_i-1) * 4 ) * 4 + (z_i-1)] 
                            - data[( x_i + (y_i-1) * 4 ) * 4 + (z_i+1)] 
                            - data[( x_i + (y_i+1) * 4 ) * 4 + (z_i-1)] 
                            + data[( x_i + (y_i+1) * 4 ) * 4 + (z_i+1)] ), 
                         vc->over_d_xyz[2] * vc->over_d_xyz[2] * /* zz */
                         (  +  data[( x_i + y_i * 4 ) * 4 + (z_i-1)] 
                            -2*data[( x_i + y_i * 4 ) * 4 + (z_i  )] 
                            +  data[( x_i + y_i * 4 ) * 4 + (z_i+1)] ) ); 
    }

}

static void fill_cell_grad_5pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, v3d_t *g) {

    int i;

    for (i=0; i<8; i++) {

        /* WARN: ASSUMES THAT THE CELL IS PRECISELY THE 5PT STENCIL CELL. 
           can be made generic by replacing: 
           --> 2 with vc->bp_xyz[0, 1 or 2] as appropriate 
           --> 6 with vc->n_xyz[0, 1 or 2] as appropriate 
        */

        const int x_i = 2 + cube_order[i][0];
        const int y_i = 2 + cube_order[i][1];
        const int z_i = 2 + cube_order[i][2];

        if ( (is2d) && (i>=4)) return;

        v3d_init(g+i, 
                     (1.0/12.0)*vc->over_d_xyz[0] * 
                     (+   data[( (x_i-2) + y_i * 6 ) * 6  + z_i] 
                      - 8*data[( (x_i-1) + y_i * 6 ) * 6  + z_i] 
                      + 8*data[( (x_i+1) + y_i * 6 ) * 6  + z_i] 
                      -   data[( (x_i+2) + y_i * 6 ) * 6  + z_i] ), 
                     is2d ? 0 :
                     (1.0/12.0)*vc->over_d_xyz[1] * 
                     (+   data[( x_i + (y_i-2) * 6 ) * 6  + z_i] 
                      - 8*data[( x_i + (y_i-1) * 6 ) * 6  + z_i] 
                      + 8*data[( x_i + (y_i+1) * 6 ) * 6  + z_i] 
                      -   data[( x_i + (y_i+2) * 6 ) * 6  + z_i] ), 
                     (1.0/12.0)*vc->over_d_xyz[2] * 
                     (+   data[( x_i + y_i * 6 ) * 6  + (z_i-2)] 
                      - 8*data[( x_i + y_i * 6 ) * 6  + (z_i-1)] 
                      + 8*data[( x_i + y_i * 6 ) * 6  + (z_i+1)] 
                      -   data[( x_i + y_i * 6 ) * 6  + (z_i+2)] ) );
    }

}



static void fill_cell_hess_5pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h) {

    int i;

    for (i=0; i<8; i++) {

        /* WARN: ASSUMES THAT THE CELL IS PRECISELY THE 5PT STENCIL CELL. 
           can be made generic by replacing: 
           --> 2 with vc->bp_xyz[0, 1 or 2] as appropriate 
           --> 6 with vc->n_xyz[0, 1 or 2] as appropriate 
        */

        const int x_i = 2 + cube_order[i][0];
        const int y_i = 2 + cube_order[i][1];
        const int z_i = 2 + cube_order[i][2];

        if ( (is2d) && (i>=4)) return;

        m3d_sym_init(h+i, 
                         (1.0/12.0)*vc->over_d_xyz[0] * vc->over_d_xyz[0] * /* xx */
                         (-    data[( (x_i-2) + y_i * 6 ) * 6 + z_i] 
                          + 16*data[( (x_i-1) + y_i * 6 ) * 6 + z_i] 
                          - 30*data[( (x_i  ) + y_i * 6 ) * 6 + z_i] 
                          + 16*data[( (x_i+1) + y_i * 6 ) * 6 + z_i] 
                          -    data[( (x_i+2) + y_i * 6 ) * 6 + z_i] ), 
                         is2d ? 0 :
                         (1.0/144.0)*vc->over_d_xyz[0] * vc->over_d_xyz[1] * /* xy */
                         (  +  8 * (+ data[( (x_i + 1) + (y_i-2) * 6 ) * 6 + z_i] 
                                    + data[( (x_i + 2) + (y_i-1) * 6 ) * 6 + z_i] 
                                    + data[( (x_i - 2) + (y_i+1) * 6 ) * 6 + z_i]
                                    + data[( (x_i - 1) + (y_i+2) * 6 ) * 6 + z_i] )
                            -  8 * (+ data[( (x_i - 1) + (y_i-2) * 6 ) * 6 + z_i] 
                                    + data[( (x_i - 2) + (y_i-1) * 6 ) * 6 + z_i] 
                                    + data[( (x_i + 2) + (y_i+1) * 6 ) * 6 + z_i]
                                    + data[( (x_i + 1) + (y_i+2) * 6 ) * 6 + z_i] )
                            +      (- data[( (x_i + 2) + (y_i-2) * 6 ) * 6 + z_i] 
                                    - data[( (x_i - 2) + (y_i+2) * 6 ) * 6 + z_i] 
                                    + data[( (x_i - 2) + (y_i-2) * 6 ) * 6 + z_i]
                                    + data[( (x_i + 2) + (y_i+2) * 6 ) * 6 + z_i] )
                            + 64 * (- data[( (x_i - 1) + (y_i-1) * 6 ) * 6 + z_i] 
                                    - data[( (x_i + 1) + (y_i+1) * 6 ) * 6 + z_i] 
                                    + data[( (x_i + 1) + (y_i-1) * 6 ) * 6 + z_i]
                                    + data[( (x_i - 1) + (y_i+1) * 6 ) * 6 + z_i] )
                            ), 
                         (1.0/144.0)*vc->over_d_xyz[0] * vc->over_d_xyz[2] * /* xz */
                         (  +  8 * (+ data[( (x_i+1) + y_i * 6 ) * 6 + (z_i-2)] 
                                    + data[( (x_i+2) + y_i * 6 ) * 6 + (z_i-1)] 
                                    + data[( (x_i-2) + y_i * 6 ) * 6 + (z_i+1)]
                                    + data[( (x_i-1) + y_i * 6 ) * 6 + (z_i+2)] )
                            -  8 * (+ data[( (x_i-1) + y_i * 6 ) * 6 + (z_i-2)] 
                                    + data[( (x_i-2) + y_i * 6 ) * 6 + (z_i-1)] 
                                    + data[( (x_i+2) + y_i * 6 ) * 6 + (z_i+1)]
                                    + data[( (x_i+1) + y_i * 6 ) * 6 + (z_i+2)] )
                            +      (- data[( (x_i+2) + y_i * 6 ) * 6 + (z_i-2)] 
                                    - data[( (x_i-2) + y_i * 6 ) * 6 + (z_i+2)] 
                                    + data[( (x_i-2) + y_i * 6 ) * 6 + (z_i-2)]
                                    + data[( (x_i+2) + y_i * 6 ) * 6 + (z_i+2)] )
                            + 64 * (- data[( (x_i-1) + y_i * 6 ) * 6 + (z_i-1)] 
                                    - data[( (x_i+1) + y_i * 6 ) * 6 + (z_i+1)] 
                                    + data[( (x_i+1) + y_i * 6 ) * 6 + (z_i-1)]
                                    + data[( (x_i-1) + y_i * 6 ) * 6 + (z_i+1)] )
                            ), 
                         is2d ? 0 :
                         (1.0/12.0)*vc->over_d_xyz[1] * vc->over_d_xyz[1] * /* yy */
                         (-    data[( x_i + (y_i-2) * 6 ) * 6 + z_i] 
                          + 16*data[( x_i + (y_i-1) * 6 ) * 6 + z_i] 
                          - 30*data[( x_i + (y_i  ) * 6 ) * 6 + z_i] 
                          + 16*data[( x_i + (y_i+1) * 6 ) * 6 + z_i] 
                          -    data[( x_i + (y_i+2) * 6 ) * 6 + z_i] ), 
                         is2d ? 0 :
                         (1.0/144.0)*vc->over_d_xyz[1] * vc->over_d_xyz[2] * /* yz */
                         (  +  8 * (+ data[( x_i + (y_i-2) * 6 ) * 6 + (z_i+1)] 
                                    + data[( x_i + (y_i-1) * 6 ) * 6 + (z_i+2)] 
                                    + data[( x_i + (y_i+1) * 6 ) * 6 + (z_i-2)]
                                    + data[( x_i + (y_i+2) * 6 ) * 6 + (z_i-1)] )
                            -  8 * (+ data[( x_i + (y_i-2) * 6 ) * 6 + (z_i-1)] 
                                    + data[( x_i + (y_i-1) * 6 ) * 6 + (z_i-2)] 
                                    + data[( x_i + (y_i+1) * 6 ) * 6 + (z_i+2)]
                                    + data[( x_i + (y_i+2) * 6 ) * 6 + (z_i+1)] )
                            +      (- data[( x_i + (y_i-2) * 6 ) * 6 + (z_i+2)] 
                                    - data[( x_i + (y_i+2) * 6 ) * 6 + (z_i-2)] 
                                    + data[( x_i + (y_i-2) * 6 ) * 6 + (z_i-2)]
                                    + data[( x_i + (y_i+2) * 6 ) * 6 + (z_i+2)] )
                            + 64 * (- data[( x_i + (y_i-1) * 6 ) * 6 + (z_i-1)] 
                                    - data[( x_i + (y_i+1) * 6 ) * 6 + (z_i+1)] 
                                    + data[( x_i + (y_i-1) * 6 ) * 6 + (z_i-1)]
                                    + data[( x_i + (y_i+1) * 6 ) * 6 + (z_i+1)] )
                            ), 
                         (1.0/12.0)*vc->over_d_xyz[2] * vc->over_d_xyz[2] * /* zz */
                         (-    data[(x_i + y_i * 6 ) * 6 + (z_i-2)] 
                          + 16*data[(x_i + y_i * 6 ) * 6 + (z_i-1)] 
                          - 30*data[(x_i + y_i * 6 ) * 6 + (z_i  )] 
                          + 16*data[(x_i + y_i * 6 ) * 6 + (z_i+1)] 
                          -    data[(x_i + y_i * 6 ) * 6 + (z_i+2)] ) );
    }

}


static void fill_cell_grad_7pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, v3d_t *g) {

    int i;

    for (i=0; i<8; i++) {
        /* WARN: ASSUMES THAT THE CELL IS PRECISELY THE 5PT STENCIL CELL. 
           can be made generic by replacing: 
           --> 3 with vc->bp_xyz[0, 1 or 2] as appropriate 
           --> 8 with vc->n_xyz[0, 1 or 2] as appropriate 
        */

        const int x_i = 3 + cube_order[i][0];
        const int y_i = 3 + cube_order[i][1];
        const int z_i = 3 + cube_order[i][2];

        if ( (is2d) && (i>=4)) return;

        v3d_init(g+i, 
                     (1.0/60.0)*vc->over_d_xyz[0] * 
                     (-    data[( (x_i-3) + y_i * 8 ) * 8 + z_i] 
                      +  9*data[( (x_i-2) + y_i * 8 ) * 8 + z_i] 
                      - 15*data[( (x_i-1) + y_i * 8 ) * 8 + z_i] 
                      + 15*data[( (x_i+1) + y_i * 8 ) * 8 + z_i] 
                      -  9*data[( (x_i+2) + y_i * 8 ) * 8 + z_i] 
                      +    data[( (x_i+3) + y_i * 8 ) * 8 + z_i] ), 
                     is2d ? 0 :
                     (1.0/60.0)*vc->over_d_xyz[1] * 
                     (-    data[( x_i + (y_i-3) * 8 ) * 8 + z_i] 
                      +  9*data[( x_i + (y_i-2) * 8 ) * 8 + z_i] 
                      - 15*data[( x_i + (y_i-1) * 8 ) * 8 + z_i] 
                      + 15*data[( x_i + (y_i+1) * 8 ) * 8 + z_i] 
                      -  9*data[( x_i + (y_i+2) * 8 ) * 8 + z_i] 
                      +    data[( x_i + (y_i+3) * 8 ) * 8 + z_i] ), 
                     (1.0/60.0)*vc->over_d_xyz[2] * 
                     (-    data[( x_i + y_i * 8 ) * 8 + (z_i-3)] 
                      +  9*data[( x_i + y_i * 8 ) * 8 + (z_i-2)] 
                      - 15*data[( x_i + y_i * 8 ) * 8 + (z_i-1)] 
                      + 15*data[( x_i + y_i * 8 ) * 8 + (z_i+1)] 
                      -  9*data[( x_i + y_i * 8 ) * 8 + (z_i+2)] 
                      +    data[( x_i + y_i * 8 ) * 8 + (z_i+3)] ) );

    }

}

static void fill_cell_hess_7pt(const se_vel_approx_cell_t *vc, const float *data, int is2d, m3d_t *h) {

    // who uses a hessian anyway?!?!
    fill_cell_hess_3pt(vc, data, is2d, h);

    
    return;
}
