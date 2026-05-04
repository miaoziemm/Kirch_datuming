
/**
 * Tracing status.
 */
typedef struct ode_stat_s {

    double dt; /**< Current ODE time step */
    int dir;

} ode_stat_t;


/**
 * Implementation of an Euler ODE integrator for dU/dt=f(U).
 *
 * \param[in] ode_cfg the ODE configuration
 * \param[in] v_app the velocity approximation object
 * \param[in] U_0 the initial U
 * \param[in] T the final time at which to stop ODE integration
 * \param[out] U_T U at T
 * \param[out] ode_stat the ODE status at T
 *
 */
static void ode_euler(const ode_cfg_t *ode_cfg,
                      se_vel_approx_t* v_app,
                      const ode_U_t* U_0,
                      double T,
                      ode_U_t* U_T,
                      ode_stat_t *ode_stat);

/**
 * Implementation of a Runge-Kutta 23 ODE integrator for dU/dt=f(U).
 *
 * \param[in] ode_cfg the ODE configuration
 * \param[in] dt_init the initial time step
 * \param[in] v_app the velocity approximation object
 * \param[in] U_0 the initial U
 * \param[in] T the final time at which to stop ODE integration
 * \param[out] U_T the Gaussian beam at t_fin
 * \param[out] ode_stat the ODE status at T
 *
 */
static void ode_rk23(const ode_cfg_t *ode_cfg,
                     double dt_init,
                     se_vel_approx_t* v_app,
                     const ode_U_t* gb,
                     double t_fin,
                     ode_U_t* gb_fin,
                     ode_stat_t *ode_stat);


ode_return_t ode_trace_tbl_fct_name(const ode_cfg_t *ode_cfg,
                                        se_vel_approx_t* v_app,
                                        const ode_U_bnds_t *bnds, 
                                        const ode_U_t *U, double to_t,
                                        ode_U_t *U_tbl,
                                        size_t tbl_len,
                                        size_t *tbl_min_i, size_t *tbl_max_i)
{
    int dir;
    double tbl_t, tbl_dt;
    int64_t tbl_i=0;
    ode_stat_t ode_stat;
    ode_return_t exit_code;

    tbl_dt = fabs(to_t - ode_U_get_t(U) ) / ( (tbl_len>1) ? (tbl_len-1):1 ) ;
    tbl_t = ode_U_get_t(U);

    if ( ode_U_get_t(U) > to_t ) {
        dir = -1;
        tbl_i=tbl_len-1;
    } else {
        dir = 1;
        tbl_i = 0;
    }

    *tbl_max_i = tbl_i;
    *tbl_min_i = tbl_i;

    U_tbl[tbl_i] = *U;
    tbl_t += dir*tbl_dt;
    tbl_i += dir;

    switch (ode_cfg->solver_type) {
        case ODE_TYPE_RK23:
            ode_stat.dt = ode_cfg->max_dt;
            break;
        case ODE_TYPE_EULER:
            ode_stat.dt = ode_cfg->fixed_dt;
            break;
        default:
            exit_code = ODE_UNKNOWN;
            goto ode_cleanup;
    }
    ode_stat.dir = dir;

    ode_U_enforce_constraints(U_tbl + tbl_i - dir, v_app);

    while ( (tbl_i<(int64_t)tbl_len) && (tbl_i>=0) ) {

        switch (ode_cfg->solver_type) {
        case ODE_TYPE_RK23:

            ode_rk23(ode_cfg, ode_stat.dt, v_app, U_tbl+tbl_i-dir, tbl_t, U_tbl+tbl_i, &ode_stat);

            break;

        case ODE_TYPE_EULER:
            
            ode_euler(ode_cfg, v_app, U_tbl+tbl_i-dir, tbl_t, U_tbl+tbl_i, &ode_stat);

            break;

        default:
            exit_code = ODE_UNKNOWN;
            goto ode_cleanup;
            break;
        }
        
        if (dir<0) {
            *tbl_min_i=tbl_i;
        } else {
            *tbl_max_i=tbl_i;
        }

        if ( ( bnds ) && 
             ( !ode_U_is_in_bnds(bnds, ode_U_get_xyz(U_tbl+tbl_i)) ) ) {
            
            exit_code = ODE_ATBNDRY;
            goto ode_cleanup;

        }
        
        tbl_t += dir*tbl_dt;
        tbl_i += dir;

    }

    exit_code = ODE_SUCCESS;

 ode_cleanup:

    return exit_code;
}


static void ode_euler(const ode_cfg_t *ode_cfg,
                      se_vel_approx_t* v_app,
                      const ode_U_t* U_0,
                      double T,
                      ode_U_t* U_T,
                      ode_stat_t *ode_stat)
{

    int dir;
    double dt;
    ode_U_t U_t, U_n, rhs;
    int done;

    U_t = *U_0;

    dt = ode_cfg->fixed_dt;

    ode_stat->dt = dt;
    dir = ode_stat->dir;

    done = 0;

    while ( !done ) {

        if ( dir * ( T - ode_U_get_t(&U_t) ) < dt) { 
            dt = dir * ( T - ode_U_get_t(&U_t) );
            done = 1;
        }

        ode_eval_rhs(ode_cfg, v_app, &U_t, &rhs);
        ode_U_add_scaled(&U_n, 1, &U_t, dir*dt, &rhs);
        ode_U_enforce_constraints(&U_n, v_app);
        
        U_t = U_n;
        
    }

    *U_T = U_t;

}

/**
 * Runge-Kutta 23 inner computation. The next value is computed in two
 * ways and the time step is decreased if the difference between them is
 * too large.
 \f{eqnarray*}{
 y^{\text{higher}}_{n+1} &=& y_n + \frac{2}{ 9}k_0 + \frac{1}{3}k_1 + \frac{4}{9}k_2 \\
 y^{\text{lower }}_{n+1} &=& y_n + \frac{7}{24}k_0 + \frac{1}{4}k_1 + \frac{1}{3}k_2 + \frac{1}{8}k_3 \\
 k_0 &=& dt * rhs\left(y_n\right) \\
 k_1 &=& dt * rhs\left(y_n + \frac{1}{2}*k_0 \right) \\
 k_2 &=& dt * rhs\left(y_n + \frac{3}{4}*k_1 \right) \\
 k_3 &=& dt * rhs\left(y_n + \frac{2}{9}*k_0 + \frac{1}{3}*k_1 + \frac{4}{9}k_2 \right)
 \f}
 *
 */
static void ode_rk23_inner(const ode_cfg_t *ode_cfg,
                           double *dt, int dir,
                           se_vel_approx_t *v_app,
                           const ode_U_t *U_t,
                           ode_U_t *U_n,
                           int *repeat_step)
{
    ode_U_t rhs_k[4];
    ode_U_t U_lo;
    double err;
    double dir_dt_72 = dir * (*dt) / 72.0;

    /* k0 */
    ode_eval_rhs(ode_cfg, v_app, U_t, rhs_k+0);

    /* k1 */
    ode_U_add_scaled(&U_lo, 1, U_t, dir_dt_72*36 /* 1.0/2.0 */, rhs_k+0);
    ode_eval_rhs(ode_cfg, v_app, &U_lo, rhs_k+1);

    /* k2 */
    ode_U_add_scaled(&U_lo, 1, U_t, dir_dt_72*54 /* 3.0/4.0 */, rhs_k+1);
    ode_eval_rhs(ode_cfg, v_app, &U_lo, rhs_k+2);

    /* y_n+1 higher */
    ode_U_add_scaled(  U_n, 1, U_t, dir_dt_72*16 /* 2.0/9.0 */, rhs_k+0);
    ode_U_accum_scaled(U_n,         dir_dt_72*24 /* 3.0/9.0 */, rhs_k+1);
    ode_U_accum_scaled(U_n,         dir_dt_72*32 /* 4.0/9.0 */, rhs_k+2);

    /* k3 */
    ode_eval_rhs(ode_cfg, v_app, U_n, rhs_k+3);

    /* y_n+1 lower */
    ode_U_add_scaled(&U_lo, 1, U_t, dir_dt_72*21 /* 7.0/24.0 */, rhs_k+0);
    ode_U_accum_scaled(&U_lo,       dir_dt_72*18 /* 6.0/24.0 */, rhs_k+1);
    ode_U_accum_scaled(&U_lo,       dir_dt_72*24 /* 8.0/24.0 */, rhs_k+2);
    ode_U_accum_scaled(&U_lo,       dir_dt_72*9  /* 3.0/24.0 */, rhs_k+3);

    ode_U_accum_scaled(&U_lo, -1, U_n);

    err = ode_U_norm_rel(&U_lo, U_n);

    if ( err > ode_cfg->max_err_tol ) {
        if ( (*dt) > (ode_cfg->min_dt) ) {
            (*repeat_step) = 1;
            (*dt) = fmax((*dt)*cbrt(.5*(ode_cfg->min_err_tol + ode_cfg->max_err_tol)/err),
                             ode_cfg->min_dt);
        } else {
            (*repeat_step)=0;
            ode_U_enforce_constraints(U_n, v_app);
        }
    } else {
        (*repeat_step)=0;
        if (err < ode_cfg->min_err_tol) {
            (*dt) = fmin((*dt)*cbrt(.5*(ode_cfg->min_err_tol + ode_cfg->max_err_tol)/err),
                             ode_cfg->max_dt);
        }
        ode_U_enforce_constraints(U_n, v_app);
    }

}

static void ode_rk23(const ode_cfg_t *ode_cfg,
                     double dt_init,
                     se_vel_approx_t* v_app,
                     const ode_U_t* U_0,
                     double T,
                     ode_U_t* U_T,
                     ode_stat_t *ode_stat)
{
    int dir;
    double dt;
    ode_U_t U_t, U_n;
    int done;
        
    U_t = *U_0;

    dir = ode_stat->dir;

    dt = dt_init;

    ode_stat->dt = dt;

    done = 0;
    while ( !done ) {

        int repeat_step;

        do {
            if ( dir * ( T - ode_U_get_t(&U_t) ) < dt ) {
                dt = dir * ( T - ode_U_get_t(&U_t) );
                done = 1;
            } else {
                ode_stat->dt = dt;
                done = 0;
            }

            ode_rk23_inner(ode_cfg, &dt, dir, v_app, &U_t, &U_n, &repeat_step);

        } while (repeat_step);

        U_t = U_n;

    }

    *U_T = U_t;

}
