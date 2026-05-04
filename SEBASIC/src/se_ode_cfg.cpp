#include "../include/se_ode_cfg.h"


void ode_cfg_get_def(ode_cfg_t* p)
{
    se_module_par_desc_set_def(p->par_desc);

    p->solver_type = ode_cfg_solver_type_from_id(p->solver_type_id);
    p->wf_weights[1] = 1 - p->wf_weights[0];
}

int ode_cfgs_are_equiv(const ode_cfg_t *A, const ode_cfg_t *B)
{
    return se_module_par_desc_pars_are_equiv(A->par_desc, B->par_desc);
}

void ode_cfg_validate(const ode_cfg_t *cfg)
{

    se_module_par_desc_verify_required(cfg->par_desc);

    switch (cfg->solver_type) {
    case ODE_TYPE_RK23:

        if ( cfg->min_dt < 0 ) 
            ERROR(("ODE minimum step must be positive."));
        if ( cfg->min_err_tol < 0 ) 
            ERROR(("ODE minimum error tolerance must be non-negative."));
        if ( cfg->max_err_tol < cfg->min_err_tol ) 
            ERROR(("ODE minimum error tolerance must not be bigger than the maximum error tolerance."));
        if ( cfg->max_dt < cfg->min_dt ) 
            ERROR(("ODE minimum step must not be bigger than the maximum step."));
        break;

    case ODE_TYPE_EULER:

        if ( cfg->fixed_dt < 0 ) 
            ERROR(("ODE step must be positive."));
        break;

    default:

        ERROR(("Unknown ray tracer type."));

    }

    if ( cfg->use_wf ) {

        if ( ( cfg->wf_threshold < 0 ) || ( cfg->wf_threshold > 1 ) ) 
            ERROR(("ODE wavefront threshold must be between 0 and 1."));

        if ( cfg->wf_size <= 0 ) 
            ERROR(("ODE wavefront size must be greater than 0."));

        if ( cfg->wf_sets < 1 ) 
            ERROR(("number of ODE wavefront concentric sets must be greater than or equal to 1."));

        if ( ( cfg->wf_weights[0] > 1 ) || ( cfg->wf_weights[0] < 0 ) )
            ERROR(("ODE ray straightening must be between 0 and 1."));
    }
}
//options
void ode_cfg_set_rk23_opts(ode_cfg_t *cfg,
                            double min_err_tol, double max_err_tol, 
                            double min_dt, double max_dt) 
{
    cfg->min_err_tol=min_err_tol;
    cfg->max_err_tol=max_err_tol;

    cfg->min_dt=min_dt;
    cfg->max_dt=max_dt;

    cfg->solver_type = ODE_TYPE_RK23;
    cfg->solver_type_id = ode_cfg_solver_type_id(cfg->solver_type);
}

void ode_cfg_set_euler_opts(ode_cfg_t *cfg, double dt) 
{
    cfg->fixed_dt=dt;

    cfg->solver_type=ODE_TYPE_EULER;
    cfg->solver_type_id = ode_cfg_solver_type_id(cfg->solver_type);
}

void ode_cfg_set_wf_opts(ode_cfg_t *cfg, 
                          int use_wf, double wf_threshold, double wf_size, int wf_sets,
                          double center_pt_weight) 
{
    cfg->use_wf=use_wf;

    if (use_wf) {

        cfg->wf_threshold=wf_threshold;
        cfg->wf_size=wf_size;
        cfg->wf_sets=wf_sets;
        cfg->wf_weights[0]=center_pt_weight;
        cfg->wf_weights[1]=1-center_pt_weight;

    } else {

        cfg->wf_size=1;
        cfg->wf_threshold=1;
        cfg->wf_sets=1;
        cfg->wf_weights[0]=0;
        cfg->wf_weights[1]=1;

    }
}


static const char * ode_get_tracer_name(const ode_cfg_t *cfg) {
    
    switch (cfg->solver_type) {
    case ODE_TYPE_RK23:
        return "Adaptive step";
        break;
    case ODE_TYPE_EULER:
        return "Fixed step";
        break;
    default:
        return "Unknown tracer type";
        break;
    }
}


void ode_cfg_report(const ode_cfg_t *cfg)
{
    INFOV((0," "));
    INFOV((0,"** ODE SOLVER CONFIGURATION: **"));
    INFOV((0,"--> %s", ode_get_tracer_name(cfg)));
    switch (cfg->solver_type) {
    case ODE_TYPE_RK23:
        INFOV((0,"--> ODE Configuration: min/max dt ( %g , %g )", 
               cfg->min_dt,cfg->max_dt));
        INFOV((0,"                       min/max error tolerance ( %g , %g )", 
               cfg->min_err_tol, cfg->max_err_tol));
        break;
    case ODE_TYPE_EULER:
        INFOV((0,"--> ODE Configuration: dt is %g", cfg->fixed_dt));
        break;
    default:
        INFOV((0,"UNKNOWN ODE TRACER TYPE\n"));
        break;
    }
    if (cfg->use_wf) {
        INFOV((0,"--> using wavefront approximation near salt"));
        INFOV((0,"          velocity threshold %g ", cfg->wf_threshold));
        INFOV((0,"          wavefront size %g with %d sets", 
               cfg->wf_size, cfg->wf_sets));
        INFOV((0,"          ray straightening %g", cfg->wf_weights[0])); 
    }
}


ode_cfg_t *ode_cfg_create(void)
{

    ode_cfg_t *cfg;

    cfg = (ode_cfg_t *)malloc(sizeof(ode_cfg_t));
    memset(cfg, 0, sizeof(ode_cfg_t));

    cfg->par_desc = (se_module_par_desc_t *)malloc(sizeof(se_module_par_desc_t));
    memset(cfg->par_desc, 0, sizeof(se_module_par_desc_t));


    se_module_par_desc_add_i(cfg->par_desc, &cfg->solver_type_id, (char *)"ode_solver", 0, 0);

    // for rk23
    se_module_par_desc_add_d(cfg->par_desc, &cfg->min_dt,      (char *)"ode_min_dt",      1e-8, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->max_dt,      (char *)"ode_max_dt",      2e-3, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->min_err_tol, (char *)"ode_min_err_tol", 1e-8, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->max_err_tol, (char *)"ode_max_err_tol", 1e-6, 0);

    // for euler
    se_module_par_desc_add_d(cfg->par_desc, &cfg->fixed_dt, (char *)"ode_dt", 2e-3, 0);

    // for wf near salt
    se_module_par_desc_add_i(cfg->par_desc, &cfg->use_wf,        (char *)"ode_use_wf",       0, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->wf_threshold,  (char *)"ode_wf_threshold", 0.25, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->wf_size,       (char *)"ode_wf_size",      25, 0);
    se_module_par_desc_add_i(cfg->par_desc, &cfg->wf_sets,       (char *)"ode_wf_sets",      1, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->wf_weights[0], (char *)"ode_wf_ctr_pt_w",  0, 0);

    return cfg;

}

void ode_cfg_init(ode_cfg_t* cfg)
{
    if (!cfg) return;
    memset(cfg, 0, sizeof(ode_cfg_t));

    cfg->par_desc = (se_module_par_desc_t *)malloc(sizeof(se_module_par_desc_t));
    memset(cfg->par_desc, 0, sizeof(se_module_par_desc_t));

    se_module_par_desc_add_i(cfg->par_desc, &cfg->solver_type_id, (char *)"ode_solver", 0, 0);

    // for rk23
    se_module_par_desc_add_d(cfg->par_desc, &cfg->min_dt,      (char *)"ode_min_dt",      1e-8, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->max_dt,      (char *)"ode_max_dt",      2e-3, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->min_err_tol, (char *)"ode_min_err_tol", 1e-8, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->max_err_tol, (char *)"ode_max_err_tol", 1e-6, 0);

    // for euler
    se_module_par_desc_add_d(cfg->par_desc, &cfg->fixed_dt, (char *)"ode_dt", 2e-3, 0);

    // for wf near salt
    se_module_par_desc_add_i(cfg->par_desc, &cfg->use_wf,        (char *)"ode_use_wf",       0, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->wf_threshold,  (char *)"ode_wf_threshold", 0.25, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->wf_size,       (char *)"ode_wf_size",      25, 0);
    se_module_par_desc_add_i(cfg->par_desc, &cfg->wf_sets,       (char *)"ode_wf_sets",      1, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->wf_weights[0], (char *)"ode_wf_ctr_pt_w",  0, 0);
}
