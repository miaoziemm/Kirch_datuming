#include "../include/se_vel_cfg.h"


void vel_cfg_get_def(vel_cfg_t* p)
{
    se_module_par_desc_set_def(p->par_desc);
}

int vel_cfgs_are_equiv(const vel_cfg_t *A, const vel_cfg_t *B)
{
    return se_module_par_desc_pars_are_equiv(A->par_desc, B->par_desc);
}

void vel_cfg_validate(const vel_cfg_t *p)
{

    se_module_par_desc_verify_required(p->par_desc);

    switch (p->vel_stencil) {
    case 3:
    case 5:
    case 7:
        break;
    default: 
        ERROR(("The interpolation stencil ( interpolation_stencil= ) must be 3, 5, or 7."));
    } 

}

void vel_cfg_report(const vel_cfg_t *p)
{

    INFOV((0,""));
    INFOV((0,"** VELOCITY: **"));
    INFOV((0,"--> Velocity-> [%s]", p->vel_file));
    if (p->aniso) {
        INFOV((0,"--> %s Anisotropic migration", p->tti?"TTI":"VTI"));
        INFOV((0,"--> Epsilon-> [%s]", p->eps_file));
        INFOV((0,"--> Delta-> [%s]", p->del_file));
        if (p->tti) {
            INFOV((0,"--> Theta X-> [%s]", p->tetx_file));
            if (p->tety_file) INFOV((0,"--> Theta Y-> [%s]", p->tety_file));
            INFOV((0,"--> Theta files contain %s",p->tet_dip?"Dips":"Angles (in degrees)"));
        }
    }
    INFOV((0,"--> Interpolation stencil size is %d", p->vel_stencil));
    if (p->detect_vel_bdry) {
        INFOV((0,"--> Detecting boundries, threshold %g", p->vel_bdry_frac));
    }

}

vel_cfg_t *vel_cfg_create(void)
{

    vel_cfg_t *cfg;

    cfg = (vel_cfg_t *)malloc(sizeof(vel_cfg_t));
    memset(cfg, 0, sizeof(vel_cfg_t));

    cfg->par_desc = (se_module_par_desc_t*)malloc(sizeof(se_module_par_desc_t));
    memset(cfg->par_desc, 0, sizeof(se_module_par_desc_t));

    se_module_par_desc_add_s(cfg->par_desc, &cfg->vel_file,  (char*)"vel",         NULL, 1);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->aniso,     (char*)"anisotropic", 0,    0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->eps_file,  (char*)"epsilon",     NULL, 0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->del_file,  (char*)"delta",       NULL, 0);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->tti,       (char*)"tti",         0,    0);
    se_module_par_desc_add_i(cfg->par_desc, &cfg->tet_dip,   (char*)"theta_dip",   1,    0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->tetx_file, (char*)"theta_x",     NULL, 0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->tety_file, (char*)"theta_y",     NULL, 0);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->force_3d, (char*)"force_3d", 0, 0);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->vel_stencil,     (char*)"interpolation_stencil",      3, 0);
    se_module_par_desc_add_i(cfg->par_desc, &cfg->detect_vel_bdry, (char*)"detect_velocity_boundaries", 1, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->vel_bdry_frac,   (char*)"velocity_boundary_fraction", 0.25, 0);

    return cfg;

}

void vel_cfg_init(vel_cfg_t* cfg)
{
    if (!cfg) return;
    memset(cfg, 0, sizeof(vel_cfg_t));

    cfg->par_desc = (se_module_par_desc_t*)malloc(sizeof(se_module_par_desc_t));
    memset(cfg->par_desc, 0, sizeof(se_module_par_desc_t));

    se_module_par_desc_add_s(cfg->par_desc, &cfg->vel_file,  (char*)"vel",         NULL, 1);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->aniso,     (char*)"anisotropic", 0,    0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->eps_file,  (char*)"epsilon",     NULL, 0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->del_file,  (char*)"delta",       NULL, 0);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->tti,       (char*)"tti",         0,    0);
    se_module_par_desc_add_i(cfg->par_desc, &cfg->tet_dip,   (char*)"theta_dip",   1,    0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->tetx_file, (char*)"theta_x",     NULL, 0);
    se_module_par_desc_add_s(cfg->par_desc, &cfg->tety_file, (char*)"theta_y",     NULL, 0);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->force_3d, (char*)"force_3d", 0, 0);

    se_module_par_desc_add_i(cfg->par_desc, &cfg->vel_stencil,     (char*)"interpolation_stencil",      3, 0);
    se_module_par_desc_add_i(cfg->par_desc, &cfg->detect_vel_bdry, (char*)"detect_velocity_boundaries", 1, 0);
    se_module_par_desc_add_d(cfg->par_desc, &cfg->vel_bdry_frac,   (char*)"velocity_boundary_fraction", 0.25, 0);
}
