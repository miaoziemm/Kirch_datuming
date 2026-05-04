#include "../include/se_ode_cfg_sep.h"


void ode_cfg_get_from_sep(ode_cfg_t* p, const sep_t *sep)
{
    se_module_par_desc_get_from_sep(p->par_desc, sep);

    p->solver_type = ode_cfg_solver_type_from_id(p->solver_type_id);
    p->wf_weights[1] = 1 - p->wf_weights[0];
}

void ode_cfg_save_in_sep(sep_t* sep, ode_cfg_t* p, int save_defaults)
{
    p->solver_type_id = ode_cfg_solver_type_id(p->solver_type);

    se_module_par_desc_save_in_sep(sep, p->par_desc, save_defaults);
}


void ode_cfg_get_from_par(ode_cfg_t* p)
{
    se_module_par_desc_get(p->par_desc);

    p->solver_type = ode_cfg_solver_type_from_id(p->solver_type_id);
    p->wf_weights[1] = 1 - p->wf_weights[0];
}

