#ifndef SE_ODE_CFG_SEP_H
#define SE_ODE_CFG_SEP_H
#include <SEBASIC/include/se_basic.h>
#include "se_fs_sep.h"
#include "se_module_par_desc_sep.h"


void ode_cfg_get_from_sep(ode_cfg_t* p, const sep_t *sep);

void ode_cfg_save_in_sep(sep_t* sep, ode_cfg_t* p, int save_defaults);

void ode_cfg_get_from_par(ode_cfg_t* p);

#endif