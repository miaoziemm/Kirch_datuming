#ifndef SE_VEL_CFG_SEP_H
#define SE_VEL_CFG_SEP_H
#include <SEBASIC/include/se_basic.h>
#include "se_fs_sep.h"
#include "se_module_par_desc_sep.h"

void vel_cfg_get_from_sep(vel_cfg_t* p, const sep_t *sep);

void vel_cfg_save_in_sep(sep_t* sep, const vel_cfg_t* p, int save_defaults);


void vel_cfg_get_from_par(vel_cfg_t* p);

void vel_init_and_load_from_sep(se_vel_t *vel, 
                                    const vel_cfg_t *vel_cfg); 

#endif