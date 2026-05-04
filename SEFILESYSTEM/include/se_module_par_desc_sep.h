#ifndef SE_MODULE_PAR_DESC_SEP_H
#define SE_MODULE_PAR_DESC_SEP_H
#include <SEBASIC/include/se_basic.h>
#include "se_fs_sep.h"


void se_module_par_desc_get_from_sep(se_module_par_desc_t *mpd,
                                      const sep_t *sep);

void se_module_par_desc_save_in_sep(sep_t *sep,
                                     const se_module_par_desc_t *mpd,
                                     int save_defaults);

void se_module_par_desc_get(se_module_par_desc_t *mpd);  

#endif