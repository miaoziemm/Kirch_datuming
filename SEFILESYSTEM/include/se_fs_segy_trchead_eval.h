#ifndef SE_FS_SEGY_TRCHEAD_EVAL_H
#define SE_FS_SEGY_TRCHEAD_EVAL_H

#include "se_fs_segy_header.h"

typedef void* trchead_eval_handle;

trchead_eval_handle init_trchead_eval( const char* expression );

void destroy_trchead_eval( trchead_eval_handle eh );

void cleanup_trchead_eval_globals( void );

double trchead_eval( trchead_eval_handle eh, const trchead* trch );

double trchead_evalx( trchead_eval_handle eh, const trchead* trch,
                          double lsx, double lsy, double lgx, double lgy);

#endif