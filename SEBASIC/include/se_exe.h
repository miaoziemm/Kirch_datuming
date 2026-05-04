#ifndef SE_EXE_H
#define SE_EXE_H

#include "se_type.h"
#include "se_alloc.h"
#include "se_util.h"
#include "se_assert.h"

int se_using_exe_from_path(void);
const char* se_get_exe_name(void);
const char* se_get_exe_path(void);
const char* se_get_exe_absolute_path(void);

void se_set_exe(const int   exe_from_path,
                 const char* exe_name,
                 const char* exe_path,
                 const char* exe_absolute_path);

void se_exe_reset(void);


#endif