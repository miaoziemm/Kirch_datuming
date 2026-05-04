#include "../include/se_exe.h"

static int   exe_from_path = 0;
static char* exe_name = NULL;
static char* exe_path = NULL;
static char* exe_absolute_path = NULL;

int se_using_exe_from_path(void)
{
    return exe_from_path;
}

const char* se_get_exe_absolute_path(void)
{
    return exe_absolute_path;
}

const char* se_get_exe_name(void)
{
    return exe_name;
}

const char* se_get_exe_path(void)
{
    return exe_path;
}

void se_exe_reset(void)
{
    exe_from_path = 0;

    if(exe_name) {
        free(exe_name);
        exe_name = NULL;
    }

    if(exe_path) {
        free(exe_path);
        exe_path = NULL;
    }

    if(exe_absolute_path) {
        free(exe_absolute_path);
        exe_absolute_path = NULL;
    }
}

void se_set_exe(const int   _exe_from_path,
                 const char* _exe_name,
                 const char* _exe_path,
                 const char* _exe_absolute_path)
{
    se_exe_reset();

    exe_from_path = _exe_from_path;
    exe_name = _exe_name ? se_strdup(_exe_name) : NULL;
    exe_path = _exe_path ? se_strdup(_exe_path) : NULL;
    exe_absolute_path = _exe_absolute_path ? se_strdup(_exe_absolute_path) : NULL;
}
