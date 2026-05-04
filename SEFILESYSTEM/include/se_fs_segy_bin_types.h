#ifndef SE_FS_SEGY_BIN_TYPES_H
#define SE_FS_SEGY_BIN_TYPES_H

#include <stdlib.h>
#include <string.h>

typedef enum {
    SEGY_BIN_TYPE_LAST = 0,
    SEGY_BIN_TYPE_CLOSEST = 1,
    SEGY_BIN_TYPE_AVERAGE = 2
} segy_bin_type_t;

const char* segy_bin_type_to_string(segy_bin_type_t type);

// segy_bin_type_t segy_string_to_bin_type(const char* name);

// segy_bin_type_t segy_bin_type_from_par(const char* parname, const char* pargroup);

#endif