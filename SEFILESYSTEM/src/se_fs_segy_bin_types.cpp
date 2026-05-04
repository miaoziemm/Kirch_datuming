#include "../include/se_fs_segy_bin_types.h"


static const char* SEGY_BIN_TYPE_NAME_LAST    = "last";
static const char* SEGY_BIN_TYPE_NAME_CLOSEST = "closest";
static const char* SEGY_BIN_TYPE_NAME_AVERAGE = "average";

const char* segy_bin_type_to_string(segy_bin_type_t type)
{
    switch(type) {
    case SEGY_BIN_TYPE_LAST:
        return SEGY_BIN_TYPE_NAME_LAST;
    case SEGY_BIN_TYPE_CLOSEST:
        return SEGY_BIN_TYPE_NAME_CLOSEST;
    case SEGY_BIN_TYPE_AVERAGE:
        return SEGY_BIN_TYPE_NAME_AVERAGE;
    }
    return "unknown";
}


