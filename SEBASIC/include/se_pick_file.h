#ifndef SE_PICK_FILE_H
#define SE_PICK_FILE_H
#include "se_array.h"
#include "se_horizon_model.h"
#include "se_hash.h"
#include "se_util.h"


/* .pick file attribute names */
extern const char* const PICKKEY_NHZ;
extern const char* const PICKKEY_OILINE;
extern const char* const PICKKEY_OXLINE;
extern const char* const PICKKEY_OZ;
extern const char* const PICKKEY_O4;
extern const char* const PICKKEY_AZIMUTH;
extern const char* const PICKKEY_COLS;
extern const char* const PICKKEY_AXES;

/* .pick file horizon attribute names */
extern const char* const PICKHKEY_LOAD;
extern const char* const PICKHKEY_COLOR;
extern const char* const PICKHKEY_FREEZE;
extern const char* const PICKHKEY_EXTENDED;

/**
 * Calls parse_pick_file, printing any error or warning messages.
 */
hrz_model_t* parse_pick_file_print(const char* filename);


/**
 * Parses a .pick horizon pick file and creates a 'hrz_model_t' object.
 * If not NULL, 'warns' may be populated with warning messages (caller fn should free it).
 *
 * On success returns the new model;
 * On failure sets 'error' (the caller should free it) and returns NULL.
 *
 * Note that if the .pick file has 4 columns, the order of columns 3 and 4 is swapped.
 */
hrz_model_t* parse_pick_file(const char* filename, char** error, array* warns);


/**
 * Writes a .pick horizon pick file from a 'hrz_model_t' model.
 *
 * On success returns OK;
 * On failure sets 'error' (the caller should free it) and returns ERROR.
 */
int write_pick_file(const hrz_model_t* model, const char* filename, char** error);



#endif