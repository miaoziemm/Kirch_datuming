#ifndef SESI_MAP_H
#define SESI_MAP_H

#include <SEBASIC/include/se_basic.h>
#include "se_fs_sep.h"
#include "se_par_sep.h"

typedef void* sesi_map_t;

sesi_map_t sesi_map_create(const char* file, const axa_t* axis, int naxis,
                             int include_all, int grow_by);
sesi_map_t sesi_map_open(const char* file, int read_only);

void sesi_map_set_cache(sesi_map_t _zsm, int64_t max_cache_size, size_t hashsize);

void sesi_map_close(sesi_map_t _zsm);
int64_t sesi_map_add_from_coords(sesi_map_t _zsm, int idx, const double* coords);
int64_t sesi_map_add_from_index(sesi_map_t _zsm, int idx, int64_t lidx);

int64_t sesi_map_set_indexes_from_coords(sesi_map_t _zsm, const double* coords, 
                                          const int* idx, int nidx);
int64_t sesi_map_set_indexes_from_index(sesi_map_t _zsm,int64_t lidx,
                                         const int* idx, int nidx);

int64_t sesi_map_add_indexes_from_index(sesi_map_t _zsm,int64_t lidx,
                                        const int* idx, int nidx);

void sesi_map_compact(sesi_map_t _zsm);

int sesi_map_get_naxis(sesi_map_t _zsm);
axa_t sesi_map_get_axis(sesi_map_t _zsm, int i);

int64_t sesi_map_get_number_of_non_empty_cells(sesi_map_t _zsm);

int64_t sesi_map_get_nelements(sesi_map_t _zsm);
int sesi_map_have_hdr_int(sesi_map_t _zsm, const char* name);
int sesi_map_get_hdr_int(sesi_map_t _zsm, const char* name, int def);
int sesi_map_have_hdr_float(sesi_map_t _zsm, const char* name);
double sesi_map_get_hdr_float(sesi_map_t _zsm, const char* name, double def);
int sesi_map_have_hdr(sesi_map_t _zsm, const char* name);
char* sesi_map_get_hdr(sesi_map_t _zsm, const char* name, const char* def);

void sesi_map_set_header(sesi_map_t _zsm, const char* name, const char* value);
void sesi_map_set_header_int(sesi_map_t _zsm, const char* name, int value);
void sesi_map_set_header_float(sesi_map_t _zsm, const char* name, double value);

void sesi_map_set_header_survey_description(sesi_map_t _zsm,
                                             const survey_description_t* sd,
                                             const char* prefix,
                                             int use_only_set_fields);
void sesi_map_get_header_survey_description(sesi_map_t _zsm,
                                             survey_description_t* sd,
                                             const char* prefix,
                                             const survey_description_t* project_sd);

int sesi_map_get_fold_from_coords(sesi_map_t _zsm, 
                                   const double* coords);
int sesi_map_get_fold_from_index(sesi_map_t _zsm, 
                                  int64_t lidx);

void sesi_map_get_indexes_from_index(sesi_map_t _zsm, 
                                      int64_t lidx,
                                      int** list, int* n, int* list_size);
void sesi_map_get_indexes_from_coords(sesi_map_t _zsm, 
                                       const double* coords,
                                       int** list, int* n, int* list_size);

void sesi_map_extract_piece(sesi_map_t _zsm, const char* new_map,
                             int* limits);

void sesi_map_print_cache_stat(sesi_map_t _zsm, int v);

int sesi_map_is_good_map_header(const char* file);

#endif