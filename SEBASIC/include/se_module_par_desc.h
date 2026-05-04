#ifndef SE_MODULE_PAR_DESC_H
#define SE_MODULE_PAR_DESC_H
#include "se_alloc.h"
#include "se_type.h"
#include "se_log.h"
#include "se_basic_math.h"
#include "se_util.h"


typedef struct se_d_par_desc_s {
    double *par_val_ptr;
    char *par_name;
    double par_def;
    int is_set;
    int is_required;
} se_d_par_desc_t;

typedef struct se_i_par_desc_s{
    int *par_val_ptr;
    char *par_name;
    int par_def;
    int is_set;
    int is_required;
} se_i_par_desc_t;

typedef struct se_i64_par_desc_s{
    int64_t *par_val_ptr;
    char *par_name;
    int64_t par_def;
    int is_set;
    int is_required;
} se_i64_par_desc_t;

typedef struct se_s_par_desc_s{
    char **par_val_ptr;
    char *par_name;
    char *par_def;
    int is_set;
    int is_required;
} se_s_par_desc_t;


typedef struct se_module_par_desc_s {

    se_i_par_desc_t *i_pars;
    int n_i_pars;

    se_i64_par_desc_t *i64_pars;
    int n_i64_pars;

    se_d_par_desc_t *d_pars;
    int n_d_pars;

    se_s_par_desc_t *s_pars;
    int n_s_pars;

} se_module_par_desc_t;

void se_module_par_desc_add_s(se_module_par_desc_t *mpd,
                               char **par_val_ptr, char *par_name, char* par_def, int is_required);
void se_module_par_desc_add_i(se_module_par_desc_t *mpd,
                               int *par_val_ptr, char *par_name, int par_def, int is_required);
void se_module_par_desc_add_i64(se_module_par_desc_t *mpd,
                                 int64_t *par_val_ptr, char *par_name, int64_t par_def, int is_required);
void se_module_par_desc_add_d(se_module_par_desc_t *mpd,
                               double *par_val_ptr, char *par_name, double par_def, int is_required);

int se_module_par_desc_pars_are_equiv(se_module_par_desc_t *A, se_module_par_desc_t *B);

void se_module_par_desc_verify_required(se_module_par_desc_t *mpd);

void se_module_par_desc_set_def(se_module_par_desc_t *mpd);

void se_module_par_desc_destroy(se_module_par_desc_t *mpd);



#endif