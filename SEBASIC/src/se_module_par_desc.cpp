#include "../include/se_module_par_desc.h"


void se_module_par_desc_add_s(se_module_par_desc_t *mpd,
                               char **par_val_ptr, char *par_name, char* par_def, 
                               int is_required)
{
    mpd->s_pars = realloc1type(se_s_par_desc_t,mpd->s_pars, mpd->n_s_pars+1);

    mpd->s_pars[mpd->n_s_pars].par_val_ptr = par_val_ptr;
    mpd->s_pars[mpd->n_s_pars].par_name = par_name;
    mpd->s_pars[mpd->n_s_pars].par_def = par_def;
    mpd->s_pars[mpd->n_s_pars].is_required = is_required;
    mpd->n_s_pars++;
}

void se_module_par_desc_add_i(se_module_par_desc_t *mpd,
                               int *par_val_ptr, char *par_name, int par_def,
                               int is_required)
{
    mpd->i_pars = realloc1type(se_i_par_desc_t,mpd->i_pars, mpd->n_i_pars+1);



    mpd->i_pars[mpd->n_i_pars].par_val_ptr = par_val_ptr;
    mpd->i_pars[mpd->n_i_pars].par_name = par_name;
    mpd->i_pars[mpd->n_i_pars].par_def = par_def;
    mpd->i_pars[mpd->n_i_pars].is_required = is_required;
    mpd->n_i_pars++;
}

void se_module_par_desc_add_i64(se_module_par_desc_t *mpd,
                                 int64_t *par_val_ptr, char *par_name, int64_t par_def,
                                 int is_required)
{
    mpd->i64_pars = realloc1type(se_i64_par_desc_t,mpd->i64_pars, mpd->n_i64_pars+1 );

    mpd->i64_pars[mpd->n_i64_pars].par_val_ptr = par_val_ptr;
    mpd->i64_pars[mpd->n_i64_pars].par_name = par_name;
    mpd->i64_pars[mpd->n_i64_pars].par_def = par_def;
    mpd->i64_pars[mpd->n_i64_pars].is_required = is_required;
    mpd->n_i64_pars++;
}

void se_module_par_desc_add_d(se_module_par_desc_t *mpd,
                               double *par_val_ptr, char *par_name, double par_def,
                               int is_required)
{
    mpd->d_pars = realloc1type(se_d_par_desc_t,mpd->d_pars, mpd->n_d_pars+1);

    mpd->d_pars[mpd->n_d_pars].par_val_ptr = par_val_ptr;
    mpd->d_pars[mpd->n_d_pars].par_name = par_name;
    mpd->d_pars[mpd->n_d_pars].par_def = par_def;
    mpd->d_pars[mpd->n_d_pars].is_required = is_required;
    mpd->n_d_pars++;
}

int se_module_par_desc_pars_are_equiv(se_module_par_desc_t *A, se_module_par_desc_t *B)
{
    int i;

    if (A->n_i_pars != B->n_i_pars) return 0;
    for (i=0; i < A->n_i_pars; i++) {
        if ( strcmp(   A->i_pars[i].par_name,         B->i_pars[i].par_name     ) ) return 0;
        if ( !      ( (*A->i_pars[i].par_val_ptr) == (*B->i_pars[i].par_val_ptr) ) ) return 0;
    }

    if (A->n_i64_pars != B->n_i64_pars) return 0;
    for (i=0; i < A->n_i64_pars; i++) {
        if (  strcmp(   A->i64_pars[i].par_name,         B->i64_pars[i].par_name     ) ) return 0;
        if ( !      ( (*A->i64_pars[i].par_val_ptr) == (*B->i64_pars[i].par_val_ptr) ) ) return 0;
    }

    if (A->n_d_pars != B->n_d_pars) return 0;
    for (i=0; i < A->n_d_pars; i++) {
        if (  strcmp(             A->d_pars[i].par_name,              B->d_pars[i].par_name     ) ) return 0;
        if ( !fequal((float)(*A->d_pars[i].par_val_ptr), (float)(*B->d_pars[i].par_val_ptr) ) ) return 0;
    }

    if (A->n_s_pars != B->n_s_pars) return 0;
    for (i=0; i < A->n_s_pars; i++) {
        if (  strcmp(   A->s_pars[i].par_name,          B->s_pars[i].par_name     ) ) return 0;
        if (       ( !(*A->s_pars[i].par_val_ptr) &&  (*B->s_pars[i].par_val_ptr) ) ) return 0;
        if (       (  (*A->s_pars[i].par_val_ptr) && !(*B->s_pars[i].par_val_ptr) ) ) return 0;
    }

    return 1;

}

void se_module_par_desc_verify_required(se_module_par_desc_t *mpd)
{
    int i;
    for (i=0; i<mpd->n_i_pars; i++) {
        if (mpd->i_pars[i].is_required && !mpd->i_pars[i].is_set) {
            ERROR(("Parameter %s must be specified!", mpd->i_pars[i].par_name));
        }
    }
    for (i=0; i<mpd->n_i64_pars; i++) {
        if (mpd->i64_pars[i].is_required && !mpd->i64_pars[i].is_set) {
            ERROR(("Parameter %s must be specified!", mpd->i64_pars[i].par_name));
        }
    }
    for (i=0; i<mpd->n_d_pars; i++) {
        if (mpd->d_pars[i].is_required && !mpd->d_pars[i].is_set) {
            ERROR(("Parameter %s must be specified!", mpd->d_pars[i].par_name));
        }
    }
    for (i=0; i<mpd->n_s_pars; i++) {
        if (mpd->s_pars[i].is_required && !mpd->s_pars[i].is_set) {
            ERROR(("Parameter %s must be specified!", mpd->s_pars[i].par_name));
        }
    }
}

void se_module_par_desc_set_def(se_module_par_desc_t *mpd)
{
    int i;
    for (i=0; i<mpd->n_i_pars; i++) {
        (*mpd->i_pars[i].par_val_ptr) = mpd->i_pars[i].par_def;
        mpd->i_pars[i].is_set = 0;
    }
    for (i=0; i<mpd->n_i64_pars; i++) {
        (*mpd->i64_pars[i].par_val_ptr) = mpd->i64_pars[i].par_def;
        mpd->i64_pars[i].is_set = 0;
    }
    for (i=0; i<mpd->n_d_pars; i++) {
        (*mpd->d_pars[i].par_val_ptr) = mpd->d_pars[i].par_def;
        mpd->d_pars[i].is_set = 0;
    }
    for (i=0; i<mpd->n_s_pars; i++) {
        if (mpd->s_pars[i].par_def) {
            (*mpd->s_pars[i].par_val_ptr) = se_strdup(mpd->s_pars[i].par_def);
        } else {
            (*mpd->s_pars[i].par_val_ptr) = NULL;
        }
        mpd->s_pars[i].is_set = 0;
    }
}

void se_module_par_desc_destroy(se_module_par_desc_t *mpd)
{
    if (mpd) {
        // Free string parameter values first
        int i;
        for (i = 0; i < mpd->n_s_pars; i++) {
            if (mpd->s_pars[i].par_val_ptr && *mpd->s_pars[i].par_val_ptr) {
                free1char(*mpd->s_pars[i].par_val_ptr);
                *mpd->s_pars[i].par_val_ptr = NULL;
            }
        }
        
        free(mpd->i_pars);
        free(mpd->i64_pars);
        free(mpd->d_pars);
        free(mpd->s_pars);

        free(mpd);
    }
}
