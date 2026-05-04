#include "../include/se_module_par_desc_sep.h"


void se_module_par_desc_get_from_sep(se_module_par_desc_t *mpd,
                                      const sep_t *sep)
{
    int i;
    for (i=0; i<mpd->n_i_pars; i++) {
        if ( sep_have_hdr( sep, mpd->i_pars[i].par_name )) {
            (*mpd->i_pars[i].par_val_ptr) = sep_get_hdr_int(sep, mpd->i_pars[i].par_name, mpd->i_pars[i].par_def);
            mpd->i_pars[i].is_set = 1;
        } else {
            (*mpd->i_pars[i].par_val_ptr) = mpd->i_pars[i].par_def;
            mpd->i_pars[i].is_set = 0;
        }
    }
    for (i=0; i<mpd->n_i64_pars; i++) {
        if ( sep_have_hdr( sep, mpd->i64_pars[i].par_name )) {
            (*mpd->i64_pars[i].par_val_ptr) = sep_get_hdr_int(sep, mpd->i64_pars[i].par_name, (int)(mpd->i64_pars[i].par_def));
            mpd->i64_pars[i].is_set = 1;
        } else {
            (*mpd->i64_pars[i].par_val_ptr) = mpd->i64_pars[i].par_def;
            mpd->i64_pars[i].is_set = 0;
        }
    }
    for (i=0; i<mpd->n_d_pars; i++) {
        if ( sep_have_hdr( sep, mpd->d_pars[i].par_name )) {
            (*mpd->d_pars[i].par_val_ptr) = sep_get_hdr_float(sep, mpd->d_pars[i].par_name, mpd->d_pars[i].par_def);
            mpd->d_pars[i].is_set = 1;
        } else {
            (*mpd->d_pars[i].par_val_ptr) = mpd->d_pars[i].par_def;
            mpd->d_pars[i].is_set = 0;
        }
    }
    for (i=0; i<mpd->n_s_pars; i++) {
        if ( sep_have_hdr( sep, mpd->s_pars[i].par_name )) {
            (*mpd->s_pars[i].par_val_ptr) = sep_get_hdr(sep, mpd->s_pars[i].par_name, mpd->s_pars[i].par_def);
            mpd->s_pars[i].is_set = 1;
        } else {
            if (mpd->s_pars[i].par_def) {
                (*mpd->s_pars[i].par_val_ptr) = se_strdup(mpd->s_pars[i].par_def);
            } else {
                (*mpd->s_pars[i].par_val_ptr) = NULL;
            }
            mpd->s_pars[i].is_set = 0;
        }
    }
}


void se_module_par_desc_save_in_sep(sep_t *sep,
                                     const se_module_par_desc_t *mpd,
                                     int save_defaults)
{

    int i;
    for (i=0; i<mpd->n_i_pars; i++) {
        if (save_defaults || mpd->i_pars[i].is_set) {
            sep_set_header_int(sep,
                                  mpd->i_pars[i].par_name,
                                  *mpd->i_pars[i].par_val_ptr);
        }
    }
    for (i=0; i<mpd->n_i64_pars; i++) {
        if (save_defaults || mpd->i64_pars[i].is_set) {
            sep_set_header_int(sep,
                                  mpd->i64_pars[i].par_name,
                                  (int)(*mpd->i64_pars[i].par_val_ptr));
        }
    }
    for (i=0; i<mpd->n_d_pars; i++) {
        if (save_defaults || mpd->d_pars[i].is_set) {
            sep_set_header_float(sep,
                                    mpd->d_pars[i].par_name,
                                    *mpd->d_pars[i].par_val_ptr);
        }
    }
    for (i=0; i<mpd->n_s_pars; i++) {
        if (save_defaults || mpd->s_pars[i].is_set) {
            if (*mpd->s_pars[i].par_val_ptr) {
                sep_set_header(sep,
                                  mpd->s_pars[i].par_name,
                                  *mpd->s_pars[i].par_val_ptr);
            }
        }
    }
    
}


void se_module_par_desc_get(se_module_par_desc_t *mpd)
{
    int i;
    for (i=0; i<mpd->n_i_pars; i++) {
        if ( se_have_par( mpd->i_pars[i].par_name )) {
            (*mpd->i_pars[i].par_val_ptr) = se_get_par_int( mpd->i_pars[i].par_name );
            mpd->i_pars[i].is_set = 1;
        } else {
            (*mpd->i_pars[i].par_val_ptr) = mpd->i_pars[i].par_def;
            mpd->i_pars[i].is_set = 0;
        }
    }
    for (i=0; i<mpd->n_i64_pars; i++) {
        if ( se_have_par( mpd->i64_pars[i].par_name )) {
            (*mpd->i64_pars[i].par_val_ptr) = se_get_par_int64( mpd->i64_pars[i].par_name );
            mpd->i64_pars[i].is_set = 1;
        } else {
            (*mpd->i64_pars[i].par_val_ptr) = mpd->i64_pars[i].par_def;
            mpd->i64_pars[i].is_set = 0;
        }
    }
    for (i=0; i<mpd->n_d_pars; i++) {
        if ( se_have_par( mpd->d_pars[i].par_name )) {
            (*mpd->d_pars[i].par_val_ptr) = se_get_par_double( mpd->d_pars[i].par_name );
            mpd->d_pars[i].is_set = 1;
        } else {
            (*mpd->d_pars[i].par_val_ptr) = mpd->d_pars[i].par_def;
            mpd->d_pars[i].is_set = 0;
        }
    }
    for (i=0; i<mpd->n_s_pars; i++) {
        if ( se_have_par( mpd->s_pars[i].par_name )) {
            (*mpd->s_pars[i].par_val_ptr) = se_get_par_str( mpd->s_pars[i].par_name );
            mpd->s_pars[i].is_set = 1;
        } else {
            if (mpd->s_pars[i].par_def) {
                (*mpd->s_pars[i].par_val_ptr) = strdup(mpd->s_pars[i].par_def);
            } else {
                (*mpd->s_pars[i].par_val_ptr) = NULL;
            }
            mpd->s_pars[i].is_set = 0;
        }
    }
}
