#include "../include/se_vel_cfg_sep.h"

void vel_cfg_get_from_sep(vel_cfg_t* p, const sep_t *sep)
{
    se_module_par_desc_get_from_sep(p->par_desc, sep);
}

void vel_cfg_save_in_sep(sep_t* sep, const vel_cfg_t* p, int save_defaults)
{
    se_module_par_desc_save_in_sep(sep, p->par_desc, save_defaults);
}

void vel_cfg_get_from_par(vel_cfg_t* p)
{
    se_module_par_desc_get(p->par_desc);
}



static void velio_load_aniso(const vel_cfg_t *vel_cfg, 
                                  const grid3d_t* velgrid, 
                                  float **del_data, float **eta_data, 
                                  float **ttix_data, float **ttiy_data,
                                  int is3d);

void vel_init_and_load_from_sep(se_vel_t *vel, 
                                    const vel_cfg_t *vel_cfg)
{

    sep_t* vel_sep = NULL;
    grid3d_t* velgrid;

    float* vel_data = NULL;
    float* del_data = NULL;
    float* eta_data = NULL;
    float* ttix_data = NULL;
    float* ttiy_data = NULL;

    int is3d;

    vel_sep = sep_open(vel_cfg->vel_file, SEP_READ, 1);

    velgrid = init_grid3d_from_sepheaders(vel_sep);

    if (vel_sep->headers->n[2]>1 || vel_cfg->force_3d ) {
        is3d=1;
    } else {
        is3d=0;
    }

    sep_close(vel_sep);

    INFOV((0,"Reading velocity ... "));

    vel_data = sep_read_whole_file(vel_cfg->vel_file);

    if (vel_cfg->aniso) {
        velio_load_aniso(vel_cfg, velgrid, &del_data, &eta_data, &ttix_data, &ttiy_data, is3d);
    } 

    se_vel_init(vel, vel_data, del_data, eta_data, ttix_data, ttiy_data, velgrid, is3d); 

    destroy_grid3d(velgrid);

}


static int velio_sep_has_different_cube(const grid3d_t* velgrid, const char* sep_filename) 
{
    sep_t* tmp_sep = NULL;
    grid3d_t* tmp_grid;
    int ret = 0;

    tmp_sep = sep_open(sep_filename, SEP_READ, 1);
    tmp_grid = init_grid3d_from_sepheaders(tmp_sep);
    if (!grid3d_equal(velgrid, tmp_grid)) ret=1;
    destroy_grid3d(tmp_grid);
    sep_close(tmp_sep);
    
    return ret;
}



static void velio_load_aniso(const vel_cfg_t *vel_cfg, 
                                  const grid3d_t* velgrid, 
                                  float **del_data, float **eta_data, 
                                  float **ttix_data, float **ttiy_data,
                                  int is3d) 
{

    int64_t i, n;
    
    n = ((int64_t)velgrid->x.n)*((int64_t)velgrid->y.n)*((int64_t)velgrid->z.n);

    if ( velio_sep_has_different_cube(velgrid, vel_cfg->del_file) ) 
        ERROR(("Velocity and Delta input SEP files do not have the same grid."));
    INFOV((1,"Reading delta ... "));
    (*del_data) = sep_read_whole_file(vel_cfg->del_file);
    
    if ( velio_sep_has_different_cube(velgrid, vel_cfg->eps_file) )
        ERROR(("Velocity and Epsilon input SEP files do not have the same grid."));
    INFOV((1,"Reading epsilon ... "));
    (*eta_data) = sep_read_whole_file(vel_cfg->eps_file);
    for (i=0; i<n; i++) {
        (*eta_data)[i] -= (*del_data)[i];
    }

    if (vel_cfg->tti) {
        if ( velio_sep_has_different_cube(velgrid, vel_cfg->tetx_file) )
            ERROR(("Velocity and crossline tilt input SEP files do not have the same grid."));
        INFOV((1,"Reading crossline tilt ... "));
        (*ttix_data) = sep_read_whole_file(vel_cfg->tetx_file);
        if (is3d) {
            if ( velio_sep_has_different_cube(velgrid, vel_cfg->tety_file) )
                ERROR(("Velocity and inline tilt input SEP files do not have the same grid."));
            INFOV((1,"Reading inline tilt ... "));
            (*ttiy_data) = sep_read_whole_file(vel_cfg->tety_file);
        } else {
            (*ttiy_data)=NULL;
        }

        if (!vel_cfg->tet_dip) {
            if (*ttix_data) {
                for (i=0; i<n; i++) {
                    (*ttix_data)[i] = (float)tan(-(*ttix_data)[i]*M_PI/180.0);
                }
            }
            if (*ttiy_data) {
                for (i=0; i<n; i++) {
                    (*ttiy_data)[i] = (float)tan(-(*ttiy_data)[i]*M_PI/180.0);
                }
            }
        }
        
    } else {
        (*ttix_data)=NULL;
        (*ttiy_data)=NULL;
    }
}


