#ifndef SE_VEL_CFG_H
#define SE_VEL_CFG_H
#include "se_module_par_desc.h"


/**
 * Beam depth migration configuration.
 */
typedef struct vel_cfg_s {

    char *vel_file; /**< Velocity filename. */
    char *eps_file; /**< Epsilon filename. */
    char *del_file; /**< Delta filename. */
    char *tetx_file; /**< Theta X filename (for TTI). */ //Theta
    char *tety_file; /**< Theta Y filename (for TTI). */  //Theta

    int aniso; /**< Flag for anisotropic migration. */  
    int tti; /**< Flag for TTI anisotropic migration. */  
    int tet_dip; /**< Flag for TTI Theta files contain dips or angles (degrees). */  

    int force_3d; /**< Flag for forcing 3D. */

    int vel_stencil; /**< Velocity interpolation stencil half size. */ 
    int detect_vel_bdry; /**< Flag to detect large velocity changes. */
    double vel_bdry_frac; /**< Fractional threshold for detecting boundaries (max_vel-min_vel)/max_vel */ //分数阈值 探测边界

    se_module_par_desc_t *par_desc; 

} vel_cfg_t;

vel_cfg_t *vel_cfg_create(void);

// Initialize a stack-allocated vel_cfg_t. Allocates and sets up par_desc.
// Use this when vel_cfg_t is embedded in another struct (no heap allocation for vel_cfg_t itself).
void vel_cfg_init(vel_cfg_t* p);

void vel_cfg_get_def(vel_cfg_t* p);

int vel_cfgs_are_equiv(const vel_cfg_t *A, const vel_cfg_t *B);

void vel_cfg_validate(const vel_cfg_t *p);

void vel_cfg_report(const vel_cfg_t *p);

inline void vel_cfg_cleanup(vel_cfg_t* p) 
{
    if (p) {
        free(p->vel_file);
        free(p->eps_file);
        free(p->del_file);
        free(p->tetx_file);
        free(p->tety_file);
        
        se_module_par_desc_destroy(p->par_desc);
        
        // Don't free p itself - for stack allocated objects
        p->vel_file = NULL;
        p->eps_file = NULL;
        p->del_file = NULL;
        p->tetx_file = NULL;
        p->tety_file = NULL;
        p->par_desc = NULL;
    }
}

inline void vel_cfg_destroy(vel_cfg_t* p) 
{
    if (p) {
        vel_cfg_cleanup(p);
        free(p);
    }
}



#endif