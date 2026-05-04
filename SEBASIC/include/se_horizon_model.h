#ifndef SE_HORIZON_MODEL_H
#define SE_HORIZON_MODEL_H
#include "se_sparse_interpolation.h"
#include "se_hash.h"
#include "se_array.h"
#include "se_rsf.h"
#include "se_grid.h"
#include "se_interpolation.h"

/* horizon line */
typedef struct
{
    char*   name;
    int     npoints;
    double* coords;
} hrz_line_t;

void init_hrz_line (hrz_line_t* l);
void clear_hrz_line(hrz_line_t* l);


/* horizon surface */
typedef struct
{
    char*     name;
    se_hash  attrs;
    array lines;
} hrz_surface_t;

void init_hrz_surface (hrz_surface_t* s);
void clear_hrz_surface(hrz_surface_t* s);


/* horizon model */
typedef struct
{
    char*     name;
    se_hash  attrs;
    int       ndim;
    array surfaces;
} hrz_model_t;

void init_hrz_model (hrz_model_t* m);
void clear_hrz_model(hrz_model_t* m);

// return a value that is out of bounds of ax for use as sparse interpolation no_value
float hrz_no_value_from_axa(const axa_t *ax);

/* Return value and set nearest grid indices for a point with coordinates from hrz_line_t
 * Called internally by rsf_from_hrz
 */
double rsf_get_grid_indices(int mod_ndim, int axis_indep, int gathers, const double *coords,
                                const grid3d_t *grid, int64_t *gz, int64_t *gx, int64_t *gy);
int rsf_get_grid_dim(const grid3d_t *grid);

/* Fill rsf with values from from the horizon with index idx in mod.
 * Allocates rsf->prvt, fills with no_value
 */
void rsf_from_hrz(const hrz_model_t* mod, int idx, const grid3d_t *grid, rsf3d_t *rsf, int axis_indep, int gathers, float no_value);

/* layer model */
typedef struct {
    double ztop, zbottom;
    rsf3d_t *rsftop, *rsfbottom;
} hrz_layer_t;

/* At a given (x,y) location, this sets *top and *bottom to the value
   of the top most and bottom most layer limit depth. If there are
   more than one layers, all layers lie between *top and *bottom. */
void hrz_layer_get_limits(const hrz_layer_t* layers, int N, double x, double y, double *top, double *bottom);

/* Returns true if the point (x,y,z) is within one or more of the N layers, and false otherwise. */
int hrz_layer_check_pt(const hrz_layer_t* layers, int N, double z, double x, double y);

/* Creates rsftop from first horizon in top, rsfbottom from first horizon in bottom.
 * Assumes layer is allocated, but not initialized; does not set ztop, zbottom.
 */
void hrz_layer_from_first_hz(const hrz_model_t *top, const hrz_model_t *bottom, const grid3d_t* grid, hrz_layer_t *layer);

/* Create a layer between horizons topfile,bottomfile on grid, with min and max depth ztop,zbottom. */
hrz_layer_t *hrz_single_constraining_layer(const char *topfile, const char *bottomfile, double ztop, double zbottom, const grid3d_t* grid);

void clear_hrz_layer(hrz_layer_t *l);

inline double hrz_get_coord(const hrz_model_t* m,
                                    const hrz_line_t* l,
                                    int point_idx,
                                    int coord_idx)
{
    ASSERT(point_idx >= 0 && point_idx < l->npoints);
    ASSERT(coord_idx >= 0 && coord_idx < m->ndim);

    return l->coords[point_idx * m->ndim + coord_idx];
}

inline void hrz_set_coord(const hrz_model_t* m,
                                  hrz_line_t* l,
                                  int point_idx,
                                  int coord_idx,
                                  double val)
{
    ASSERT(point_idx >= 0 && point_idx < l->npoints);
    ASSERT(coord_idx >= 0 && coord_idx < m->ndim);

    l->coords[point_idx * m->ndim + coord_idx] = val;
}

const char* hrz_get_attr  (se_hash attrs, const char* key);
int         hrz_get_attr_i(se_hash attrs, const char* key);
double      hrz_get_attr_d(se_hash attrs, const char* key);



#endif