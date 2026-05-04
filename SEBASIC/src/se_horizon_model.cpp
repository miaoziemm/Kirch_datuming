#include "../include/se_horizon_model.h"
#include "../include/se_pick_file.h"

double rsf_get_grid_indices(int mod_ndim, int axis_indep, int gathers, const double *coords,
                                const grid3d_t *grid, int64_t *gz, int64_t *gx, int64_t *gy) {
    const int ndim=rsf_get_grid_dim(grid);
    double x,y=grid->y.o,z; //z,h,x,y or z,x,y from image/gathers
    double v=0.0,vx=0.0,vy=0.0; //x,y,z from grid
    *gz = 0;

    if (gathers) {
        double h=coords[0], vz=0.0;
        x=coords[1];
        if (mod_ndim==4) {
            y=coords[2];
            z=coords[3];
        }
        else {
            z=coords[2];
        }
        switch (axis_indep) {
        case 1: // z
            v=z; vx=x; vy=y; vz=h;
            break;
        case 2: // h
            v=h; vx=x; vy=y; vz=z;
            break;
        case 3: // x
            v=x; vx=h; vy=y; vz=z;
            break;
        case 4: // y
            v=y; vx=h; vy=x; vz=z;
            break;
        default:
            ERROR(("Invalid value for axis_indep=%d",axis_indep));
            break;
        }
        *gz=(int64_t)round((vz-grid->z.o)/grid->z.d);
    }
    else {
        x=coords[0];
        z=coords[2];
        if (ndim==2 || mod_ndim==3) {
            y=coords[1];
        }
        switch (axis_indep) {
        case 1: // z
            v=z; vx=x; vy=y;
            break;
        case 2: // x
            v=x; vx=z; vy=y;
            break;
        case 3: // y
            v=y; vx=z; vy=x;
            break;
        default:
            ERROR(("Invalid value for axis_indep=%d",axis_indep));
            break;
        }
    }
    *gx=(int64_t)round((vx-grid->x.o)/grid->x.d);
    *gy=(int64_t)round((vy-grid->y.o)/grid->y.d);

    return v;
}

int rsf_get_grid_dim(const grid3d_t *grid) {
    return 3-((grid->x.n>1 || !FISZERO(grid->x.o))?0:1)
        -((grid->y.n>1 || !FISZERO(grid->y.o))?0:1)
        -((grid->z.n>1 || !FISZERO(grid->z.o))?0:1);
}

void rsf_from_hrz(const hrz_model_t* mod, int idx, 
                      const grid3d_t *grid, rsf3d_t *rsf, 
                      int axis_indep, int gathers, float no_value)
{
    int i;
    const hrz_surface_t *s = (hrz_surface_t *)array_get_at(mod->surfaces, idx);
    const int ndim=rsf_get_grid_dim(grid);
    int64_t points_outside_boundary, total_points;

    rsf_set3d_defmem_from_grid(grid, NULL, rsf);

    {
        size_t j;
        for(j=0; j < ((size_t)grid->x.n*grid->y.n)*grid->z.n; ++j) {
            *((float*)rsf->prvt+j) = no_value;
        }
    }

    // check validity of independent axis
    if (axis_indep < 1 || axis_indep > mod->ndim || axis_indep > ndim+1) {
        ERROR(("Invalid value for axis_indep=%d",axis_indep));
    }

    INFOV((10,"axis_indep=%d, ndim=%d",axis_indep,ndim));
    INFOV((10,"grid.x: (%d,%g,%g) grid.y: (%d,%g,%g) grid.z: (%d,%g,%g)",
           grid->x.n,grid->x.d,grid->x.o,
           grid->y.n,grid->y.d,grid->y.o,
           grid->z.n,grid->z.d,grid->z.o));

    total_points = points_outside_boundary = 0;
    // for each line
    for (i = 0; i < array_size(s->lines); ++i) {
        int j;
        const hrz_line_t* l = (hrz_line_t *)array_get_at(s->lines, i);

        // for each point
        for (j = 0; j < l->npoints; ++j) {
            int64_t gx, gy, gz; //x,y,z from grid
            double v;

            v=rsf_get_grid_indices(mod->ndim, axis_indep, gathers, l->coords + j*mod->ndim,
                                       grid, &gz, &gx, &gy);

            // check the indices
            if (gx < 0 || gx >= grid->x.n || gy < 0 || gy >= grid->y.n || gz<0 || gz >= grid->z.n) {
                INFOV((10, "Horizon \"%s\" line %d point %d outside grid boundary.",mod->name,i,j));
                points_outside_boundary += 1;
            }
            else {
                INFOV((10,"inserted v=%g at (%lld,%lld,%lld)",(float)v,gx,gy,gz));
                *(rsf->v(rsf, gx, gy, gz))=(float)v; // insert into grid
            }

            total_points += 1;
        }
    }

    if(points_outside_boundary > 0)
        WARN(("Horizon \"%s\" has %Ld (%.3f%%) point(s) outside grid boundary.",
              mod->name, (long long int)points_outside_boundary,
              (double)points_outside_boundary*100.0/(double)total_points));
}

void hrz_layer_from_first_hz(const hrz_model_t *top, 
                                 const hrz_model_t *bottom, 
                                 const grid3d_t* grid, 
                                 hrz_layer_t *layer)
{
    float no_value=hrz_no_value_from_axa(&grid->z);

    if (top != NULL) {
        layer->rsftop=alloc1type(rsf3d_t, 1);
        rsf_from_hrz(top, 0, grid, layer->rsftop, 1, 0, no_value);
        sparse_linear_interpolator(layer->rsftop, no_value, 1.0e-5f, 2.0, 0, -1.0);
    }

    if (bottom != NULL) {
        layer->rsfbottom=alloc1type(rsf3d_t, 1);
        rsf_from_hrz(bottom, 0, grid, layer->rsfbottom, 1, 0, no_value);
        sparse_linear_interpolator(layer->rsfbottom, no_value, 1.0e-5f, 2.0, 0, -1.0);
    }
}

hrz_layer_t *hrz_single_constraining_layer(const char *topfile, 
                                                   const char *bottomfile, 
                                                   double ztop, double zbottom, 
                                                   const grid3d_t* grid)
{
    hrz_layer_t *layer = alloc1type(hrz_layer_t, 1);
    hrz_model_t *top=NULL, *bottom=NULL;

    if (ztop > zbottom) {
        ERROR(("top=%g, bottom=%g constraints incompatible",ztop,zbottom));
    }

    layer->ztop = ztop;
    layer->zbottom = zbottom;
    layer->rsftop=NULL;
    layer->rsfbottom=NULL;

    // process top horizon
    if (topfile && topfile[0] != 0) {
        top=parse_pick_file_print(topfile);

        if (top->ndim != 3) {
            ERROR(("Horizon .pick files must have 3 columns; %s has %d.",topfile,top->ndim));
        }
    }

    // process bottom horizon
    if (bottomfile && bottomfile[0] != 0) {
        bottom = parse_pick_file_print(bottomfile);

        if (bottom->ndim != 3) {
            ERROR(("Horizon .pick files must have 3 columns; %s has %d.",bottomfile,bottom->ndim));
        }
    }

    /* interpolate */

    if (top != NULL || bottom != NULL) {
        INFOV((1,"Interpolating horizons"));
        hrz_layer_from_first_hz(top, bottom, grid, layer);
    }

    if (top != NULL) {
        clear_hrz_model(top);
        free(top);
    }

    if (bottom != NULL) {
        clear_hrz_model(bottom);
        free(bottom);
    }

    return layer;
}

void hrz_layer_get_limits(const hrz_layer_t* layers, int N, double x, double y, double *top, double *bottom) {

    int i;

    (*top) = INFINITY;
    (*bottom) = -INFINITY;

    for(i=0;i<N;i++) {
        double ztop=layers[i].ztop, zbottom=layers[i].zbottom;
        rsf3d_t *rsftop=layers[i].rsftop, *rsfbottom=layers[i].rsfbottom;

        if (rsftop != NULL) {
            ztop=maxr(se_rsf3d_value_binterp(x, y, rsftop->g.z.o, rsftop),ztop);
        }

        if (rsfbottom != NULL) {
            zbottom=minr(se_rsf3d_value_binterp(x, y, rsfbottom->g.z.o, rsfbottom),zbottom);
        }

        (*top) = minr(*top, ztop);
        (*bottom) = maxr(*bottom, zbottom);

    }


}

int hrz_layer_check_pt(const hrz_layer_t* layers, int N, double z, double x, double y) {
    int i;

    for(i=0;i<N;i++) {
        double ztop=layers[i].ztop, zbottom=layers[i].zbottom;
        rsf3d_t *rsftop=layers[i].rsftop, *rsfbottom=layers[i].rsfbottom;

        if (rsftop != NULL) {
            ztop=maxr(se_rsf3d_value_binterp(x, y, rsftop->g.z.o, rsftop),ztop);
        }

        if (rsfbottom != NULL) {
            zbottom=minr(se_rsf3d_value_binterp(x, y, rsfbottom->g.z.o, rsfbottom),zbottom);
        }

        if (z >= ztop && z <= zbottom) {
            return 1;
        }
    }

    return 0;
}

void clear_hrz_layer(hrz_layer_t *l) {
    if (l->rsftop != NULL) {
        rsf_destroy3d_and_data_defmem(l->rsftop);
    }
    if (l->rsfbottom != NULL) {
        rsf_destroy3d_and_data_defmem(l->rsfbottom);
    }
}

void init_hrz_line(hrz_line_t* l)
{
    l->name = strdup("<no-name>");
    l->npoints = 0;
    l->coords = NULL;
}

void clear_hrz_line(hrz_line_t* l)
{
    if (l != NULL) {
        free(l->name);
        free(l->coords);
    }
}

void init_hrz_surface(hrz_surface_t* s)
{
    s->name = strdup("<no-name>");
    s->attrs = create_hash();
    s->lines = create_array();
}

void clear_hrz_surface(hrz_surface_t* s)
{
    if (s != NULL) {
        free(s->name);
        destroy_hash_and_entries(s->attrs, 1);
        array_process(s->lines, (process_value_fn) clear_hrz_line);
        destroy_array(s->lines, 1);
    }
}

void init_hrz_model(hrz_model_t* m)
{
    m->name = strdup("<no-name>");
    m->attrs = create_hash();
    m->ndim = 0;
    m->surfaces = create_array();
}

void clear_hrz_model(hrz_model_t* m)
{
    if (m != NULL) {
        free(m->name);
        destroy_hash_and_entries(m->attrs, 1);
        array_process(m->surfaces, (process_value_fn) clear_hrz_surface);
        destroy_array(m->surfaces, 1);
    }
}

float hrz_no_value_from_axa(const axa_t *ax)
{
    real z_min=axa_min_x(ax), z_range=axa_max_x(ax)-z_min;
    return (float)( z_min - 0.1*(z_range+1.0) );
}

const char* hrz_get_attr(se_hash attrs, const char* key)
{
    return (const char*) ht_get(attrs, key);
}

int hrz_get_attr_i(se_hash attrs, const char* key)
{
    const char* sval = (const char*) ht_get(attrs, key);

    int ival;
    if (parse_int(sval, &ival) != 0) {
        WARN(("Cannot parse horizon model attribute '%s=%s': not an integer", key, sval));
    }
    return ival;
}

double hrz_get_attr_d(se_hash attrs, const char* key)
{
    const char* sval = (const char*) ht_get(attrs, key);

    double dval;
    if (parse_double(sval, &dval) != 0) {
        WARN(("Cannot parse horizon model attribute '%s=%s': not a number", key, sval));
    }
    return dval;
}
