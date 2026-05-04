#include "../include/se_ray.h"


double ray_tbl_get_average_slowness(const ray_t* ray_tbl, 
                                        size_t tbl_min_i, 
                                        size_t tbl_max_i)
{
    double l, t;
    double s;
    
    l = ray_tbl_get_length(ray_tbl, tbl_min_i, tbl_max_i);
    t = ray_tbl_get_time(ray_tbl, tbl_min_i, tbl_max_i);

    if (!FISZERO(l)) {
        s = t / l;
    } else {
        s = v3d_norm(&ray_tbl[tbl_min_i].pqr);
    }

    return s;
}

double ray_tbl_get_time(const ray_t* ray_tbl, 
                            size_t tbl_min_i, 
                            size_t tbl_max_i)
{
    return (ray_tbl[tbl_max_i].t - ray_tbl[tbl_min_i].t);
}

double ray_tbl_get_length(const ray_t* ray_tbl, 
                              size_t tbl_min_i, 
                              size_t tbl_max_i)
{
    size_t i;
    double len;

    len = 0;
    for (i=tbl_min_i+1; i <= tbl_max_i; i++) {
        v3d_t diff;
        v3d_subtract(&diff, &ray_tbl[i].xyz, &ray_tbl[i-1].xyz);
        len += v3d_norm(&diff);
    }

    return len;
}

double ray_tbl_get_max_vel(const ray_t* ray_tbl, 
                               size_t tbl_min_i, 
                               size_t tbl_max_i)
{
    size_t i;
    double vel;

    vel = 0;
    for (i=tbl_min_i; i <= tbl_max_i; i++) {
        double s;
        s = v3d_norm(&ray_tbl[i].pqr);
        if (!FISZERO(s)) {
            vel = maxd(vel, 1/s);
        }

    }

    return vel;
}

double ray_get_offset(const ray_t* ray_A,
                          const ray_t* ray_B)
{
    v3d_t diff;
    
    v3d_subtract(&diff, 
                     &ray_A->xyz, 
                     &ray_B->xyz);

    return v3d_norm(&diff);

}

double ray_tbl_get_slowness(const ray_t* ray_tbl, 
                                size_t tbl_i)
{
    return v3d_norm(&ray_tbl[tbl_i].pqr);
}

double ray_tbl_get_z(const ray_t* ray_tbl, 
                         size_t tbl_i)
{
    return ray_tbl[tbl_i].xyz.v[2];
}

void ray_tbl_get_min_max_velocity(const ray_t* ray_tbl, 
                                      size_t tbl_min_i, 
                                      size_t tbl_max_i,
                                      double *min_vel, double *max_vel)
{
    size_t i;

    (*min_vel) = INFINITY;
    (*max_vel) = 0;

    for (i=tbl_min_i+1; i <= tbl_max_i; i++) {
        double v;
        v = 1/v3d_norm(&ray_tbl[i].pqr);
        if ( v < (*min_vel) ) (*min_vel) = v;
        if ( v > (*max_vel) ) (*max_vel) = v;
    }

}


double ray_tbl_get_cos_to_v3d(const ray_t* ray_tbl, 
                                  size_t tbl_i,
                                  const v3d_t *v)
{
    return v3d_cos(&ray_tbl[tbl_i].pqr,v);
}
