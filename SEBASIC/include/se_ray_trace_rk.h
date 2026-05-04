#ifndef SE_RAY_TRACE_RK_H
#define SE_RAY_TRACE_RK_H

#include "se_grid.h"
#include "se_vel.h"

/* Ray types */

/* one step along ray */
typedef struct ray_step_rk_s
{
    float t;              /* time */
    float a;              /* angle (dip) */
    float b;              /* azimuth */
    float x, y, z;        /* x,y,z coordinates */
    float q1, p1, q2, p2; /* Cerveny's dynamic ray tracing solution */
    int kmah;             /* KMAH index */
    float c, s;           /* cos(angle) and sin(angle) */
    float cb, sb;         /* cos(azimuth) and sin(azimuth) */
    float px, py, pz;     /* slowness vector components */
    float v, dvdx, dvdy, dvdz;  /* velocity and its derivatives */
} ray_step_rk_t;

/* circle for efficiently finding nearest ray step */
typedef struct ray_circle_rk_s
{
    int irsf; /* index of first raystep in circle */
    int irsl; /* index of last raystep in circle */
    float x, y, z;  /* center of circle */
    float r;  /* radius of circle */
} ray_circle_rk_t;

/* one ray */
typedef struct ray_path_rk_s
{
    int nrs;        /* number of ray steps */
    ray_step_rk_t *rs; /* array[nrs] of ray steps */
    int nc;         /* number of circles */
    int ic;         /* index of circle containing nearest step */
    ray_circle_rk_t *c;      /* array[nc] of circles */
} ray_path_rk_t;

/* Ray initialization parameters */
typedef struct ray_init_rk_s {
    float x0, y0, z0;
    float angle;   /* dip */
    float azimuth; /* azimuth */
    int nt;
    float dt;
    float ft;
} ray_init_rk_t;

/* Function Prototypes */

/**
 * @brief Traces a ray through a 2D velocity model.
 */
ray_path_rk_t *ray_trace_make_2d_rk4(ray_init_rk_t *init, grid2d_t *model_grid, se_vel_approx_t *vel_approx);

/**
 * @brief Traces a ray through a 3D velocity model.
 */
ray_path_rk_t *ray_trace_make_3d_rk4(ray_init_rk_t *init, grid3d_t *model_grid, se_vel_approx_t *vel_approx);

/**
 * @brief Frees the memory allocated for a ray (2D version).
 */
void ray_trace_free_2d_rk4(ray_path_rk_t *ray);

/**
 * @brief Frees the memory allocated for a ray (3D version).
 */
void ray_trace_free_3d_rk4(ray_path_rk_t *ray);

/**
 * @brief Finds the index of the ray step nearest to a given point (x, z).
 */
int ray_trace_nearest_step_2d_rk4(ray_path_rk_t *ray, float x, float z);

/**
 * @brief Finds the index of the ray step nearest to a given point (x, y, z).
 */
int ray_trace_nearest_step_3d_rk4(ray_path_rk_t *ray, float x, float y, float z);

#endif
