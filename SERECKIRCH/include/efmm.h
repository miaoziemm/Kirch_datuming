#ifndef EFMM_H
#define EFMM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

typedef struct efmm_s
{
    int nx_work;
    int nz_work;
    float dx;
    float dz;

    float src_x;
    float src_z;

    float *vel;
    float *tt;
} efmm_t;


int efmm_init(efmm_t *efmm, int nx_work, int nz_work, float dx, float dz, float src_x, float src_z);

int efmm_set_src(efmm_t *efmm, float src_x, float src_z);
int efmm_set_vel(efmm_t *efmm, float *vel);

int efmm_free(efmm_t *efmm);

int efmm_solver(efmm_t *efmm);

void efmm_get_bottom_tt(efmm_t *efmm, float *bottom_tt);
void efmm_get_top_tt(efmm_t *efmm, float *top_tt);

void efmm_eikods_init (int n3,int n2,int n1);
/*< Initialize data dimensions >*/


void efmm_eikods (float* time                /* time */, 
         float* v                   /* slowness squared */, 
         int* in                    /* in/front/out flag */, 
         bool* plane                /* if plane source */, 
         int   n3,  int n2,  int n1 /* dimensions */, 
         float o3,float o2,float o1 /* origin */, 
         float d3,float d2,float d1 /* sampling */, 
         float s3,float s2,float s1 /* source */, 
         int   b3,  int b2,  int b1 /* box around the source */, 
         int order                  /* accuracy order (1,2,3) */, 
         int l                      /* direction of source perturbation */, 
         float* dl1, float* ds1     /* first-order derivatives */, 
         float* dl2, float* ds2     /* second-order derivatives */);
/*< Run fast marching eikonal solver >*/


#endif

