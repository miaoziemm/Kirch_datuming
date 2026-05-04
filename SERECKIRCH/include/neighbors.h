#ifndef NEIGHBORS_H
#define NEIGHBORS_H
#include <SEBASIC/include/se_basic.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "pqueue.h"
struct Upd {
    double stencil, value;
    double delta;
};

void sf_neighbors_init (int *in1     /* status flag [n[0]*n[1]*n[2]] */, 
            float *rdx1  /* grid sampling [3] */, 
            int *n1      /* grid samples [3] */, 
            int order1   /* accuracy order */, 
            float *time1 /* traveltime [n[0]*n[1]*n[2]] */);

int  sf_neighbours(int i) ;

int  sf_neighbours2(int i) ;
int sf_neighbors_distance(int np         /* number of points */,
			  float *vv1     /* slowness squared */,
			  float **points /* point coordinates[np][3] */,
			  float *d       /* grid sampling [3] */,
			  float *o       /* grid origin [3] */);
int sf_neighbors_nearsource(float* xs   /* source location [3] */, 
			    int* b      /* constant-velocity box around it [3] */, 
			    float* d    /* grid sampling [3] */, 
			    float* vv1  /* slowness [n[0]*n[1]*n[2]] */, 
			    bool *plane /* if plane-wave source */);
int sf_neighbors_nearsource_rtp(float* xs   /* source location [3] */, 
			    int* b      /* constant-velocity box around it [3] */, 
			    float* d    /* grid sampling [3] */, 
			    float* vv1  /* slowness [n[0]*n[1]*n[2]] */, 
			    bool *plane /* if plane-wave source */);
int sf_neighbors_surface(float* vv1  /* slowness [n[0]*n[1]*n[2]] */,
			 float* tt0  /* surface traveltime [n[1]*n[2]] */,
			 bool forw /* forward or backward continuation */);
int sf_neighbors_mask(float* vv1  /* slowness [n[0]*n[1]*n[2]] */,
		      float* tref /* reference traveltime [n[0]*n[1]*n[2]] */,
		      bool* known /* where known [n[0]*n[1]*n[2]] */,
		      bool forw   /* forward or backward continuation */);
              



#endif