#ifndef EIKODS_H
#define EIKODS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "pqueue.h"
#include "neighbors.h"
#include "decart.h"
void eikods_init (int n3,int n2,int n1);
/*< Initialize data dimensions >*/


void eikods (float* time                /* time */, 
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

void eikods_close (void);
/*< Free allocated storage >*/


#endif