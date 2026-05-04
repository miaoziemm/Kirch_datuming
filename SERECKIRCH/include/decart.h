#ifndef DECART_H
#define DECART_H
#include <SEBASIC/include/se_basic.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "pqueue.h"
/* if off_t, fseeko, and ftello not defined */
/* redefine these items back to the int, fseek, ftell */
#ifndef off_t
#define off_t long
#endif


void sf_line2cart(int dim         /* number of dimensions */, 
		  const int* nn /* box size [dim] */, 
		  int i         /* line coordinate */, 
		  int* ii       /* cartesian coordinates [dim] */);

int sf_cart2line(int dim         /* number of dimensions */, 
		 const int* nn /* box size [dim] */, 
		 const int* ii /* cartesian coordinates [dim] */) ;

int sf_first_index (int i          /* dimension [0...dim-1] */, 
		    int j        /* line coordinate */, 
		    int dim        /* number of dimensions */, 
		    const int *n /* box size [dim] */, 
		    const int *s /* step [dim] */);
void sf_large_line2cart(int dim         /* number of dimensions */, 
			const off_t* nn /* box size [dim] */, 
			off_t i         /* line coordinate */, 
			off_t* ii       /* cartesian coordinates [dim] */);
off_t sf_large_cart2line(int dim         /* number of dimensions */, 
			 const off_t* nn /* box size [dim] */, 
			 const off_t* ii /* cartesian coordinates [dim] */);

off_t sf_large_first_index (int i          /* dimension [0...dim-1] */, 
			    off_t j        /* line coordinate */, 
			    int dim        /* number of dimensions */, 
			    const off_t *n /* box size [dim] */, 
			    const off_t *s /* step [dim] */);

#endif