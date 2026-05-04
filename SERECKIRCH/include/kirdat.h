#ifndef KIRDAT_H
#define KIRDAT_H

#include <SEBASIC/include/se_basic.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

void filt_init(float dt0 /* time sampling */,
	       float length /* filter length */);

void filt_close(void);

void filt_set(float tau /* time delay */);

float kirdat_pick(float delta /* sample position */,
	   float* trace /* input trace */,
	   int shift /* sample shift */);


#endif