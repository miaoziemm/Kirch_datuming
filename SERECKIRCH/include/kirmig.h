#ifndef KIRMIG_H
#define KIRMIG_H

#include <SEBASIC/include/se_basic.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


void doubint(int nt, float *trace);

void kirmig_init(int n1, float d1, float o1);

void kirmig_pick(bool adj, float ti, float deltat, float *sample, float *trace);



#endif