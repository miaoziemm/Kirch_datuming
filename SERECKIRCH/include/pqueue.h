#ifndef PQUEUE_H
#define PQUEUE_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
enum {SF_IN, SF_FRONT, SF_OUT};

void sf_pqueue_init (int n);
void sf_pqueue_start (void);

void sf_pqueue_close (void);

void sf_pqueue_insert (float* v);

void sf_pqueue_insert2 (float* v);

float* sf_pqueue_extract (void);

float* sf_pqueue_extract2 (void);

void sf_pqueue_update (float **v);



#endif