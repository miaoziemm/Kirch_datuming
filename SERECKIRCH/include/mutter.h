#ifndef MUTTER_H
#define MUTTER_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

void mutter_init (int n1 , float o1, float d1 /* time axis */, 
          bool abs1                   /* if twosided */,
          bool inner1                 /* inner mute */, 
          bool hyper1                 /* hyperbolic mute*/);

void mutter (float tp     /* time step */,
            float slope0 /* first slope */, 
            float slopep /* second slope */, 
            float x      /* offset */, 
            float *data  /* trace */,
            bool  nan   /* nan instaed of zeros*/);

            

#endif