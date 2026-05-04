#ifndef TPOW_H
#define TPOW_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <vector>
#include <cstddef>

/*
 * Apply time-power and space-power gain to seismic traces.
 *
 * data layout:
 *   data[((i2 * nx + ix) * nt) + it]
 *
 * nt   : number of time samples, corresponding to n1
 * nx   : number of traces in the second dimension, corresponding to n2
 * nblk : number of higher-dimensional blocks, corresponding to sf_leftsize(in,2)
 *
 * tpow : power of time, default in original Madagascar code is 2.0
 * xpow : power of space, default in original Madagascar code is 0.0
 *
 * dt, t0 : time sampling interval and origin
 * dx, x0 : space sampling interval and origin
 *
 * The formula follows the original code:
 *   tgain[it] = pow(t0 + (it + 1) * dt, tpow)
 *   xgain[ix] = pow(x0 + (ix + 1) * dx, xpow)
 */
void tpow_gain(
    float *data,
    int nt,
    int nx,
    int nblk,
    float tpow,
    float xpow,
    float dt,
    float t0,
    float dx,
    float x0
);


/*
 * Apply time-power gain to a single trace.
 *
 * trace layout:
 *   trace[it], it = 0, 1, ..., nt-1
 */
void tpow_gain_trace(
    float *trace,
    int nt,
    float tpow,
    float dt,
    float t0
);

#endif