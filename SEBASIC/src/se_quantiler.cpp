#include "../include/se_quantiler.h"


real se_quantiler_update_array(se_quantiler_t* zqt, const real* f, size_t n)
{
    int i;
    for(i = 0; i < (int)n; ++i)
        se_quantiler_update(zqt, f[i]);

    return se_quantiler_estimate(zqt);
}

real se_quantiler_update_farray(se_quantiler_t* zqt, const float* f, size_t n)
{
    int i;
    for(i = 0; i < (int)n; ++i)
        se_quantiler_update(zqt, f[i]);

    return se_quantiler_estimate(zqt);
}


void _se_quantiler_init_one(se_quantiler_t* zqt, real f)
{
    /* If ignoring null samples, check for null. */
    if (zqt->ignore_null && FEQUAL(f, zqt->fnull) )
        return;

    if(zqt->isminmax) {
        zqt->q2 = f;
        zqt->initialized = 1;
        return;
    }

    /* If fewer than 5 (non-null) samples, may not complete initialization. */
    if(zqt->m0 < 0.0) {
      zqt->m0 = 0.0;
      zqt->q0 = f;
    } else if( FEQUAL(zqt->m1, 0.0) ) {
        zqt->m1 = 1.0;
        zqt->q1 = f;
    } else if( FEQUAL(zqt->m2, 0.0) ) {
        zqt->m2 = 2.0;
        zqt->q2 = f;
    } else if( FEQUAL(zqt->m3, 0.0) ) {
        zqt->m3 = 3.0;
        zqt->q3 = f;
    } else if( FEQUAL(zqt->m4, 0.0) ) {
        zqt->m4 = 4.0;
        zqt->q4 = f;
    }

    if( FEQUAL(zqt->m4, 0.0) )
        return;

    /* Initialize marker heights to five samples sorted. */
    {
        real y[] = {zqt->q0, zqt->q1, zqt->q2, zqt->q3, zqt->q4};
        int i;
        for( i = 1; i < 5; ++i ) {
            int j;
            for( j = i; j > 0 && y[j-1] > y[j]; --j ) {
                real ytemp = y[j-1];
                y[j-1] = y[j];
                y[j] = ytemp;
            }
        }

        zqt->q0 = y[0];
        zqt->q1 = y[1];
        zqt->q2 = y[2];
        zqt->q3 = y[3];
        zqt->q4 = y[4];
    }


    /* Initialize desired marker positions. */
    zqt->f1 = 2.0 * zqt->q;
    zqt->f2 = 4.0 * zqt->q;
    zqt->f3 = 2.0 + 2.0*zqt->q;

    /* Compute increments in desired marker positions. */
    zqt->d1 = zqt->q / 2.0;
    zqt->d2 = zqt->q;
    zqt->d3 = (1.0 + zqt->q)/2.0;

    zqt->initialized = 1;
}
