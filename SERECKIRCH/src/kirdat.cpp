#include "../include/kirdat.h"

#ifdef SE_USE_OMP
#include <omp.h>
#endif

static float dt;
static int nsam, mts;
static float** filt;

void filt_init(float dt0 /* time sampling */,
	       float length /* filter length */)
/*< initialize filter >*/
{
    int mts;

    dt = dt0;
    nsam = (int)(length/dt)+2;

#ifdef SE_USE_OMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    filt = alloc2float(nsam,mts);
}

void filt_close(void)
/*< close >*/
{
    int its;

    for (its=0; its < mts; its++) free(filt[its]);

    free(filt);
}

void filt_set(float tau /* time delay */)
/*< set up filter >*/
{
    int its, isam;

#ifdef SE_USE_OMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    /* value */
    filt[its][0] = 0.;

    for (isam=1; isam < nsam; isam++) {
	filt[its][isam] = sqrtf(powf((tau+isam*dt)/tau,2.)-1.);
    }

    /* first derivative */
    for (isam=0; isam < nsam-1; isam++) {
	filt[its][isam] = filt[its][isam+1]-filt[its][isam];
    }

    /* second derivative */
    for (isam=nsam-2; isam > 0; isam--) {
	filt[its][isam] = filt[its][isam]-filt[its][isam-1];
    }
}

float kirdat_pick(float delta /* sample position */,
	   float* trace /* input trace */,
	   int shift /* sample shift */)
/*< filter input trace for one output sample >*/
{
    float value=0.;
    int its, isam;    

#ifdef SE_USE_OMP
    its = omp_get_thread_num();
#else
    its = 0;
#endif

    for (isam=0; (isam < nsam-1) && (shift-isam >= 0); isam++) {
	value += ((1.-delta)*trace[shift-isam]+delta*trace[shift-isam+1])
	    *filt[its][isam]/dt;
    }

    return value;
}
