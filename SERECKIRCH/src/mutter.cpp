
#include "../include/mutter.h"

static float dt, t0;
static int nt;
static bool abs0, inner, hyper;              

void mutter_init (int n1 , float o1, float d1 /* time axis */, 
		  bool abs1                   /* if twosided */,
		  bool inner1                 /* inner mute */, 
		  bool hyper1                 /* hyperbolic mute */)
/*< Initialize >*/
{
    nt = n1;
    t0 = o1;
    dt = d1;
    abs0 = abs1;
    inner = inner1;
    hyper = hyper1;
}

void mutter (float tp     /* time step */, 
	     float slope0 /* first slope */, 
	     float slopep /* second slope */, 
	     float x      /* offset */, 
	     float *data  /* trace */,
	     bool  nan   /* nan instaed of zeros*/)
/*< Mute >*/
{
    int it;
    float wt, t;

    if (abs0) x = fabsf(x);

    for (it=0; it < nt; it++) {
	t = t0+it*dt;
	if (hyper) t *= t;
	wt = t - x * slope0;
	if ((inner && wt > 0.) || (!inner && wt < 0.)) {
	    if (nan)
#ifdef NAN
	        data[it]= NAN;
#else
	    	data[it]= 0.0 /0.0;
#endif
	    else
	    	data[it] = 0.;
	    	
	} else {
	    wt = t - tp - x * slopep;
	    if ((inner && wt >=0.) || (!inner && wt <= 0.)) {
		wt = sinf(0.5 * M_PI * 
			  (t-x*slope0)/(tp+x*(slopep-slope0)));
		data[it] *= (wt*wt);
	    } 
	}
    }
}

/* 	$Id$	 */
