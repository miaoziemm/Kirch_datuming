#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


int main(int argc, char* argv[])
{
	se_par_init(argc, argv);
    bool verb;
    int it, nt, ih, nh, is, ns, nsg, nrg, left, right, ic, aper, shift, c, cc, hh;
    int ir, nr, jump, sleft, sright, tap;
    float sdatum, rdatum, length, t0, dt, h0, dh, s0, ds, sg0, dsg, rg0, drg, dist, tau, delta;
    float r, dr, s, h, coef;
    float ***tr_in, ***tr_out, **stable, **rtable;
	sep_t *in, *out, *sgreen, *rgreen, *interm;
	char *in_f = NULL, *out_f = NULL, *sgreen_f = NULL, *rgreen_f = NULL, *interm_f = NULL;

	if(!se_have_par("input_file")) ERROR(("Need input_file=")); else in_f = se_get_par_str("input_file");
	if(!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");
	if(!se_have_par("sgreen_file")) ERROR(("Need sgreen_file=")); else sgreen_f = se_get_par_str("sgreen_file");
	if(!se_have_par("rgreen_file")) ERROR(("Need rgreen_file=")); else rgreen_f = se_get_par_str("rgreen_file");
	

    in = sf_input("in");
    out = sf_output("out");



    
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getfloat("sdatum",&sdatum)) sf_error("Need sdatum=");
    /* source datum depth */

    if (!sf_getfloat("rdatum",&rdatum)) sf_error("Need rdatum=");
    /* receiver datum depth */

    if (!sf_getint("aperture",&aper)) aper=50;
    /* aperture (number of traces) */

    if (!sf_getint("taper",&tap)) tap=10;
    /* taper (number of traces) */

    if (!sf_getfloat("length",&length)) length=0.025;
    /* filter length (in seconds) */

    /* read input */
    if (!sf_histint(in,"n1",&nt)) sf_error("No nt=");
    if (!sf_histint(in,"n2",&nh)) sf_error("No nh=");
    if (!sf_histint(in,"n3",&ns)) sf_error("No ns=");

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No t0=");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No dt=");

    if (!sf_histfloat(in,"o2",&h0)) sf_error("No h0=");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No dh=");

    if (!sf_histfloat(in,"o3",&s0)) sf_error("No s0=");
    if (!sf_histfloat(in,"d3",&ds)) sf_error("No ds=");

    tr_in = sf_floatalloc3(nt,nh,ns);
    sf_floatread(tr_in[0][0],nt*nh*ns,in);

    /* allocate memory for output */
    tr_out = sf_floatalloc3(nt,nh,ns);

    /* read Green's function (source) */
    sgreen = sf_input("sgreen");

    if (!sf_histint(sgreen,"n1",&nsg)) sf_error("No nsg=");
    if (!sf_histfloat(sgreen,"o1",&sg0)) sf_error("No sg0=");
    if (!sf_histfloat(sgreen,"d1",&dsg)) sf_error("No dsg=");

    stable = sf_floatalloc2(nsg,nsg);
    sf_floatread(stable[0],nsg*nsg,sgreen);
    sf_fileclose(sgreen);

    /* read Green's function (receiver) */
    rgreen = sf_input("rgreen");

    if (!sf_histint(rgreen,"n1",&nrg)) sf_error("No nrg=");
    if (!sf_histfloat(rgreen,"o1",&rg0)) sf_error("No rg0=");
    if (!sf_histfloat(rgreen,"d1",&drg)) sf_error("No drg=");

    rtable = sf_floatalloc2(nrg,nrg);
    sf_floatread(rtable[0],nrg*nrg,rgreen);
    sf_fileclose(rgreen);

    /* output intermediate traces */
    if (NULL != sf_getstring("interm")) {
	interm = sf_output("interm");
    } else {
	interm = NULL;
    }

    /* initialize */
    filt_init(dt,length);

    /* common-shot gather */
#ifdef _OPENMP
#pragma omp parallel for private(ih,c,left,right,ic,cc,coef,tau,dist,shift,it,delta)
#endif
    for (is=0; is < ns; is++) {
	if (verb) sf_warning("Processing common-shot gather %d of %d.",is+1,ns);

	for (ih=0; ih < nh; ih++) {

	    c = (s0+is*ds+h0+ih*dh-rg0)/drg+0.5;
	    if (c < 0 || c > nrg-1) sf_error("Receiver table too small.");

	    /* aperture */
	    left  = (ih-aper < 0)?    0:    ih-aper;
	    right = (ih+aper > nh-1)? nh-1: ih+aper;
	    
	    for (ic=left; ic <= right; ic++) {
		
		cc = (s0+is*ds+h0+ic*dh-rg0)/drg+0.5;
		if (cc < 0 || cc > nrg-1) sf_error("Receiver table too small.");

		/* taper coefficient */
		coef = 1.;
		coef *= (ic-left  >= tap)? 1.: (ic-left)/tap;
		coef *= (right-ic >= tap)? 1.: (right-ic)/tap;

		/* time delay */
		tau = rtable[cc][c];

		/* distance */
		dist = rdatum*rdatum+(ic-ih)*dh*(ic-ih)*dh;

		/* filter (tau dependent) */
		filt_set(tau);

		shift = 0;
		delta = 0.;
		for (it=0; it < nt; it++) {
		    if (((float)it)*dt < tau) 
			continue;
		    else if (shift == 0)
			delta = (((float)it*dt)-tau)/dt;

		    tr_out[is][ih][it] += coef/SF_PI
			*dh*rdatum*tau/dist
			*pick(delta,tr_in[is][ic],shift);
		    shift++;
		}
	    }

	}
    }

    if (NULL != interm) sf_floatwrite(tr_out[0][0],nt*nh*ns,interm);

    /* zero input */
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < nt; it++) {
		tr_in[is][ih][it] = 0.;
	    }
	}
    }

    /* acquisition */
    s = fabsf((ns-1)*ds);
    h = fabsf((nh-1)*dh);

    if (fabsf(ds) >= fabsf(dh)) {
	dr = fabsf(dh);
	jump = 1;
    } else {
	dr = fabsf(ds);
	jump = dh/ds+0.5;
    }
    
    nr = (s+h)/dr+1.5;

    /* common-receiver gather */
#ifdef _OPENMP
#pragma omp parallel for private(r,sleft,sright,is,c,ih,left,right,ic,cc,hh,coef,tau,dist,shift,it,delta)
#endif
    for (ir=0; ir < nr; ir++) {
	if (verb) sf_warning("Processing common-receiver gather %d of %d.",ir+1,nr);

	r = ir*dr+((ds<=0.)?-1.:0.)*s+((dh<=0.)?-1.:0.)*h;
	
	/* source receiver reciprocity */
	sleft  = (ir*dr+((ds<=0.)?-1.:0.)*s+((ds<=0.)?0.:-1.)*h)/ds+0.5;
	sright = (ir*dr+((ds<=0.)?-1.:0.)*s+((ds<=0.)?-1.:0.)*h)/ds+0.5;

	if (sleft < 0) sleft = 0;
	if (sright > ns-1) sright = ns-1;
	
	/* in case of fabsf(ds)>=fabsf(dh) */
	left = (r-sleft*ds)/dh+0.5;
	if (left < 0 || left > nh-1) sleft++;

	right = (r-sright*ds)/dh+0.5;
	if (right < 0 || right > nh-1) sright--;

	for (is=sleft; is <= sright; is=is+jump) {
	    
	    c = (s0+is*ds-sg0)/dsg+0.5;
	    if (c < 0 || c > nsg-1) sf_error("Source table too small.");

	    ih = (r-is*ds)/dh+0.5;
	    
	    /* aperture */
	    left  = (is-jump*aper < sleft)?  sleft:  is-jump*aper;
	    right = (is+jump*aper > sright)? sright: is+jump*aper;
	    
	    for (ic=left; ic <= right; ic=ic+jump) {
		
		cc = (s0+ic*ds-sg0)/dsg+0.5;
		if (cc < 0 || cc > nsg-1) sf_error("Source table too small.");

		hh = (r-ic*ds)/dh+0.5;

		/* taper coefficient */
		coef = 1.;
		coef *= (ic-left  >= tap)? 1.: (ic-left)/jump/tap;
		coef *= (right-ic >= tap)? 1.: (right-ic)/jump/tap;
		
		/* time delay */
		tau = stable[cc][c];
		
		/* distance */
		dist = sdatum*sdatum+(ic-is)*ds*(ic-is)*ds;
		
		/* filter (tau dependent) */
		filt_set(tau);
		
		shift = 0;
		delta = 0.;
		for (it=0; it < nt; it++) {
		    if (((float)it)*dt < tau) 
			continue;
		    else if (shift == 0)
			delta = (((float)it*dt)-tau)/dt;
		    
		    tr_in[is][ih][it] += coef/SF_PI
			*ds*sdatum*tau/dist
			*pick(delta,tr_out[ic][hh],shift);
		    shift++;
		}
	    }
	}	
    }

    /* write output */
    sf_floatwrite(tr_in[0][0],nt*nh*ns,out);

    exit(0);
}
