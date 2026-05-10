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
    int verb;
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
	

	in = sep_open(in_f, SEP_READ, 0);
	out = sep_open(out_f, SEP_WRITE, 0);
	
    
	if(!se_have_par("verb")) verb = 1; else verb = se_get_par_int("verb");
    /* verbosity flag */
	if(!se_have_par("sdatum")) ERROR(("Need sdatum=")); else sdatum = se_get_par_float("sdatum");
	if(!se_have_par("rdatum")) ERROR(("Need rdatum=")); else rdatum = se_get_par_float("rdatum");

	if(!se_have_par("aperture")) aper=50; else aper = se_get_par_int("aperture");
    /* aperture (number of traces) */
	if(!se_have_par("taper")) tap=10; else tap = se_get_par_int("taper");
    /* taper (number of traces) */
	
	if(!se_have_par("length")) length=0.025; else length = se_get_par_float("length");
    /* filter length (in seconds) */

    /* read input */
	if(in->headers->ndim < 3) ERROR(("Input must be 3D."));
	nt = in->headers->n[0];
	nh = in->headers->n[1];
	ns = in->headers->n[2];
	t0 = in->headers->o[0];
	dt = in->headers->d[0];
	h0 = in->headers->o[1];
	dh = in->headers->d[1];
	s0 = in->headers->o[2];
	ds = in->headers->d[2];

    tr_in = alloc3float(nt,nh,ns);
	se_fsio_read_float(in->data->io, tr_in[0][0], nt*nh*ns);

    /* allocate memory for output */
    tr_out = alloc3float(nt,nh,ns);

	sgreen = sep_open(sgreen_f, SEP_READ, 0);
	
    /* read Green's function (source) */
	if(sgreen->headers->ndim < 2) ERROR(("Source Green's function must be 2D."));
	nsg = sgreen->headers->n[0];
	sg0 = sgreen->headers->o[0];
	dsg = sgreen->headers->d[0];

    stable = alloc2float(nsg,nsg);
    se_fsio_read_float(sgreen->data->io, stable[0], nsg*nsg);
    sep_close(sgreen);

    /* read Green's function (receiver) */
    rgreen = sep_open(rgreen_f, SEP_READ, 0);
	if(rgreen->headers->ndim < 2) ERROR(("Receiver Green's function must be 2D."));
	nrg = rgreen->headers->n[0];
	rg0 = rgreen->headers->o[0];
	drg = rgreen->headers->d[0];


    rtable = alloc2float(nrg,nrg);
    se_fsio_read_float(rgreen->data->io, rtable[0], nrg*nrg);
    sep_close(rgreen);

    /* output intermediate traces */
	if(!se_have_par("interm")) interm_f = NULL; else interm_f = se_get_par_str("interm");
	if(interm_f != NULL) {
		interm = sep_open(interm_f, SEP_WRITE, 0);
	} 


    /* initialize */
    filt_init(dt,length);

    /* common-shot gather */
#ifdef _OPENMP
#pragma omp parallel for private(ih,c,left,right,ic,cc,coef,tau,dist,shift,it,delta)
#endif
    for (is=0; is < ns; is++) {
	if (verb) WARN(("Processing common-shot gather %d of %d.",is+1,ns));

	for (ih=0; ih < nh; ih++) {

	    c = (s0+is*ds+h0+ih*dh-rg0)/drg+0.5;
	    if (c < 0 || c > nrg-1) ERROR(("Receiver table too small."));

	    /* aperture */
	    left  = (ih-aper < 0)?    0:    ih-aper;
	    right = (ih+aper > nh-1)? nh-1: ih+aper;
	    
	    for (ic=left; ic <= right; ic++) {
		
		cc = (s0+is*ds+h0+ic*dh-rg0)/drg+0.5;
		if (cc < 0 || cc > nrg-1) ERROR(("Receiver table too small."));

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

		    tr_out[is][ih][it] += coef/M_PI
			*dh*rdatum*tau/dist
			*kirdat_pick(delta,tr_in[is][ic],shift);
		    shift++;
		}
	    }

	}
    }

    if (NULL != interm_f) se_fsio_write_float(interm->data->io, tr_out[0][0], nt*nh*ns);

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
	if (verb) WARN(("Processing common-receiver gather %d of %d.",ir+1,nr));

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
	    if (c < 0 || c > nsg-1) ERROR(("Source table too small."));

	    ih = (r-is*ds)/dh+0.5;
	    
	    /* aperture */
	    left  = (is-jump*aper < sleft)?  sleft:  is-jump*aper;
	    right = (is+jump*aper > sright)? sright: is+jump*aper;
	    
	    for (ic=left; ic <= right; ic=ic+jump) {
		
		cc = (s0+ic*ds-sg0)/dsg+0.5;
		if (cc < 0 || cc > nsg-1) ERROR(("Source table too small."));

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
		    
		    tr_in[is][ih][it] += coef/M_PI
			*ds*sdatum*tau/dist
			*kirdat_pick(delta,tr_out[ic][hh],shift);
		    shift++;
		}
	    }
	}	
    }

    /* write output */
	se_fsio_write_float(out->data->io, tr_in[0][0], nt*nh*ns);

    sep_close(in);
	sep_close(out);

	free3float(tr_in);
	free3float(tr_out);
	free2float(stable);
	free2float(rtable);

	return 0;
}
