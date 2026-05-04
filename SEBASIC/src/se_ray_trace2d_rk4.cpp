#include "se_ray_trace_rk.h"
#include "se_alloc.h"
#include "se_su_fft.h"
#include "se_basic_math.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <complex>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef NINT
#define NINT(x) ((int)((x) > 0.0 ? (x) + 0.5 : (x)-0.5))
#endif

/* ======================================================================
 * Ray Tracing Functions
 * ====================================================================== */

/* File-local helpers, scoped via rt2d_ prefix to keep visibility contained. */
static inline void rt2d_eval_medium(float x, float z, se_vel_approx_t *vel_approx,
						  float &v, float &dvdx, float &dvdz,
						  float &ddvdxdx, float &ddvdxdz, float &ddvdzdz)
{
	v3d_t pt, gx;
	m3d_t hv;
	double v_val = 0.0;

	pt.v[0] = x;
	pt.v[1] = 0.0;
	pt.v[2] = z;

	if (!se_vel_get_vel_grad_hess_at_pt(vel_approx, &pt, NULL, &v_val, &gx, NULL, &hv))
	{
		v = 0.0f;
		dvdx = dvdz = 0.0f;
		ddvdxdx = ddvdxdz = ddvdzdz = 0.0f;
		return;
	}

	v = static_cast<float>(v_val);
	dvdx = static_cast<float>(gx.v[0]);
	dvdz = static_cast<float>(gx.v[2]);
	ddvdxdx = static_cast<float>(hv.m[0]);
	ddvdxdz = static_cast<float>(hv.m[2]);
	ddvdzdz = static_cast<float>(hv.m[8]);
}

static inline float rt2d_ddvdndn(float c, float s, float ddvdxdx, float ddvdxdz, float ddvdzdz)
{
	/* second derivative along normal direction */
	return ddvdxdx * c * c + ddvdzdz * s * s - 2.0f * ddvdxdz * s * c;
}

ray_path_rk_t *ray_trace_make_2d_rk4(ray_init_rk_t *init, grid2d_t *model_grid, se_vel_approx_t *vel_approx)
{
	int it, kmah;
	float t, x, z, a, c, s, p1, q1, p2, q2,
		lx, lz,
		v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz,
		vv, ddvdndn;
	ray_path_rk_t *ray;
	ray_step_rk_t *rs;

    /* Unpack parameters */
    float x0 = init->x0;
    float z0 = init->z0;
    float a0 = init->angle;
    int nt = init->nt;
    float dt = init->dt;
    float ft = init->ft;

    int nx = model_grid->x.n;
    float dx = model_grid->x.d;
    float fx = model_grid->x.o;
    int nz = model_grid->y.n;
    float dz = model_grid->y.d;
    float fz = model_grid->y.o;

	/* last x and z in velocity model */
	lx = fx + (nx - 1) * dx;
	lz = fz + (nz - 1) * dz;

	/* ensure takeoff point is within model */
	if (x0 < fx || x0 > lx || z0 < fz || z0 > lz)
		return NULL;

	/* allocate space for ray and raysteps */
	ray = (ray_path_rk_t *)alloc1(1, sizeof(ray_path_rk_t));
	rs = (ray_step_rk_t *)alloc1(nt, sizeof(ray_step_rk_t));

	/* cosine and sine of takeoff angle */
	c = cos(a0);
	s = sin(a0);

	/* velocity and derivatives at takeoff point */
	rt2d_eval_medium(x0, z0, vel_approx, v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz);

	ddvdndn = 0;
	vv = v * v;

	/* first ray step */
	rs[0].t = t = ft;
	rs[0].a = a = a0;
    rs[0].b = 0.0f; /* 2D assumption */
	rs[0].x = x = x0;
    rs[0].y = 0.0f; /* 2D assumption */
	rs[0].z = z = z0;
	rs[0].q1 = q1 = 1.0;
	rs[0].p1 = p1 = 0.0;
	rs[0].q2 = q2 = 0.0;
	rs[0].p2 = p2 = 1.0;
	rs[0].kmah = kmah = 0;
	rs[0].c = c;
	rs[0].s = s;
    rs[0].cb = 1.0f; /* cos(0) */
    rs[0].sb = 0.0f; /* sin(0) */
    rs[0].px = s / v;
    rs[0].py = 0.0f;
    rs[0].pz = c / v;
	rs[0].v = v;
	rs[0].dvdx = dvdx;
    rs[0].dvdy = 0.0f;
	rs[0].dvdz = dvdz;

	const float h = dt;
	const float hhalf = dt * 0.5f;
	const float hsixth = dt / 6.0f;

	/* loop over time steps */
	for (it = 1; it < nt; ++it)
	{
		/* variables used for Runge-Kutta integration */
		float q2old, xt, zt, at, p1t, q1t, p2t, q2t,
			dx_step, dz_step, da, dp1, dq1, dp2, dq2,
			dxt, dzt, dat, dp1t, dq1t, dp2t, dq2t,
			dxm, dzm, dam, dp1m, dq1m, dp2m, dq2m;

		/* if ray is out of bounds, break */
		if (x < fx || x > lx || z < fz || z > lz)
			break;

		/* remember old q2 */
		q2old = q2;

		/* step 1 of 4th-order Runge-Kutta */
		dx_step = v * s;
		dz_step = v * c;
		da = dvdz * s - dvdx * c;
		dp1 = -ddvdndn * q1 / v;
		dq1 = vv * p1;
		dp2 = -ddvdndn * q2 / v;
		dq2 = vv * p2;
		xt = x + hhalf * dx_step;
		zt = z + hhalf * dz_step;
		at = a + hhalf * da;
		p1t = p1 + hhalf * dp1;
		q1t = q1 + hhalf * dq1;
		p2t = p2 + hhalf * dp2;
		q2t = q2 + hhalf * dq2;

		rt2d_eval_medium(xt, zt, vel_approx, v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz);
        
		ddvdndn = rt2d_ddvdndn(c, s, ddvdxdx, ddvdxdz, ddvdzdz);
		vv = v * v;

		/* step 2 of 4th-order Runge-Kutta */
		dxt = v * sin(at);
		dzt = v * cos(at);
		dat = dvdz * sin(at) - dvdx * cos(at);
		dp1t = -ddvdndn * q1t / v;
		dq1t = vv * p1t;
		dp2t = -ddvdndn * q2t / v;
		dq2t = vv * p2t;
		xt = x + hhalf * dxt;
		zt = z + hhalf * dzt;
		at = a + hhalf * dat;
		p1t = p1 + hhalf * dp1t;
		q1t = q1 + hhalf * dq1t;
		p2t = p2 + hhalf * dp2t;
		q2t = q2 + hhalf * dq2t;

		rt2d_eval_medium(xt, zt, vel_approx, v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz);

		ddvdndn = rt2d_ddvdndn(cos(at), sin(at), ddvdxdx, ddvdxdz, ddvdzdz);
		vv = v * v;

		/* step 3 of 4th-order Runge-Kutta */
		dxm = v * sin(at);
		dzm = v * cos(at);
		dam = dvdz * sin(at) - dvdx * cos(at);
		dp1m = -ddvdndn * q1t / v;
		dq1m = vv * p1t;
		dp2m = -ddvdndn * q2t / v;
		dq2m = vv * p2t;
		xt = x + h * dxm;
		zt = z + h * dzm;
		at = a + h * dam;
		p1t = p1 + h * dp1m;
		q1t = q1 + h * dq1m;
		p2t = p2 + h * dp2m;
		q2t = q2 + h * dq2m;

		rt2d_eval_medium(xt, zt, vel_approx, v, dvdx, dvdz, ddvdxdx, ddvdxdz, ddvdzdz);

		ddvdndn = rt2d_ddvdndn(cos(at), sin(at), ddvdxdx, ddvdxdz, ddvdzdz);
		vv = v * v;

		/* step 4 of 4th-order Runge-Kutta */
		x += hsixth * (dx_step + 2.0 * (dxt + dxm) + v * sin(at));
		z += hsixth * (dz_step + 2.0 * (dzt + dzm) + v * cos(at));
		a += hsixth * (da + 2.0 * (dat + dam) + dvdz * sin(at) - dvdx * cos(at));
		p1 += hsixth * (dp1 + 2.0 * (dp1t + dp1m) - ddvdndn * q1t / v);
		q1 += hsixth * (dq1 + 2.0 * (dq1t + dq1m) + vv * p1t);
		p2 += hsixth * (dp2 + 2.0 * (dp2t + dp2m) - ddvdndn * q2t / v);
		q2 += hsixth * (dq2 + 2.0 * (dq2t + dq2m) + vv * p2t);

		/* update ray step */
		rs[it].t = t + h;
		t = rs[it].t;
		rs[it].a = a;
        rs[it].b = 0.0f;
		rs[it].x = x;
        rs[it].y = 0.0f;
		rs[it].z = z;
		rs[it].q1 = q1;
		rs[it].p1 = p1;
		rs[it].q2 = q2;
		rs[it].p2 = p2;
		rs[it].c = cos(a);
		rs[it].s = sin(a);
        rs[it].cb = 1.0f;
        rs[it].sb = 0.0f;
        rs[it].px = rs[it].s / v;
        rs[it].py = 0.0f;
        rs[it].pz = rs[it].c / v;
		rs[it].v = v;
		rs[it].dvdx = dvdx;
        rs[it].dvdy = 0.0f;
		rs[it].dvdz = dvdz;

		/* check for KMAH index change */
		if (q2 * q2old < 0.0)
			kmah++;
		rs[it].kmah = kmah;
	}

	/* set number of ray steps */
	ray->nrs = it;
	ray->rs = rs;
	ray->nc = 0;
	ray->c = NULL;
	ray->ic = 0;

	return ray;
}

void ray_trace_free_2d_rk4(ray_path_rk_t *ray)
{
	if (ray->c != NULL)
		free1((void *)ray->c);
	free1((void *)ray->rs);
	free1((void *)ray);
}

static ray_circle_rk_t *ray_trace_make_circles_2d_rk4(int nc, int nrs, ray_step_rk_t *rs)
{
	int nrsc, ic, irsf, irsl, irs;
	float xmin, xmax, zmin, zmax, x, z, r;
	ray_circle_rk_t *c;

	/* allocate space for circles */
	c = (ray_circle_rk_t *)alloc1(nc, sizeof(ray_circle_rk_t));

	/* determine typical number of ray steps per circle */
	nrsc = 1 + (nrs - 1) / nc;

	/* loop over circles */
	for (ic = 0; ic < nc; ++ic)
	{

		/* index of first and last raystep */
		irsf = ic * nrsc;
		irsl = irsf + nrsc - 1;
		if (irsf >= nrs)
			irsf = nrs - 1;
		if (irsl >= nrs)
			irsl = nrs - 1;

		/* coordinate bounds of ray steps */
		xmin = xmax = rs[irsf].x;
		zmin = zmax = rs[irsf].z;
		for (irs = irsf + 1; irs <= irsl; ++irs)
		{
			if (rs[irs].x < xmin)
				xmin = rs[irs].x;
			if (rs[irs].x > xmax)
				xmax = rs[irs].x;
			if (rs[irs].z < zmin)
				zmin = rs[irs].z;
			if (rs[irs].z > zmax)
				zmax = rs[irs].z;
		}

		/* center and radius of circle */
		x = 0.5 * (xmin + xmax);
		z = 0.5 * (zmin + zmax);
		r = sqrt((x - xmin) * (x - xmin) + (z - zmin) * (z - zmin));

		/* set circle */
		c[ic].irsf = irsf;
		c[ic].irsl = irsl;
		c[ic].x = x;
        c[ic].y = 0.0f;
		c[ic].z = z;
		c[ic].r = r;
	}

	return c;
}

int ray_trace_nearest_step_2d_rk4(ray_path_rk_t *ray, float x, float z)
{
	int nrs = ray->nrs, ic = ray->ic, nc = ray->nc;
	ray_step_rk_t *rs = ray->rs;
	ray_circle_rk_t *c = (ray_circle_rk_t *)ray->c;
	int irs, irsf, irsl, irsmin = 0, update, jc;
	float dsmin, ds, dx, dz, dmin, rdmin;

	/* if necessary, make circles localizing ray steps */
	if (c == NULL)
	{
		ray->ic = ic = 0;
		ray->nc = nc = sqrt((float)nrs);
		ray->c = c = ray_trace_make_circles_2d_rk4(nc, nrs, rs);
	}

	/* determine nearest ray step */
	/* first check the circle containing the nearest step from previous call */
	irsf = c[ic].irsf;
	irsl = c[ic].irsl;
	dx = x - rs[irsf].x;
	dz = z - rs[irsf].z;
	dsmin = dx * dx + dz * dz;
	irsmin = irsf;
	for (irs = irsf + 1; irs <= irsl; ++irs)
	{
		dx = x - rs[irs].x;
		dz = z - rs[irs].z;
		ds = dx * dx + dz * dz;
		if (ds < dsmin)
		{
			dsmin = ds;
			irsmin = irs;
		}
	}

	/* if nearest step is on boundary of circle, check all circles */
	/* initialize dmin from dsmin so it's defined in all branches */
	dmin = sqrt(dsmin);
	if (irsmin == irsf || irsmin == irsl)
	{
		update = 1;
	}
	else
	{
		/* otherwise, check nearest circle to see if we can avoid checking all */
		rdmin = c[ic].r + dmin;
		dx = x - c[ic].x;
		dz = z - c[ic].z;
		if (dx * dx + dz * dz > rdmin * rdmin)
			update = 1;
		else
			update = 0;
	}

	/* if necessary, search all circles */
	if (update)
	{
		for (jc = 0; jc < nc; ++jc)
		{
			dx = x - c[jc].x;
			dz = z - c[jc].z;
			ds = dx * dx + dz * dz;
			rdmin = c[jc].r + dmin;
			if (ds < rdmin * rdmin)
			{
				irsf = c[jc].irsf;
				irsl = c[jc].irsl;
				for (irs = irsf; irs <= irsl; ++irs)
				{
					dx = x - rs[irs].x;
					dz = z - rs[irs].z;
					ds = dx * dx + dz * dz;
					if (ds < dsmin)
					{
						dsmin = ds;
						irsmin = irs;
						ic = jc;
						dmin = sqrt(dsmin);
					}
				}
			}
		}
		ray->ic = ic;
	}

	return irsmin;
}
