// 3D Gaussian beam ray tracing using 4th order Runge-Kutta, patterned after the
// existing 2D implementation.

#include "se_ray_trace_rk.h"
#include "se_alloc.h"
#include "se_basic_math.h"
#include <math.h>
#include <stdlib.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef NINT
#define NINT(x) ((int)((x) > 0.0 ? (x) + 0.5 : (x)-0.5))
#endif

struct rt3d_ray_state
{
	float x, y, z;
	float px, py, pz;
	float p1, q1, p2, q2;
};

struct rt3d_medium_eval
{
	float v;
	float dvdx, dvdy, dvdz;
	float ddvdn1, ddvdn2;
};

static inline float rt3d_project_hessian(const m3d_t &hv, float nx, float ny, float nz)
{
	return static_cast<float>(hv.m[0] * nx * nx + hv.m[4] * ny * ny + hv.m[8] * nz * nz +
								2.0 * (hv.m[1] * nx * ny + hv.m[2] * nx * nz + hv.m[5] * ny * nz));
}

static inline void rt3d_build_normals_from_slowness(float px, float py, float pz,
									   float &n1x, float &n1y, float &n1z,
									   float &n2x, float &n2y, float &n2z)
{
	v3d_t t;
	t.v[0] = px;
	t.v[1] = py;
	t.v[2] = pz;

	const float nrm = sqrtf(px * px + py * py + pz * pz);
	if (nrm < 1e-8f)
	{
		n1x = 1.0f;
		n1y = n1z = 0.0f;
		n2x = 0.0f;
		n2y = 1.0f;
		n2z = 0.0f;
		return;
	}

	v3d_t t_unit;
	v3d_assign_scaled(&t_unit, 1.0 / nrm, &t);

	v3d_t u, w;
	v3d_get_orthogonal_pair(&u, &w, &t_unit);

	n1x = static_cast<float>(u.v[0]);
	n1y = static_cast<float>(u.v[1]);
	n1z = static_cast<float>(u.v[2]);
	n2x = static_cast<float>(w.v[0]);
	n2y = static_cast<float>(w.v[1]);
	n2z = static_cast<float>(w.v[2]);
}

static inline void rt3d_evaluate_medium(const rt3d_ray_state &state, se_vel_approx_t *vel_approx, rt3d_medium_eval &me)
{
	v3d_t pt;
	pt.v[0] = state.x;
	pt.v[1] = state.y;
	pt.v[2] = state.z;

	v3d_t gx;
	m3d_t hv;
	double v_val = 0.0;
	if (!se_vel_get_vel_grad_hess_at_pt(vel_approx, &pt, NULL, &v_val, &gx, NULL, &hv))
	{
		me.v = 0.0f;
		me.dvdx = me.dvdy = me.dvdz = 0.0f;
		me.ddvdn1 = me.ddvdn2 = 0.0f;
		return;
	}

	me.v = static_cast<float>(v_val);
	me.dvdx = static_cast<float>(gx.v[0]);
	me.dvdy = static_cast<float>(gx.v[1]);
	me.dvdz = static_cast<float>(gx.v[2]);

	float n1x, n1y, n1z, n2x, n2y, n2z;
	rt3d_build_normals_from_slowness(state.px, state.py, state.pz, n1x, n1y, n1z, n2x, n2y, n2z);
	me.ddvdn1 = rt3d_project_hessian(hv, n1x, n1y, n1z);
	me.ddvdn2 = rt3d_project_hessian(hv, n2x, n2y, n2z);
}

static inline void rt3d_compute_derivs(const rt3d_ray_state &state, const rt3d_medium_eval &me, rt3d_ray_state &dst)
{
	const float vv = me.v * me.v;
	dst.x = vv * state.px;
	dst.y = vv * state.py;
	dst.z = vv * state.pz;
	dst.px = -me.dvdx / me.v;
	dst.py = -me.dvdy / me.v;
	dst.pz = -me.dvdz / me.v;
	dst.p1 = -me.ddvdn1 * state.q1 / me.v;
	dst.q1 = vv * state.p1;
	dst.p2 = -me.ddvdn2 * state.q2 / me.v;
	dst.q2 = vv * state.p2;
}

static inline rt3d_ray_state rt3d_rs_add(const rt3d_ray_state &lhs, const rt3d_ray_state &rhs)
{
	rt3d_ray_state out = lhs;
	out.x += rhs.x;
	out.y += rhs.y;
	out.z += rhs.z;
	out.px += rhs.px;
	out.py += rhs.py;
	out.pz += rhs.pz;
	out.p1 += rhs.p1;
	out.q1 += rhs.q1;
	out.p2 += rhs.p2;
	out.q2 += rhs.q2;
	return out;
}

static inline rt3d_ray_state rt3d_rs_scale(const rt3d_ray_state &rhs, float scale)
{
	rt3d_ray_state out = rhs;
	out.x *= scale;
	out.y *= scale;
	out.z *= scale;
	out.px *= scale;
	out.py *= scale;
	out.pz *= scale;
	out.p1 *= scale;
	out.q1 *= scale;
	out.p2 *= scale;
	out.q2 *= scale;
	return out;
}

static inline rt3d_ray_state rt3d_rs_add_scaled(const rt3d_ray_state &lhs, const rt3d_ray_state &rhs, float scale)
{
	return rt3d_rs_add(lhs, rt3d_rs_scale(rhs, scale));
}

static inline bool rt3d_inside_bounds(const rt3d_ray_state &s, float fx, float lx, float fy, float ly, float fz, float lz)
{
	return (s.x >= fx && s.x <= lx && s.y >= fy && s.y <= ly && s.z >= fz && s.z <= lz);
}

static ray_circle_rk_t *rt3d_make_circles_rk4(int nc, int nrs, ray_step_rk_t *rs)
{
	int nrsc, ic, irsf, irsl, irs;
	float xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, r;
	ray_circle_rk_t *c;

	c = (ray_circle_rk_t *)alloc1(nc, sizeof(ray_circle_rk_t));
	nrsc = 1 + (nrs - 1) / nc;

	for (ic = 0; ic < nc; ++ic)
	{
		irsf = ic * nrsc;
		irsl = irsf + nrsc - 1;
		if (irsf >= nrs)
			irsf = nrs - 1;
		if (irsl >= nrs)
			irsl = nrs - 1;

		xmin = xmax = rs[irsf].x;
		ymin = ymax = rs[irsf].y;
		zmin = zmax = rs[irsf].z;

		for (irs = irsf + 1; irs <= irsl; ++irs)
		{
			xmin = MIN(xmin, rs[irs].x);
			xmax = MAX(xmax, rs[irs].x);
			ymin = MIN(ymin, rs[irs].y);
			ymax = MAX(ymax, rs[irs].y);
			zmin = MIN(zmin, rs[irs].z);
			zmax = MAX(zmax, rs[irs].z);
		}

		x = 0.5f * (xmin + xmax);
		y = 0.5f * (ymin + ymax);
		z = 0.5f * (zmin + zmax);
		r = sqrtf((x - xmin) * (x - xmin) + (y - ymin) * (y - ymin) + (z - zmin) * (z - zmin));

		c[ic].irsf = irsf;
		c[ic].irsl = irsl;
		c[ic].x = x;
		c[ic].y = y;
		c[ic].z = z;
		c[ic].r = r;
	}

	return c;
}

ray_path_rk_t *ray_trace_make_3d_rk4(ray_init_rk_t *init, grid3d_t *model_grid, se_vel_approx_t *vel_approx)
{
	const float x0 = init->x0;
	const float y0 = init->y0;
	const float z0 = init->z0;
	const float az0 = init->azimuth; // azimuth
	const float dip0 = init->angle;  // dip from vertical
	const int nt = init->nt;
	const float dt = init->dt;
	const float ft = init->ft;

	const int nx = static_cast<int>(model_grid->x.n);
	const float dx = static_cast<float>(model_grid->x.d);
	const float fx = static_cast<float>(model_grid->x.o);
	const int ny = static_cast<int>(model_grid->y.n);
	const float dy = static_cast<float>(model_grid->y.d);
	const float fy = static_cast<float>(model_grid->y.o);
	const int nz = static_cast<int>(model_grid->z.n);
	const float dz = static_cast<float>(model_grid->z.d);
	const float fz = static_cast<float>(model_grid->z.o);

	const float lx = fx + (nx - 1) * dx;
	const float ly = fy + (ny - 1) * dy;
	const float lz = fz + (nz - 1) * dz;

	if (x0 < fx || x0 > lx || y0 < fy || y0 > ly || z0 < fz || z0 > lz)
	{
		return NULL;
	}

	ray_path_rk_t *ray = (ray_path_rk_t *)alloc1(1, sizeof(ray_path_rk_t));
	ray_step_rk_t *rs = (ray_step_rk_t *)alloc1(nt, sizeof(ray_step_rk_t));

	rt3d_ray_state state{};
	state.x = x0;
	state.y = y0;
	state.z = z0;
	state.p1 = 0.0f;
	state.q1 = 1.0f;
	state.p2 = 1.0f;
	state.q2 = 0.0f;

	const float sin_dip = sinf(dip0);
	const float cos_dip = cosf(dip0);
	const float cos_az = cosf(az0);
	const float sin_az = sinf(az0);

	rt3d_medium_eval me{};
	rt3d_evaluate_medium(state, vel_approx, me);

	if (me.v <= 0.0f)
	{
		free1((void *)rs);
		free1((void *)ray);
		return NULL;
	}

	state.px = sin_dip * cos_az / me.v;
	state.py = sin_dip * sin_az / me.v;
	state.pz = cos_dip / me.v;

	// Refresh medium evaluation so ddvdn terms use the correct slowness direction.
	rt3d_evaluate_medium(state, vel_approx, me);

	rs[0].t = ft;
	rs[0].a = dip0;
	rs[0].b = az0;
	rs[0].x = x0;
	rs[0].y = y0;
	rs[0].z = z0;
	rs[0].q1 = state.q1;
	rs[0].p1 = state.p1;
	rs[0].q2 = state.q2;
	rs[0].p2 = state.p2;
	rs[0].kmah = 0;
	rs[0].c = cos_dip;
	rs[0].s = sin_dip;
	rs[0].cb = cos_az;
	rs[0].sb = sin_az;
	rs[0].px = state.px;
	rs[0].py = state.py;
	rs[0].pz = state.pz;
	rs[0].v = me.v;
	rs[0].dvdx = me.dvdx;
	rs[0].dvdy = me.dvdy;
	rs[0].dvdz = me.dvdz;

	int kmah = 0;
	float t = ft;
	int written = 1;
	const float hhalf = 0.5f * dt;
	const float hsixth = dt / 6.0f;

	for (int it = 1; it < nt; ++it)
	{
		if (state.x < fx || state.x > lx || state.y < fy || state.y > ly || state.z < fz || state.z > lz)
		{
			break;
		}

		rt3d_ray_state k1{};
		rt3d_compute_derivs(state, me, k1);

		rt3d_ray_state state2 = rt3d_rs_add_scaled(state, k1, hhalf);
		if (!rt3d_inside_bounds(state2, fx, lx, fy, ly, fz, lz))
		{
			break;
		}
		rt3d_medium_eval me2{};
		rt3d_evaluate_medium(state2, vel_approx, me2);
		rt3d_ray_state k2{};
		rt3d_compute_derivs(state2, me2, k2);

		rt3d_ray_state state3 = rt3d_rs_add_scaled(state, k2, hhalf);
		if (!rt3d_inside_bounds(state3, fx, lx, fy, ly, fz, lz))
		{
			break;
		}
		rt3d_medium_eval me3{};
		rt3d_evaluate_medium(state3, vel_approx, me3);
		rt3d_ray_state k3{};
		rt3d_compute_derivs(state3, me3, k3);

		rt3d_ray_state state4 = rt3d_rs_add_scaled(state, k3, dt);
		if (!rt3d_inside_bounds(state4, fx, lx, fy, ly, fz, lz))
		{
			break;
		}
		rt3d_medium_eval me4{};
		rt3d_evaluate_medium(state4, vel_approx, me4);
		rt3d_ray_state k4{};
		rt3d_compute_derivs(state4, me4, k4);

		rt3d_ray_state incr = rt3d_rs_add(rt3d_rs_add(rt3d_rs_add(k1, rt3d_rs_scale(k2, 2.0f)), rt3d_rs_scale(k3, 2.0f)), k4);
		state = rt3d_rs_add_scaled(state, incr, hsixth);

		const float q2old = rs[it - 1].q2;
		t += dt;

		// Store medium properties at the updated position for output.
		rt3d_evaluate_medium(state, vel_approx, me);
		if (me.v <= 0.0f)
		{
			break;
		}

		const float pxv = state.px * me.v;
		const float pyv = state.py * me.v;
		const float pzv = state.pz * me.v;
		const float horiz = sqrtf(pxv * pxv + pyv * pyv);
		const float dip = atan2f(horiz, pzv);
		const float az = atan2f(pyv, pxv);

		rs[it].t = t;
		rs[it].a = dip;
		rs[it].b = az;
		rs[it].x = state.x;
		rs[it].y = state.y;
		rs[it].z = state.z;
		rs[it].q1 = state.q1;
		rs[it].p1 = state.p1;
		rs[it].q2 = state.q2;
		rs[it].p2 = state.p2;
		rs[it].c = cosf(dip);
		rs[it].s = sinf(dip);
		rs[it].cb = cosf(az);
		rs[it].sb = sinf(az);
		rs[it].px = state.px;
		rs[it].py = state.py;
		rs[it].pz = state.pz;
		rs[it].v = me.v;
		rs[it].dvdx = me.dvdx;
		rs[it].dvdy = me.dvdy;
		rs[it].dvdz = me.dvdz;

		if (state.q2 * q2old < 0.0f)
		{
			++kmah;
		}
		rs[it].kmah = kmah;

		++written;
	}

	ray->nrs = written;
	ray->rs = rs;
	ray->nc = 0;
	ray->c = NULL;
	ray->ic = 0;

	return ray;
}

void ray_trace_free_3d_rk4(ray_path_rk_t *ray)
{
	if (ray->c != NULL)
	{
		free1((void *)ray->c);
	}
	free1((void *)ray->rs);
	free1((void *)ray);
}

int ray_trace_nearest_step_3d_rk4(ray_path_rk_t *ray, float x, float y, float z)
{
	int nrs = ray->nrs, ic = ray->ic, nc = ray->nc;
	ray_step_rk_t *rs = ray->rs;
	ray_circle_rk_t *c = (ray_circle_rk_t *)ray->c;
	int irs, irsf, irsl, irsmin = 0, update, jc, js, kc;
	float dsmin, ds, dx, dy, dz, dmin, rdmin, xrs, yrs, zrs;

	if (c == NULL)
	{
		ray->ic = ic = 0;
		ray->nc = nc = sqrtf((float)nrs);
		ray->c = c = rt3d_make_circles_rk4(nc, nrs, rs);
	}

	dx = x - c[ic].x;
	dy = y - c[ic].y;
	dz = z - c[ic].z;
	dmin = 2.0f * (sqrtf(dx * dx + dy * dy + dz * dz) + c[ic].r);
	dsmin = dmin * dmin;

	for (kc = 0, jc = ic, js = 0; kc < nc; ++kc)
	{
		dx = x - c[jc].x;
		dy = y - c[jc].y;
		dz = z - c[jc].z;
		ds = dx * dx + dy * dy + dz * dz;

		rdmin = c[jc].r + dmin;

		if (ds <= rdmin * rdmin)
		{
			irsf = c[jc].irsf;
			irsl = c[jc].irsl;
			update = 0;
			for (irs = irsf; irs <= irsl; ++irs)
			{
				xrs = rs[irs].x;
				yrs = rs[irs].y;
				zrs = rs[irs].z;
				dx = x - xrs;
				dy = y - yrs;
				dz = z - zrs;
				ds = dx * dx + dy * dy + dz * dz;
				if (ds < dsmin)
				{
					dsmin = ds;
					irsmin = irs;
					update = 1;
				}
			}

			if (update)
			{
				dmin = sqrtf(dsmin);
				ic = jc;
			}
		}

		js = (js > 0) ? -js - 1 : -js + 1;
		jc += js;
		if (jc < 0 || jc >= nc)
		{
			js = (js > 0) ? -js - 1 : -js + 1;
			jc += js;
		}
	}

	ray->ic = ic;

	return irsmin;
}
