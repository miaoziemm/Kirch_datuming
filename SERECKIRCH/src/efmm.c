#include "efmm.h"

#include <string.h>

#ifndef TT_INF
#define TT_INF 1.0e30
#endif

static efmm_t efmm_eikods_ctx;
static int efmm_eikods_ctx_ready = 0;
static int efmm_eikods_ctx_n1 = 0;
static int efmm_eikods_ctx_n2 = 0;
static int efmm_eikods_ctx_n3 = 0;

typedef struct
{
    int *iz;
    int *ix;
    double *t;
    size_t n;
    size_t cap;
} narrow_band_t;

static int band_init(narrow_band_t *b, size_t cap)
{
    b->iz = (int *)malloc(sizeof(int) * cap);
    b->ix = (int *)malloc(sizeof(int) * cap);
    b->t = (double *)malloc(sizeof(double) * cap);
    if (!b->iz || !b->ix || !b->t)
    {
        free(b->iz);
        free(b->ix);
        free(b->t);
        b->iz = NULL;
        b->ix = NULL;
        b->t = NULL;
        b->n = 0;
        b->cap = 0;
        return -1;
    }
    b->n = 0;
    b->cap = cap;
    return 0;
}

static void band_free(narrow_band_t *b)
{
    free(b->iz);
    free(b->ix);
    free(b->t);
    b->iz = NULL;
    b->ix = NULL;
    b->t = NULL;
    b->n = 0;
    b->cap = 0;
}

static int band_reserve_one(narrow_band_t *b)
{
    if (b->n < b->cap)
    {
        return 0;
    }

    size_t ncap = (b->cap == 0) ? 1024 : (2 * b->cap);
    int *niz = (int *)realloc(b->iz, sizeof(int) * ncap);
    int *nix = (int *)realloc(b->ix, sizeof(int) * ncap);
    double *nt = (double *)realloc(b->t, sizeof(double) * ncap);
    if (!niz || !nix || !nt)
    {
        free(niz);
        free(nix);
        free(nt);
        return -1;
    }

    b->iz = niz;
    b->ix = nix;
    b->t = nt;
    b->cap = ncap;
    return 0;
}

static void band_bubble_left(narrow_band_t *b, size_t pos)
{
    while (pos > 0 && b->t[pos - 1] > b->t[pos])
    {
        double tt = b->t[pos - 1];
        int zi = b->iz[pos - 1];
        int xi = b->ix[pos - 1];
        b->t[pos - 1] = b->t[pos];
        b->iz[pos - 1] = b->iz[pos];
        b->ix[pos - 1] = b->ix[pos];
        b->t[pos] = tt;
        b->iz[pos] = zi;
        b->ix[pos] = xi;
        pos--;
    }
}

static int band_insert(narrow_band_t *b, int iz, int ix, double t)
{
    if (band_reserve_one(b) != 0)
    {
        return -1;
    }

    size_t pos = b->n;
    b->iz[pos] = iz;
    b->ix[pos] = ix;
    b->t[pos] = t;
    b->n++;
    band_bubble_left(b, pos);
    return 0;
}

static int band_decrease_key(narrow_band_t *b, int iz, int ix, double tnew)
{
    for (size_t p = 0; p < b->n; ++p)
    {
        if (b->iz[p] == iz && b->ix[p] == ix)
        {
            b->t[p] = tnew;
            band_bubble_left(b, p);
            return 0;
        }
    }
    return -1;
}

static int band_pop_min(narrow_band_t *b, int *iz, int *ix, double *t)
{
    if (b->n == 0)
    {
        return -1;
    }

    *iz = b->iz[0];
    *ix = b->ix[0];
    *t = b->t[0];
    if (b->n > 1)
    {
        memmove(&b->iz[0], &b->iz[1], sizeof(int) * (b->n - 1));
        memmove(&b->ix[0], &b->ix[1], sizeof(int) * (b->n - 1));
        memmove(&b->t[0], &b->t[1], sizeof(double) * (b->n - 1));
    }
    b->n--;
    return 0;
}

static inline double min2(double a, double b)
{
    return (a < b) ? a : b;
}

static inline size_t idx_colmajor(int iz, int ix, int nz)
{
    return (size_t)iz + (size_t)ix * (size_t)nz;
}

static inline double classic_update_aniso(double tx,
                                          int has_x,
                                          double tz,
                                          int has_z,
                                          double s0,
                                          double dx,
                                          double dz)
{
    double best = TT_INF;

    if (has_x)
    {
        best = min2(best, tx + s0 * dx);
    }
    if (has_z)
    {
        best = min2(best, tz + s0 * dz);
    }

    if (has_x && has_z)
    {
        double invdx2 = 1.0 / (dx * dx);
        double invdz2 = 1.0 / (dz * dz);
        double A = invdx2 + invdz2;
        double B = -2.0 * (tx * invdx2 + tz * invdz2);
        double C = tx * tx * invdx2 + tz * tz * invdz2 - s0 * s0;
        double disc = B * B - 4.0 * A * C;
        if (disc < 0.0)
        {
            disc = 0.0;
        }
        double root = (-B + sqrt(disc)) / (2.0 * A);
        if (isfinite(root) && root >= ((tx > tz) ? tx : tz))
        {
            best = min2(best, root);
        }
    }

    return best;
}

static void recompute_eff_point_aniso(int iz,
                                      int ix,
                                      int nz,
                                      int nx,
                                      const double *aliveT,
                                      const double *aliveE,
                                      double s0,
                                      double dx,
                                      double dz,
                                      const double *dist,
                                      const double *dx_d,
                                      const double *dz_d,
                                      double *tnew,
                                      double *enew)
{
    const double epsd = 1e-12;

    size_t id = idx_colmajor(iz, ix, nz);
    double d = dist[id];
    if (d < epsd)
    {
        *enew = s0;
        *tnew = 0.0;
        return;
    }

    int has_x = 0;
    double ax = 0.0;
    double bx = 0.0;
    double Tx = -TT_INF;

    double tL = TT_INF;
    double tR = TT_INF;
    if (ix - 1 >= 0)
    {
        tL = aliveT[idx_colmajor(iz, ix - 1, nz)];
    }
    if (ix + 1 < nx)
    {
        tR = aliveT[idx_colmajor(iz, ix + 1, nz)];
    }

    if (tL < TT_INF || tR < TT_INF)
    {
        has_x = 1;
        if (tL <= tR)
        {
            double Eup = aliveE[idx_colmajor(iz, ix - 1, nz)];
            ax = d / dx + dx_d[id];
            bx = (d / dx) * Eup;
            Tx = tL;
        }
        else
        {
            double Eup = aliveE[idx_colmajor(iz, ix + 1, nz)];
            ax = d / dx - dx_d[id];
            bx = (d / dx) * Eup;
            Tx = tR;
        }
    }

    int has_z = 0;
    double az = 0.0;
    double bz = 0.0;
    double Tz = -TT_INF;

    double tU = TT_INF;
    double tD = TT_INF;
    if (iz - 1 >= 0)
    {
        tU = aliveT[idx_colmajor(iz - 1, ix, nz)];
    }
    if (iz + 1 < nz)
    {
        tD = aliveT[idx_colmajor(iz + 1, ix, nz)];
    }

    if (tU < TT_INF || tD < TT_INF)
    {
        has_z = 1;
        if (tU <= tD)
        {
            double Eup = aliveE[idx_colmajor(iz - 1, ix, nz)];
            az = d / dz + dz_d[id];
            bz = (d / dz) * Eup;
            Tz = tU;
        }
        else
        {
            double Eup = aliveE[idx_colmajor(iz + 1, ix, nz)];
            az = d / dz - dz_d[id];
            bz = (d / dz) * Eup;
            Tz = tD;
        }
    }

    if (!has_x && !has_z)
    {
        *tnew = TT_INF;
        *enew = TT_INF;
        return;
    }

    double candE[2];
    int nc = 0;
    if (has_x && has_z)
    {
        double A = ax * ax + az * az;
        double B = -2.0 * (ax * bx + az * bz);
        double C = bx * bx + bz * bz - s0 * s0;
        double disc = B * B - 4.0 * A * C;
        if (disc < 0.0)
        {
            disc = 0.0;
        }
        double sq = sqrt(disc);
        candE[nc++] = (-B + sq) / (2.0 * A);
        candE[nc++] = (-B - sq) / (2.0 * A);
    }
    else if (has_x)
    {
        if (fabs(ax) < 1e-14)
        {
            *tnew = TT_INF;
            *enew = TT_INF;
            return;
        }
        candE[nc++] = (bx + s0) / ax;
        candE[nc++] = (bx - s0) / ax;
    }
    else
    {
        if (fabs(az) < 1e-14)
        {
            *tnew = TT_INF;
            *enew = TT_INF;
            return;
        }
        candE[nc++] = (bz + s0) / az;
        candE[nc++] = (bz - s0) / az;
    }

    double bestT = TT_INF;
    double bestE = TT_INF;
    double Tref = -TT_INF;
    if (has_x)
    {
        Tref = (Tref > Tx) ? Tref : Tx;
    }
    if (has_z)
    {
        Tref = (Tref > Tz) ? Tref : Tz;
    }

    for (int k = 0; k < nc; ++k)
    {
        double E = candE[k];
        if (!isfinite(E) || E <= 0.0)
        {
            continue;
        }
        double T = d * E;
        if (T + 1e-12 < Tref)
        {
            continue;
        }
        if (T < bestT)
        {
            bestT = T;
            bestE = E;
        }
    }

    if (bestT >= TT_INF)
    {
        double tx = TT_INF;
        double tz = TT_INF;

        if (ix - 1 >= 0)
        {
            tx = min2(tx, aliveT[idx_colmajor(iz, ix - 1, nz)]);
        }
        if (ix + 1 < nx)
        {
            tx = min2(tx, aliveT[idx_colmajor(iz, ix + 1, nz)]);
        }
        if (iz - 1 >= 0)
        {
            tz = min2(tz, aliveT[idx_colmajor(iz - 1, ix, nz)]);
        }
        if (iz + 1 < nz)
        {
            tz = min2(tz, aliveT[idx_colmajor(iz + 1, ix, nz)]);
        }

        int has_xn = (tx < TT_INF);
        int has_zn = (tz < TT_INF);
        double T = classic_update_aniso(tx, has_xn, tz, has_zn, s0, dx, dz);
        if (T < TT_INF)
        {
            bestT = T;
            bestE = T / d;
        }
    }

    *tnew = bestT;
    *enew = bestE;
}


int efmm_set_src(efmm_t *efmm, float src_x, float src_z)
{
    efmm->src_x = src_x;
    efmm->src_z = src_z;
    return 0;
}

int efmm_set_vel(efmm_t *efmm, float *vel)
{
    int size = efmm->nx_work * efmm->nz_work;
    for (int i = 0; i < size; i++)
    {
        efmm->vel[i] = vel[i];
    }
    return 0;
}


int efmm_init(efmm_t *efmm, int nx_work, int nz_work, int dx, int dz, float src_x, float src_z)
{
    efmm->nx_work = nx_work;
    efmm->nz_work = nz_work;
    efmm->dx = (float)dx;
    efmm->dz = (float)dz;
    efmm->src_x = src_x;
    efmm->src_z = src_z;

    int size = nx_work * nz_work;
    efmm->vel = (float *)malloc(size * sizeof(float));
    if (efmm->vel == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for velocity model.\n");
        return -1;
    }

    efmm->tt = (float *)malloc(size * sizeof(float));
    if (efmm->tt == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for travel time field.\n");
        free(efmm->vel);
        efmm->vel = NULL;
        efmm->tt = NULL;
        return -1;
    }

    return 0;

}


int efmm_free(efmm_t *efmm){
    free(efmm->vel);
    free(efmm->tt);
    efmm->vel = NULL;
    efmm->tt = NULL;
    return 0;
}




int efmm_solver(efmm_t *efmm)
{
    if (!efmm || !efmm->vel || !efmm->tt)
    {
        return -1;
    }

    int nz = efmm->nz_work;
    int nx = efmm->nx_work;
    if (nz <= 0 || nx <= 0)
    {
        return -2;
    }

    double dx = (double)efmm->dx;
    double dz = (double)efmm->dz;
    if (dx <= 0.0 || dz <= 0.0)
    {
        return -3;
    }

    size_t n = (size_t)nz * (size_t)nx;
    const double epsd = 1e-12;

    double *slow = (double *)malloc(sizeof(double) * n);
    double *dist = (double *)malloc(sizeof(double) * n);
    double *dx_d = (double *)malloc(sizeof(double) * n);
    double *dz_d = (double *)malloc(sizeof(double) * n);
    double *aliveT = (double *)malloc(sizeof(double) * n);
    double *recpt = (double *)malloc(sizeof(double) * n);
    double *aliveE = (double *)malloc(sizeof(double) * n);
    double *recpE = (double *)malloc(sizeof(double) * n);

    if (!slow || !dist || !dx_d || !dz_d || !aliveT || !recpt || !aliveE || !recpE)
    {
        free(slow);
        free(dist);
        free(dx_d);
        free(dz_d);
        free(aliveT);
        free(recpt);
        free(aliveE);
        free(recpE);
        return -4;
    }

    double src_x = (double)efmm->src_x;
    double src_z = (double)efmm->src_z;
    double max_x = (double)(nx - 1) * dx;
    double max_z = (double)(nz - 1) * dz;

    if (src_x < 0.0)
    {
        src_x = 0.0;
    }
    else if (src_x > max_x)
    {
        src_x = max_x;
    }

    if (src_z < 0.0)
    {
        src_z = 0.0;
    }
    else if (src_z > max_z)
    {
        src_z = max_z;
    }

    for (int ix = 0; ix < nx; ++ix)
    {
        double x = (double)ix * dx;
        for (int iz = 0; iz < nz; ++iz)
        {
            size_t id = idx_colmajor(iz, ix, nz);
            float v = efmm->vel[id];
            if (v <= 0.0f)
            {
                free(slow);
                free(dist);
                free(dx_d);
                free(dz_d);
                free(aliveT);
                free(recpt);
                free(aliveE);
                free(recpE);
                return -5;
            }

            slow[id] = 1.0 / (double)v;

            double z = (double)iz * dz;
            double rx = x - src_x;
            double rz = z - src_z;
            double d = sqrt(rx * rx + rz * rz);
            dist[id] = d;
            if (d > epsd)
            {
                double inv = 1.0 / d;
                dx_d[id] = rx * inv;
                dz_d[id] = rz * inv;
            }
            else
            {
                dx_d[id] = 0.0;
                dz_d[id] = 0.0;
            }

            aliveT[id] = TT_INF;
            recpt[id] = TT_INF;
            aliveE[id] = TT_INF;
            recpE[id] = TT_INF;
        }
    }

    narrow_band_t band;
    if (band_init(&band, (size_t)(8 * nz + 8)) != 0)
    {
        free(slow);
        free(dist);
        free(dx_d);
        free(dz_d);
        free(aliveT);
        free(recpt);
        free(aliveE);
        free(recpE);
        return -6;
    }

    /* Bilinear source seeding: initialize four surrounding grid vertices. */
    double gx = src_x / dx;
    double gz = src_z / dz;
    int ix0 = (int)floor(gx);
    int iz0 = (int)floor(gz);
    if (ix0 < 0)
    {
        ix0 = 0;
    }
    if (iz0 < 0)
    {
        iz0 = 0;
    }
    if (ix0 >= nx)
    {
        ix0 = nx - 1;
    }
    if (iz0 >= nz)
    {
        iz0 = nz - 1;
    }
    int ix1 = (ix0 + 1 < nx) ? (ix0 + 1) : ix0;
    int iz1 = (iz0 + 1 < nz) ? (iz0 + 1) : iz0;

    const int seed_iz[4] = {iz0, iz0, iz1, iz1};
    const int seed_ix[4] = {ix0, ix1, ix0, ix1};
    for (int s = 0; s < 4; ++s)
    {
        int i = seed_iz[s];
        int j = seed_ix[s];
        size_t id = idx_colmajor(i, j, nz);
        double d0 = dist[id];
        double t0 = d0 * slow[id];
        if (t0 < recpt[id])
        {
            recpt[id] = t0;
            recpE[id] = (d0 > epsd) ? (t0 / d0) : slow[id];
        }
    }

    for (int ix = 0; ix < nx; ++ix)
    {
        for (int iz = 0; iz < nz; ++iz)
        {
            size_t id = idx_colmajor(iz, ix, nz);
            if (recpt[id] < TT_INF)
            {
                if (band_insert(&band, iz, ix, recpt[id]) != 0)
                {
                    band_free(&band);
                    free(slow);
                    free(dist);
                    free(dx_d);
                    free(dz_d);
                    free(aliveT);
                    free(recpt);
                    free(aliveE);
                    free(recpE);
                    return -7;
                }
            }
        }
    }

    int max_iter = 4 * nz * nx + 100;
    int iter = 0;
    while (band.n > 0 && iter < max_iter)
    {
        iter++;
        int i;
        int j;
        double tmin;
        if (band_pop_min(&band, &i, &j, &tmin) != 0)
        {
            break;
        }

        size_t id = idx_colmajor(i, j, nz);
        if (aliveT[id] != TT_INF)
        {
            continue;
        }
        if (tmin > recpt[id])
        {
            continue;
        }

        aliveT[id] = recpt[id];
        if (dist[id] > epsd)
        {
            aliveE[id] = (recpE[id] < TT_INF) ? recpE[id] : (recpt[id] / dist[id]);
        }
        else
        {
            aliveE[id] = slow[id];
        }

        const int ni[4] = {i - 1, i + 1, i, i};
        const int nj[4] = {j, j, j - 1, j + 1};
        for (int k = 0; k < 4; ++k)
        {
            int ii = ni[k];
            int jj = nj[k];
            if (ii < 0 || ii >= nz || jj < 0 || jj >= nx)
            {
                continue;
            }
            size_t nid = idx_colmajor(ii, jj, nz);
            if (aliveT[nid] != TT_INF)
            {
                continue;
            }

            double oldt = recpt[nid];
            double tnew;
            double enew;
            recompute_eff_point_aniso(ii, jj, nz, nx, aliveT, aliveE, slow[nid], dx, dz, dist, dx_d, dz_d, &tnew, &enew);
            if (tnew >= TT_INF)
            {
                continue;
            }

            if (oldt == TT_INF)
            {
                recpt[nid] = tnew;
                recpE[nid] = enew;
                if (band_insert(&band, ii, jj, tnew) != 0)
                {
                    band_free(&band);
                    free(slow);
                    free(dist);
                    free(dx_d);
                    free(dz_d);
                    free(aliveT);
                    free(recpt);
                    free(aliveE);
                    free(recpE);
                    return -8;
                }
            }
            else if (tnew < oldt)
            {
                recpt[nid] = tnew;
                recpE[nid] = enew;
                if (band_decrease_key(&band, ii, jj, tnew) != 0)
                {
                    band_free(&band);
                    free(slow);
                    free(dist);
                    free(dx_d);
                    free(dz_d);
                    free(aliveT);
                    free(recpt);
                    free(aliveE);
                    free(recpE);
                    return -8;
                }
            }
        }
    }

    for (size_t id = 0; id < n; ++id)
    {
        efmm->tt[id] = (float)recpt[id];
    }

    band_free(&band);
    free(slow);
    free(dist);
    free(dx_d);
    free(dz_d);
    free(aliveT);
    free(recpt);
    free(aliveE);
    free(recpE);

    return 0;
}

void efmm_get_bottom_tt(efmm_t *efmm, float *bottom_tt)
{
    int nx = efmm->nx_work;
    int nz = efmm->nz_work;
    for (int ix = 0; ix < nx; ++ix)
    {
        bottom_tt[ix] = efmm->tt[idx_colmajor(nz - 1, ix, nz)];
    }
}

void efmm_get_top_tt(efmm_t *efmm, float *top_tt)
{
    int nx = efmm->nx_work;
    int nz = efmm->nz_work;
    for (int ix = 0; ix < nx; ++ix)
    {
        top_tt[ix] = efmm->tt[idx_colmajor(0, ix, nz)];
    }
}

static int efmm_eikods_solve_at_source(float src_x, float src_z, float *tt_out)
{
    size_t n;

    if (!efmm_eikods_ctx_ready || tt_out == NULL)
    {
        return -1;
    }

    n = (size_t)efmm_eikods_ctx.nx_work * (size_t)efmm_eikods_ctx.nz_work;
    efmm_eikods_ctx.src_x = src_x;
    efmm_eikods_ctx.src_z = src_z;

    if (efmm_solver(&efmm_eikods_ctx) != 0)
    {
        return -1;
    }

    for (size_t i = 0; i < n; ++i)
    {
        tt_out[i] = efmm_eikods_ctx.tt[i];
    }

    return 0;
}

void efmm_eikods_init (int n3,int n2,int n1)
/*< Initialize data dimensions >*/
{
    if (n1 <= 0 || n2 <= 0 || n3 <= 0)
    {
        if (efmm_eikods_ctx_ready)
        {
            efmm_free(&efmm_eikods_ctx);
            efmm_eikods_ctx_ready = 0;
        }
        efmm_eikods_ctx_n1 = 0;
        efmm_eikods_ctx_n2 = 0;
        efmm_eikods_ctx_n3 = 0;
        return;
    }

    if (efmm_eikods_ctx_ready &&
        efmm_eikods_ctx_n1 == n1 &&
        efmm_eikods_ctx_n2 == n2 &&
        efmm_eikods_ctx_n3 == n3)
    {
        return;
    }

    if (efmm_eikods_ctx_ready)
    {
        efmm_free(&efmm_eikods_ctx);
        efmm_eikods_ctx_ready = 0;
    }

    efmm_eikods_ctx_n1 = n1;
    efmm_eikods_ctx_n2 = n2;
    efmm_eikods_ctx_n3 = n3;

    if (efmm_init(&efmm_eikods_ctx, n1, n2, 1, 1, 0.0f, 0.0f) == 0)
    {
        efmm_eikods_ctx_ready = 1;
    }
}

void efmm_eikods (float* time                /* time */, 
	     float* v                   /* slowness squared */, 
	     int* in                    /* in/front/out flag */, 
	     bool* plane                /* if plane source */, 
	     int   n3,  int n2,  int n1 /* dimensions */, 
	     float o3,float o2,float o1 /* origin */, 
	     float d3,float d2,float d1 /* sampling */, 
	     float s3,float s2,float s1 /* source */, 
	     int   b3,  int b2,  int b1 /* box around the source */, 
	     int order                  /* accuracy order (1,2,3) */, 
	     int l                      /* direction of source perturbation */, 
	     float* dl1, float* ds1     /* first-order derivatives */, 
	     float* dl2, float* ds2     /* second-order derivatives */)
/*< Run fast marching eikonal solver >*/
{
    size_t n;
    float *vel;
    float *tt0;
    float *ttp;
    float *ttm;
    float base_src_x;
    float base_src_z;
    int axis;
    float step;
    int need_first;
    int need_second;

    (void)in;
    (void)plane;
    (void)o3;
    (void)d3;
    (void)s3;
    (void)b3;
    (void)b2;
    (void)b1;
    (void)order;
    (void)l;

    if (time == NULL || v == NULL || n1 <= 0 || n2 <= 0 || n3 <= 0)
    {
        return;
    }

    if (n3 != 1)
    {
        return;
    }

    efmm_eikods_init(n3, n2, n1);
    if (!efmm_eikods_ctx_ready)
    {
        return;
    }

    efmm_eikods_ctx.dx = d1;
    efmm_eikods_ctx.dz = d2;
    base_src_x = s1 - o1;
    base_src_z = s2 - o2;

    n = (size_t)n1 * (size_t)n2;
    vel = (float *)malloc(sizeof(float) * n);
    if (vel == NULL)
    {
        return;
    }

    tt0 = (float *)malloc(sizeof(float) * n);
    if (tt0 == NULL)
    {
        free(vel);
        return;
    }

    need_first = (dl1 != NULL || ds1 != NULL);
    need_second = (dl2 != NULL || ds2 != NULL);
    ttp = NULL;
    ttm = NULL;
    if (need_first || need_second)
    {
        ttp = (float *)malloc(sizeof(float) * n);
        ttm = (float *)malloc(sizeof(float) * n);
        if (ttp == NULL || ttm == NULL)
        {
            free(vel);
            free(tt0);
            free(ttp);
            free(ttm);
            return;
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        if (v[i] > 0.0f)
        {
            vel[i] = 1.0f / sqrtf(v[i]);
        }
        else
        {
            vel[i] = 0.0f;
        }
    }

    if (efmm_set_vel(&efmm_eikods_ctx, vel) != 0 || efmm_solver(&efmm_eikods_ctx) != 0)
    {
        free(vel);
        free(tt0);
        free(ttp);
        free(ttm);
        return;
    }

    free(vel);

    for (size_t i = 0; i < n; ++i)
    {
        tt0[i] = efmm_eikods_ctx.tt[i];
        time[i] = tt0[i];
    }

    if (need_first || need_second)
    {
        axis = l;
        if (axis < 0)
        {
            axis = 0;
        }
        if (axis > 1)
        {
            axis = 1;
        }

        step = (axis == 0) ? d1 : d2;
        if (step <= 0.0f)
        {
            step = 1.0f;
        }

        if (efmm_eikods_solve_at_source(base_src_x + ((axis == 0) ? step : 0.0f),
                                        base_src_z + ((axis == 1) ? step : 0.0f),
                                        ttp) != 0 ||
            efmm_eikods_solve_at_source(base_src_x - ((axis == 0) ? step : 0.0f),
                                        base_src_z - ((axis == 1) ? step : 0.0f),
                                        ttm) != 0)
        {
            if (need_first)
            {
                if (dl1 != NULL)
                {
                    memset(dl1, 0, sizeof(float) * n);
                }
                if (ds1 != NULL)
                {
                    memset(ds1, 0, sizeof(float) * n);
                }
            }
            if (need_second)
            {
                if (dl2 != NULL)
                {
                    memset(dl2, 0, sizeof(float) * n);
                }
                if (ds2 != NULL)
                {
                    memset(ds2, 0, sizeof(float) * n);
                }
            }
            free(tt0);
            free(ttp);
            free(ttm);
            return;
        }

        if (need_first)
        {
            for (size_t i = 0; i < n; ++i)
            {
                float first = (ttp[i] - ttm[i]) / (2.0f * step);
                if (dl1 != NULL)
                {
                    dl1[i] = first;
                }
                if (ds1 != NULL)
                {
                    ds1[i] = -first;
                }
            }
        }

        if (need_second)
        {
            float inv_step2 = 1.0f / (step * step);
            for (size_t i = 0; i < n; ++i)
            {
                float second = (ttp[i] - 2.0f * tt0[i] + ttm[i]) * inv_step2;
                if (dl2 != NULL)
                {
                    dl2[i] = second;
                }
                if (ds2 != NULL)
                {
                    ds2[i] = second;
                }
            }
        }
    }

    free(tt0);
    free(ttp);
    free(ttm);
}