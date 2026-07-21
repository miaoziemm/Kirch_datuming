#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Boundary
{
    int iz;
    float z;
};

static std::vector<float> parse_depths_fallback(const char *text)
{
    std::vector<float> depths;
    if (text == NULL || text[0] == '\0')
        return depths;
    std::string s(text);
    std::replace(s.begin(), s.end(), ',', ' ');
    std::stringstream ss(s);
    float v = 0.0f;
    while (ss >> v)
        depths.push_back(v);
    return depths;
}

static std::vector<float> read_depths_param()
{
    int ndepth = 0;
    float *arr = se_get_pararray_float("depths", &ndepth);
    std::vector<float> depths;
    if (arr != NULL && ndepth > 0)
    {
        depths.assign(arr, arr + ndepth);
        free1float(arr);
    }
    else
    {
        depths = parse_depths_fallback(se_get_par_str("depths"));
    }
    return depths;
}

static int clamp_index(int i, int n)
{
    if (i < 0)
        return 0;
    if (i >= n)
        return n - 1;
    return i;
}

static float velocity_at(const std::vector<float> &vel, int nz, int nx, int iz, int ix)
{
    iz = clamp_index(iz, nz);
    ix = clamp_index(ix, nx);
    float v = vel[(size_t)ix * (size_t)nz + (size_t)iz];
    if (!(v > 0.0f) || !std::isfinite(v))
        ERROR(("Velocity model contains non-positive or invalid sample at iz=%d ix=%d", iz, ix));
    return v;
}

static float sponge_weight(int ix, int nx, int sponge, float strength)
{
    if (sponge <= 0 || strength <= 0.0f)
        return 1.0f;
    int dist = std::min(ix, nx - 1 - ix);
    if (dist >= sponge)
        return 1.0f;
    float x = (float)(sponge - dist) / (float)sponge;
    return std::exp(-strength * x * x);
}

static std::vector<std::vector<std::complex<float>>> solve_helmholtz_layer(const std::vector<float> &vel,
                                                                            int nz,
                                                                            int nx,
                                                                            int iz_top,
                                                                            int iz_bot,
                                                                            float dx,
                                                                            float dz,
                                                                            float omega,
                                                                            int sponge_n,
                                                                            float sponge_strength,
                                                                            float helmholtz_damping,
                                                                            int max_dense_unknowns,
                                                                            const std::vector<std::complex<float>> &top_field,
                                                                            const std::vector<std::complex<float>> &bottom_field)
{
    const int nzl = iz_bot - iz_top + 1;
    std::vector<std::vector<std::complex<float>>> u((size_t)nzl, std::vector<std::complex<float>>(nx));
    if (nzl <= 1)
    {
        u[0] = top_field;
        return u;
    }

    const float idx2 = 1.0f / (dx * dx);
    const float idz2 = 1.0f / (dz * dz);

    for (int iz = 0; iz < nzl; ++iz)
    {
        const float a = (float)iz / (float)(nzl - 1);
        for (int ix = 0; ix < nx; ++ix)
            u[(size_t)iz][ix] = (1.0f - a) * top_field[ix] + a * bottom_field[ix];
    }
    u.front() = top_field;
    u.back() = bottom_field;

    if (nzl <= 2 || nx <= 2)
        return u;

    const int nxint = nx - 2;
    const int nzint = nzl - 2;
    const int nunk = nxint * nzint;
    if (nunk > max_dense_unknowns)
        ERROR(("dense_direct layer has %d unknowns, exceeding max_dense_unknowns=%d. Increase MAX_DENSE_UNKNOWNS/max_dense_unknowns for this test, reduce the layer thickness, or use more closely spaced DEPTH_INDICES so each dense layer solve remains tractable.", nunk, max_dense_unknowns));

    std::vector<std::complex<float>> A((size_t)nunk * (size_t)nunk, std::complex<float>(0.0f, 0.0f));
    std::vector<std::complex<float>> b((size_t)nunk, std::complex<float>(0.0f, 0.0f));
    auto uid = [nxint](int iz_local, int ix) { return (iz_local - 1) * nxint + (ix - 1); };

    for (int iz = 1; iz < nzl - 1; ++iz)
    {
        const int giz = iz_top + iz;
        for (int ix = 1; ix < nx - 1; ++ix)
        {
            const int row = uid(iz, ix);
            const float v = velocity_at(vel, nz, nx, giz, ix);
            const float k = omega / v;
            const float side = 1.0f - sponge_weight(ix, nx, sponge_n, sponge_strength);
            const std::complex<float> k2 = std::complex<float>(k * k, (helmholtz_damping + side) * k * k);
            A[(size_t)row * nunk + row] = k2 - 2.0f * (idx2 + idz2);

            auto add_neighbor = [&](int niz, int nix, float coeff) {
                if (niz > 0 && niz < nzl - 1 && nix > 0 && nix < nx - 1)
                {
                    const int col = uid(niz, nix);
                    A[(size_t)row * nunk + col] += coeff;
                }
                else
                {
                    b[(size_t)row] -= coeff * u[(size_t)niz][nix];
                }
            };
            add_neighbor(iz, ix - 1, idx2);
            add_neighbor(iz, ix + 1, idx2);
            add_neighbor(iz - 1, ix, idz2);
            add_neighbor(iz + 1, ix, idz2);
        }
    }

    for (int krow = 0; krow < nunk; ++krow)
    {
        int piv = krow;
        float piv_abs = std::abs(A[(size_t)krow * nunk + krow]);
        for (int r = krow + 1; r < nunk; ++r)
        {
            const float a = std::abs(A[(size_t)r * nunk + krow]);
            if (a > piv_abs)
            {
                piv = r;
                piv_abs = a;
            }
        }
        if (piv_abs <= 1.0e-20f || !std::isfinite(piv_abs))
            ERROR(("dense_direct encountered singular/non-finite pivot at row %d", krow));
        if (piv != krow)
        {
            for (int c = krow; c < nunk; ++c)
                std::swap(A[(size_t)krow * nunk + c], A[(size_t)piv * nunk + c]);
            std::swap(b[(size_t)krow], b[(size_t)piv]);
        }
        const std::complex<float> pivot = A[(size_t)krow * nunk + krow];
        for (int r = krow + 1; r < nunk; ++r)
        {
            const std::complex<float> factor = A[(size_t)r * nunk + krow] / pivot;
            if (std::abs(factor) == 0.0f)
                continue;
            A[(size_t)r * nunk + krow] = 0.0f;
            for (int c = krow + 1; c < nunk; ++c)
                A[(size_t)r * nunk + c] -= factor * A[(size_t)krow * nunk + c];
            b[(size_t)r] -= factor * b[(size_t)krow];
        }
    }

    std::vector<std::complex<float>> x((size_t)nunk);
    for (int r = nunk - 1; r >= 0; --r)
    {
        std::complex<float> sum = b[(size_t)r];
        for (int c = r + 1; c < nunk; ++c)
            sum -= A[(size_t)r * nunk + c] * x[(size_t)c];
        x[(size_t)r] = sum / A[(size_t)r * nunk + r];
    }
    for (int iz = 1; iz < nzl - 1; ++iz)
        for (int ix = 1; ix < nx - 1; ++ix)
            u[(size_t)iz][ix] = x[(size_t)uid(iz, ix)];
    INFO(("Helmholtz dense_direct layer solve completed with %d unknowns", nunk));
    return u;
}

int main(int argc, char *argv[])
{
    se_par_init(argc, argv);

    if (!se_have_par("model"))
        ERROR(("Need model= 2-D velocity SEP file"));
    if (!se_have_par("in"))
        ERROR(("Need in= boundary wavefield SEP file from layer_kirchhoff_test"));
    if (!se_have_par("out"))
        ERROR(("Need out= complete wavefield SEP file"));
    if (!se_have_par("depths") && !se_have_par("depth_indices"))
        ERROR(("Need depths= boundary depths or depth_indices= boundary z sample indices matching in.n3"));
    if (!se_have_par("freq"))
        ERROR(("Need freq= frequency in Hz"));

    const char *model_f = se_get_par_str("model");
    const char *in_f = se_get_par_str("in");
    const char *out_f = se_get_par_str("out");
    const float freq = se_get_par_float("freq");
    const int sponge_n = se_have_par("sponge") ? se_get_par_int("sponge") : 0;
    const float sponge_strength = se_have_par("sponge_strength") ? se_get_par_float("sponge_strength") : 0.0f;
    const float helmholtz_damping = se_have_par("helmholtz_damping") ? se_get_par_float("helmholtz_damping") : 0.0f;
    const int max_dense_unknowns = se_have_par("max_dense_unknowns") ? se_get_par_int("max_dense_unknowns") : 20000;

    std::vector<float> depths;
    if (se_have_par("depths"))
        depths = read_depths_param();

    sep_t *model = sep_open(model_f, SEP_READ, 0);
    sep_t *bound = sep_open(in_f, SEP_READ, 0);
    sep_t *out = sep_open(out_f, SEP_WRITE, 0);

    if (model->headers->ndim < 2)
        ERROR(("model= must be a 2-D velocity model with n1=z and n2=x."));
    if (bound->headers->ndim < 3 || bound->headers->n[0] != 2)
        ERROR(("in= must have n1=2(real/imag), n2=nx, n3=number of boundary depths."));

    const int nz = model->headers->n[0];
    const int nx = model->headers->n[1];
    const float oz = (float)model->headers->o[0];
    const float dz = (float)model->headers->d[0];
    const float ox = (float)model->headers->o[1];
    const float dx = (float)model->headers->d[1];
    if (nx != bound->headers->n[1])
        ERROR(("model n2=%d does not match boundary n2=%d", nx, bound->headers->n[1]));
    std::vector<Boundary> bnd;
    if (se_have_par("depth_indices"))
    {
        int nidx = 0;
        int32_t *idx_array = se_get_pararray_int("depth_indices", &nidx);
        if (idx_array == NULL || nidx <= 0)
            ERROR(("No valid indices parsed from depth_indices=."));
        std::vector<int> indices(idx_array, idx_array + nidx);
        free1int(idx_array);
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
        if ((int)indices.size() != bound->headers->n[2] && (int)indices.size() + 1 != bound->headers->n[2])
            ERROR(("depth index count %zu must match boundary n3=%d, or be one smaller when the boundary file includes the model top", indices.size(), bound->headers->n[2]));
        if ((int)indices.size() + 1 == bound->headers->n[2])
            bnd.push_back({0, oz});
        for (int iz : indices)
        {
            if (iz < 0 || iz >= nz)
                ERROR(("Boundary depth index %d outside model range 0..%d", iz, nz - 1));
            if (!bnd.empty() && iz <= bnd.back().iz)
                ERROR(("Boundary depth indices must be strictly increasing."));
            bnd.push_back({iz, oz + iz * dz});
        }
    }
    else
    {
        if (depths.size() < 2)
            ERROR(("Need at least two boundary depths to fill a layer."));
        if ((int)depths.size() != bound->headers->n[2] && (int)depths.size() + 1 != bound->headers->n[2])
            ERROR(("depth count %zu must match boundary n3=%d, or be one smaller when the boundary file includes the model top", depths.size(), bound->headers->n[2]));

        std::vector<float> boundary_depths;
        boundary_depths.reserve((size_t)bound->headers->n[2]);
        if ((int)depths.size() + 1 == bound->headers->n[2])
            boundary_depths.push_back(oz);
        boundary_depths.insert(boundary_depths.end(), depths.begin(), depths.end());
        for (float z : boundary_depths)
        {
            int iz = (int)std::floor((z - oz) / dz + 0.5f);
            if (iz < 0 || iz >= nz)
                ERROR(("Boundary depth %g maps to iz=%d outside model range 0..%d", z, iz, nz - 1));
            if (!bnd.empty() && iz <= bnd.back().iz)
                ERROR(("Boundary depths must be strictly increasing after mapping to samples."));
            bnd.push_back({iz, oz + iz * dz});
        }
    }

    std::vector<float> vel((size_t)nz * (size_t)nx);
    se_fsio_read_float(model->data->io, vel.data(), vel.size());

    std::vector<float> raw_bound((size_t)2 * (size_t)nx * bnd.size());
    se_fsio_read_float(bound->data->io, raw_bound.data(), raw_bound.size());

    std::vector<std::vector<std::complex<float>>> boundary_fields(bnd.size(), std::vector<std::complex<float>>(nx));
    for (size_t ib = 0; ib < bnd.size(); ++ib)
    {
        const size_t base = ib * (size_t)2 * (size_t)nx;
        for (int ix = 0; ix < nx; ++ix)
            boundary_fields[ib][ix] = std::complex<float>(raw_bound[base + (size_t)2 * ix],
                                                          raw_bound[base + (size_t)2 * ix + 1]);
    }

    std::vector<float> sponge(nx, 1.0f);
    for (int ix = 0; ix < nx; ++ix)
        sponge[ix] = sponge_weight(ix, nx, sponge_n, sponge_strength);

    const int iz0 = bnd.front().iz;
    const int iz1 = bnd.back().iz;
    const int nzout = iz1 - iz0 + 1;
    out->headers->ndim = 3;
    out->headers->n[0] = 2;
    out->headers->n[1] = nx;
    out->headers->n[2] = nzout;
    out->headers->d[0] = 1.0;
    out->headers->d[1] = dx;
    out->headers->d[2] = dz;
    out->headers->o[0] = 0.0;
    out->headers->o[1] = ox;
    out->headers->o[2] = oz + iz0 * dz;

    std::vector<float> row((size_t)2 * (size_t)nx);
    const float omega = 2.0f * (float)M_PI * freq;
    for (size_t il = 0; il + 1 < bnd.size(); ++il)
    {
        const int top = bnd[il].iz;
        const int bot = bnd[il + 1].iz;
        std::vector<std::vector<std::complex<float>>> layer_field =
            solve_helmholtz_layer(vel, nz, nx, top, bot, dx, dz, omega,
                                  sponge_n, sponge_strength, helmholtz_damping, max_dense_unknowns,
                                  boundary_fields[il], boundary_fields[il + 1]);

        const int first_local = (il == 0) ? 0 : 1;
        for (int loc = first_local; loc <= bot - top; ++loc)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                const std::complex<float> &u = layer_field[(size_t)loc][ix];
                row[(size_t)2 * ix] = u.real();
                row[(size_t)2 * ix + 1] = u.imag();
            }
            se_fsio_write_float(out->data->io, row.data(), row.size());
        }
        INFO(("Filled layer %zu between depth=%g and depth=%g", il + 1, bnd[il].z, bnd[il + 1].z));
    }

    sep_close(model);
    sep_close(bound);
    sep_close(out);
    se_par_destroy();
    return 0;
}
