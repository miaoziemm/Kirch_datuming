#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct LayerBoundary
{
    int iz;
    float z;
};

static std::vector<float> parse_depths(const char *text)
{
    std::vector<float> depths;
    if (text == NULL || text[0] == '\0')
        return depths;

    std::string s(text);
    std::replace(s.begin(), s.end(), ',', ' ');
    std::stringstream ss(s);
    float z = 0.0f;
    while (ss >> z)
        depths.push_back(z);
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

static float segment_traveltime(const std::vector<float> &vel,
                                int nz,
                                int nx,
                                float oz,
                                float dz,
                                float ox,
                                float dx,
                                int iz_top,
                                int iz_bot,
                                int ix_top,
                                int ix_bot)
{
    if (iz_bot <= iz_top)
        return 0.0f;

    const float x_top = ox + ix_top * dx;
    const float x_bot = ox + ix_bot * dx;
    const float z_top = oz + iz_top * dz;
    const float z_bot = oz + iz_bot * dz;
    const float length = std::hypot(x_bot - x_top, z_bot - z_top);
    const int nstep = std::max(2, iz_bot - iz_top + std::abs(ix_bot - ix_top) + 1);

    float sum_slowness = 0.0f;
    for (int is = 0; is < nstep; ++is)
    {
        const float a = (nstep == 1) ? 0.0f : (float)is / (float)(nstep - 1);
        const float x = x_top + a * (x_bot - x_top);
        const float z = z_top + a * (z_bot - z_top);
        const int ix = (int)std::floor((x - ox) / dx + 0.5f);
        const int iz = (int)std::floor((z - oz) / dz + 0.5f);
        sum_slowness += 1.0f / velocity_at(vel, nz, nx, iz, ix);
    }
    return length * sum_slowness / (float)nstep;
}

static void build_layer_traveltime(const std::vector<float> &vel,
                                   int nz,
                                   int nx,
                                   float oz,
                                   float dz,
                                   float ox,
                                   float dx,
                                   int iz_top,
                                   int iz_bot,
                                   float **tau)
{
    for (int ix_top = 0; ix_top < nx; ++ix_top)
    {
        for (int ix_bot = 0; ix_bot < nx; ++ix_bot)
        {
            tau[ix_bot][ix_top] = segment_traveltime(vel, nz, nx, oz, dz, ox, dx,
                                                     iz_top, iz_bot, ix_top, ix_bot);
        }
    }
}

static void build_kirchhoff_amp(float **tau,
                                int nx,
                                float dx,
                                fftwf_complex *amp)
{
    for (int ix_bot = 0; ix_bot < nx; ++ix_bot)
    {
        for (int ix_top = 0; ix_top < nx; ++ix_top)
        {
            const float t = std::max(tau[ix_bot][ix_top], 1.0e-6f);
            const float scale = dx / std::sqrt(t);
            const size_t k = (size_t)ix_bot * (size_t)nx + (size_t)ix_top;
            amp[k][0] = scale;
            amp[k][1] = 0.0f;
        }
    }
}

static void propagate_direct(float **tau,
                             int nx,
                             float dx,
                             float omega,
                             const std::vector<std::complex<float>> &uin,
                             std::vector<std::complex<float>> &uout)
{
    std::fill(uout.begin(), uout.end(), std::complex<float>(0.0f, 0.0f));
    for (int ix_bot = 0; ix_bot < nx; ++ix_bot)
    {
        std::complex<float> sum(0.0f, 0.0f);
        for (int ix_top = 0; ix_top < nx; ++ix_top)
        {
            const float t = std::max(tau[ix_bot][ix_top], 1.0e-6f);
            const float phase = omega * t;
            const std::complex<float> g = (dx / std::sqrt(t)) *
                                          std::complex<float>(std::cos(phase), -std::sin(phase));
            sum += g * uin[ix_top];
        }
        uout[ix_bot] = sum;
    }
}

static void propagate_butterfly(float **tau,
                                int nx,
                                float dx,
                                float omega,
                                int bf_p,
                                int bf_leaf,
                                int bf_panel_levels,
                                float bf_amp_eps,
                                float bf_phase_tol,
                                const std::vector<std::complex<float>> &uin,
                                std::vector<std::complex<float>> &uout)
{
    fftwf_complex *amp = (fftwf_complex *)fftwf_malloc((size_t)nx * (size_t)nx * sizeof(fftwf_complex));
    fftwf_complex *bin = (fftwf_complex *)fftwf_malloc((size_t)nx * sizeof(fftwf_complex));
    fftwf_complex *bout = (fftwf_complex *)fftwf_malloc((size_t)nx * sizeof(fftwf_complex));
    if (amp == NULL || bin == NULL || bout == NULL)
        ERROR(("Out of memory for butterfly propagation."));

    build_kirchhoff_amp(tau, nx, dx, amp);
    for (int ix = 0; ix < nx; ++ix)
    {
        bin[ix][0] = uin[ix].real();
        bin[ix][1] = uin[ix].imag();
    }

    BFStrictSegmentedFactor *factor = bf1d_strict_segmented_create_phase_amp(nx, tau, amp, omega,
                                                                              bf_p, bf_leaf,
                                                                              bf_panel_levels,
                                                                              bf_amp_eps,
                                                                              bf_phase_tol);
    bf1d_strict_segmented_apply(factor, bin, bout);
    bf1d_strict_segmented_destroy(factor);

    for (int ix = 0; ix < nx; ++ix)
        uout[ix] = std::complex<float>(bout[ix][0], bout[ix][1]);

    fftwf_free(amp);
    fftwf_free(bin);
    fftwf_free(bout);
}

static void initialize_hankel_wavefield(int nx,
                                        float ox,
                                        float dx,
                                        float srcx,
                                        float srcz,
                                        float topz,
                                        float freq,
                                        float vref,
                                        std::vector<std::complex<float>> &u)
{
    const float k = 2.0f * (float)M_PI * freq / vref;
    for (int ix = 0; ix < nx; ++ix)
    {
        const float x = ox + ix * dx;
        const float r = std::max(std::hypot(x - srcx, topz - srcz), 0.25f * std::fabs(dx));
        const float kr = std::max(k * r, 1.0e-4f);
        u[ix] = std::complex<float>((float)std::cyl_bessel_j(0, kr), (float)std::cyl_neumann(0, kr));
    }
}

int main(int argc, char *argv[])
{
    se_par_init(argc, argv);

    if (!se_have_par("in"))
        ERROR(("Need in= velocity-model SEP file"));
    if (!se_have_par("out"))
        ERROR(("Need out= output wavefield SEP file"));
    if (!se_have_par("depths") && !se_have_par("depth_indices"))
        ERROR(("Need depths= output depths or depth_indices= output z sample indices"));
    if (!se_have_par("freq"))
        ERROR(("Need freq= frequency in Hz"));

    const char *in_f = se_get_par_str("in");
    const char *out_f = se_get_par_str("out");
    int ndepth = 0;
    float *depth_array = se_have_par("depths") ? se_get_pararray_float("depths", &ndepth) : NULL;
    std::vector<float> depths;
    if (depth_array != NULL && ndepth > 0)
    {
        depths.assign(depth_array, depth_array + ndepth);
        free1float(depth_array);
    }
    else if (se_have_par("depths"))
    {
        depths = parse_depths(se_get_par_str("depths"));
    }
    const float freq = se_get_par_float("freq");
    const int use_bf = se_have_par("use_bf") ? se_get_par_int("use_bf") : 1;
    const int include_top = se_have_par("include_top") ? se_get_par_int("include_top") : 1;
    const int bf_p = se_have_par("bf_p") ? se_get_par_int("bf_p") : 12;
    const int bf_leaf = se_have_par("bf_n_leaf") ? se_get_par_int("bf_n_leaf") : 16;
    const int bf_panel_levels = se_have_par("bf_panel_levels") ? se_get_par_int("bf_panel_levels") : 1;
    const float bf_amp_eps = se_have_par("bf_amp_eps") ? se_get_par_float("bf_amp_eps") : 1.0e-20f;
    const float bf_phase_tol = se_have_par("bf_phase_tol") ? se_get_par_float("bf_phase_tol") : 2.0f;

    sep_t *model = sep_open(in_f, SEP_READ, 0);
    sep_t *out = sep_open(out_f, SEP_WRITE, 0);
    if (model->headers->ndim < 2)
        ERROR(("Need a 2-D velocity model with n1=z and n2=x."));

    const int nz = model->headers->n[0];
    const int nx = model->headers->n[1];
    const float dz = (float)model->headers->d[0];
    const float dx = (float)model->headers->d[1];
    const float oz = (float)model->headers->o[0];
    const float ox = (float)model->headers->o[1];
    if (nz <= 1 || nx <= 1 || dz == 0.0f || dx == 0.0f)
        ERROR(("Invalid model geometry."));
    if (!se_have_par("depth_indices") && depths.empty())
        ERROR(("No valid depths parsed from depths=."));

    std::vector<LayerBoundary> layers;
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
        for (int iz : indices)
        {
            if (iz <= 0 || iz >= nz)
                ERROR(("Requested depth index %d outside valid propagated range 1..%d", iz, nz - 1));
            layers.push_back({iz, oz + iz * dz});
        }
    }
    else
    {
        std::sort(depths.begin(), depths.end());
        depths.erase(std::unique(depths.begin(), depths.end()), depths.end());
        for (float z : depths)
        {
            int iz = (int)std::floor((z - oz) / dz + 0.5f);
            if (iz <= 0 || iz >= nz)
                ERROR(("Requested depth %g maps to iz=%d outside valid propagated range 1..%d", z, iz, nz - 1));
            if (layers.empty() || layers.back().iz != iz)
                layers.push_back({iz, oz + iz * dz});
        }
    }

    std::vector<float> vel((size_t)nz * (size_t)nx);
    se_fsio_read_float(model->data->io, vel.data(), vel.size());

    const float srcx = se_have_par("srcx") ? se_get_par_float("srcx") : ox + 0.5f * (nx - 1) * dx;
    const float srcz = se_have_par("srcz") ? se_get_par_float("srcz") : oz;
    const float vref = se_have_par("vref") ? se_get_par_float("vref") : velocity_at(vel, nz, nx, 0, nx / 2);
    const float omega = 2.0f * (float)M_PI * freq;

    out->headers->ndim = 3;
    out->headers->n[0] = 2;
    out->headers->n[1] = nx;
    out->headers->n[2] = (int)layers.size() + (include_top ? 1 : 0);
    out->headers->d[0] = 1.0;
    out->headers->d[1] = dx;
    out->headers->d[2] = dz;
    out->headers->o[0] = 0.0;
    out->headers->o[1] = ox;
    out->headers->o[2] = include_top ? oz : layers.front().z;

    std::vector<std::complex<float>> u_top(nx), u_bot(nx);
    initialize_hankel_wavefield(nx, ox, dx, srcx, srcz, oz, freq, vref, u_top);

    int iz_top = 0;
    std::vector<float> row((size_t)2 * (size_t)nx);
    if (include_top)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            row[(size_t)ix * 2] = u_top[ix].real();
            row[(size_t)ix * 2 + 1] = u_top[ix].imag();
        }
        se_fsio_write_float(out->data->io, row.data(), row.size());
        INFO(("Wrote top wavefield at depth=%g (iz=0)", oz));
    }
    for (size_t il = 0; il < layers.size(); ++il)
    {
        float **tau = alloc2float(nx, nx);
        build_layer_traveltime(vel, nz, nx, oz, dz, ox, dx, iz_top, layers[il].iz, tau);
        if (use_bf)
            propagate_butterfly(tau, nx, dx, omega, bf_p, bf_leaf, bf_panel_levels, bf_amp_eps, bf_phase_tol, u_top, u_bot);
        else
            propagate_direct(tau, nx, dx, omega, u_top, u_bot);
        free2float(tau);

        for (int ix = 0; ix < nx; ++ix)
        {
            row[(size_t)ix * 2] = u_bot[ix].real();
            row[(size_t)ix * 2 + 1] = u_bot[ix].imag();
        }
        se_fsio_write_float(out->data->io, row.data(), row.size());
        INFO(("Wrote layer %zu wavefield at depth=%g (iz=%d)", il + 1, layers[il].z, layers[il].iz));
        u_top.swap(u_bot);
        iz_top = layers[il].iz;
    }

    sep_close(model);
    sep_close(out);
    se_par_destroy();
    return 0;
}
