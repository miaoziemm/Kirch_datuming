#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_par_sep.h>
#include <SERECKIRCH/include/efmm.h>

#include "../SEWAVE/sewave2d.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
constexpr float PI = 3.14159265358979323846f;

inline size_t index2(int iz, int ix, int nz)
{
    return static_cast<size_t>(ix) * static_cast<size_t>(nz) +
           static_cast<size_t>(iz);
}

inline size_t index3(int it, int iz, int ix, int nt, int nz, int nx)
{
    (void)nt;
    (void)nx;
    return (static_cast<size_t>(it) * static_cast<size_t>(nx) +
            static_cast<size_t>(ix)) *
               static_cast<size_t>(nz) +
           static_cast<size_t>(iz);
}

inline size_t trace_index(int it, int ix, int nx)
{
    return static_cast<size_t>(it) * static_cast<size_t>(nx) +
           static_cast<size_t>(ix);
}

int clamp_int(int value, int lower, int upper)
{
    return std::max(lower, std::min(upper, value));
}

float ricker(float time, float frequency)
{
    const float tau = time - 1.5f / frequency;
    const float value = PI * frequency * tau;
    return (1.0f - 2.0f * value * value) * std::exp(-value * value);
}

struct Block
{
    int index = 0;
    int start = 0;
    int end = 0;
    int first_time = 0;
};

struct PaddedGrid
{
    int physical_nx = 0;
    int physical_nz = 0;
    int nx = 0;
    int nz = 0;
    int nbc = 0;
    float dx = 1.0f;
    float dz = 1.0f;
    float dt = 0.001f;
    std::vector<float> velocity;
    std::vector<float> v2dt2;
    std::vector<float> damping;
};

struct Parameters
{
    int nt = 1000;
    float dt = 0.001f;
    float fdom = 20.0f;
    float sx = 0.0f;
    float sz = 0.0f;
    float rz = 0.0f;
    int ns = 1;
    float ds = 0.0f;

    int blocks = 4;
    int cycles = 3;
    float omega = 1.0f;
    float tol = 1.0e-4f;
    int double_sweep = 1;
    std::string solver = "block_gs";
    std::string initial = "ray";

    int nbc = 40;
    float absorb_alpha = 4.0f;
    float ray_scale = 1.0f;
    float causal_margin = 0.08f;
    int check_interval = 1;
    int compare_reference = 0;
    float max_memory_gb = 16.0f;
};

struct IterationStats
{
    int cycles_done = 0;
    double initial_seconds = 0.0;
    double iteration_seconds = 0.0;
    double residual_seconds = 0.0;
    double reference_seconds = 0.0;
    double relative_residual = 0.0;
    double record_error = -1.0;
    long long updated_points = 0;
};

void print_usage(const char *program)
{
    std::printf(
        "Usage:\n"
        "  %s velocity=vel.rsf output=data.rsf [parameters]\n\n"
        "This executable stores the complete space-time field U(t,x,z) and\n"
        "iterates whole block waveforms rather than storing only interface traces.\n\n"
        "Wave parameters:\n"
        "  nt=INT dt=FLOAT fdom=FLOAT sx=FLOAT sz=FLOAT rz=FLOAT\n"
        "  ns=INT ds=FLOAT\n\n"
        "Space-time iteration parameters:\n"
        "  solver=block_gs|waveform_jacobi   Default: block_gs\n"
        "  initial=ray|zero                  Default: ray\n"
        "  blocks=INT cycles=INT omega=FLOAT tol=FLOAT\n"
        "  double_sweep=0|1                  Downward only or down/up sweep\n"
        "  causal_margin=FLOAT               Seconds before first arrival to update\n"
        "  check_interval=INT                Residual check cadence\n\n"
        "Boundary and diagnostics:\n"
        "  nbc=INT alpha=FLOAT ray_scale=FLOAT\n"
        "  compare_reference=0|1 reference_output=FILE\n"
        "  ray_output=FILE residual_log=FILE snapshot=FILE snapshot_it=INT\n"
        "  max_memory_gb=FLOAT\n\n"
        "Important: block_gs is a waveform-relaxation experiment. It is not\n"
        "expected to beat one direct explicit time march on a serial CPU.\n",
        program != nullptr ? program : "space_time_block_waveform_gs");
}

std::vector<float> fmm_travel_time(
    const sewave::Grid2D &grid,
    float source_x,
    float source_z)
{
    efmm_t context;
    std::memset(&context, 0, sizeof(context));
    const int status = efmm_init(
        &context,
        grid.nz,
        grid.nx,
        grid.dz,
        grid.dx,
        source_z - grid.z0,
        source_x - grid.x0);
    if (status != 0)
        throw std::runtime_error("efmm_init failed");
    if (efmm_set_vel(&context, const_cast<float *>(grid.v.data())) != 0)
    {
        efmm_free(&context);
        throw std::runtime_error("efmm_set_vel failed");
    }
    if (efmm_solver(&context) != 0)
    {
        efmm_free(&context);
        throw std::runtime_error("efmm_solver failed");
    }
    std::vector<float> result(
        context.tt,
        context.tt + static_cast<size_t>(grid.nx) * grid.nz);
    efmm_free(&context);
    return result;
}

PaddedGrid make_padded_grid(
    const sewave::Grid2D &grid,
    const Parameters &parameters)
{
    PaddedGrid padded;
    padded.physical_nx = grid.nx;
    padded.physical_nz = grid.nz;
    padded.nbc = parameters.nbc;
    padded.nx = grid.nx + 2 * parameters.nbc;
    padded.nz = grid.nz + 2 * parameters.nbc;
    padded.dx = grid.dx;
    padded.dz = grid.dz;
    padded.dt = parameters.dt;

    const size_t plane =
        static_cast<size_t>(padded.nx) * static_cast<size_t>(padded.nz);
    padded.velocity.resize(plane);
    padded.v2dt2.resize(plane);
    padded.damping.resize(plane, 1.0f);

    for (int ix = 0; ix < padded.nx; ++ix)
    {
        const int gx = clamp_int(ix - padded.nbc, 0, grid.nx - 1);
        for (int iz = 0; iz < padded.nz; ++iz)
        {
            const int gz = clamp_int(iz - padded.nbc, 0, grid.nz - 1);
            const size_t id = index2(iz, ix, padded.nz);
            const float velocity = grid.v[index2(gz, gx, grid.nz)];
            padded.velocity[id] = velocity;
            padded.v2dt2[id] =
                velocity * velocity * parameters.dt * parameters.dt;

            const int distance = std::min(
                std::min(ix, padded.nx - 1 - ix),
                std::min(iz, padded.nz - 1 - iz));
            if (distance < padded.nbc)
            {
                const float q =
                    (padded.nbc - distance) /
                    static_cast<float>(std::max(1, padded.nbc));
                padded.damping[id] =
                    std::exp(-parameters.absorb_alpha * q * q);
            }
        }
    }
    return padded;
}

std::vector<float> make_padded_travel_time(
    const sewave::Grid2D &grid,
    const PaddedGrid &padded,
    const std::vector<float> &physical_tt)
{
    std::vector<float> tt(
        static_cast<size_t>(padded.nx) * static_cast<size_t>(padded.nz));
    for (int ix = 0; ix < padded.nx; ++ix)
    {
        const int raw_x = ix - padded.nbc;
        const int gx = clamp_int(raw_x, 0, grid.nx - 1);
        const float outside_x =
            static_cast<float>(std::abs(raw_x - gx)) * grid.dx;
        for (int iz = 0; iz < padded.nz; ++iz)
        {
            const int raw_z = iz - padded.nbc;
            const int gz = clamp_int(raw_z, 0, grid.nz - 1);
            const float outside_z =
                static_cast<float>(std::abs(raw_z - gz)) * grid.dz;
            const size_t id = index2(iz, ix, padded.nz);
            const float extra_distance =
                std::sqrt(outside_x * outside_x + outside_z * outside_z);
            tt[id] = physical_tt[index2(gz, gx, grid.nz)] +
                     extra_distance /
                         std::max(padded.velocity[id], 1.0f);
        }
    }
    return tt;
}

std::vector<Block> make_blocks(
    const PaddedGrid &padded,
    const std::vector<float> &travel_time,
    const Parameters &parameters)
{
    std::vector<Block> blocks;
    blocks.reserve(parameters.blocks);
    for (int ib = 0; ib < parameters.blocks; ++ib)
    {
        Block block;
        block.index = ib;
        block.start = static_cast<int>(
            std::llround(
                static_cast<double>(ib) * padded.nz /
                parameters.blocks));
        block.end = static_cast<int>(
            std::llround(
                static_cast<double>(ib + 1) * padded.nz /
                parameters.blocks));
        if (ib == 0)
            block.start = 0;
        if (ib + 1 == parameters.blocks)
            block.end = padded.nz;

        float minimum_time = std::numeric_limits<float>::max();
        for (int ix = 4; ix < padded.nx - 4; ++ix)
        {
            for (int iz = std::max(4, block.start);
                 iz < std::min(padded.nz - 4, block.end);
                 ++iz)
            {
                minimum_time = std::min(
                    minimum_time,
                    travel_time[index2(iz, ix, padded.nz)]);
            }
        }
        if (!std::isfinite(minimum_time))
            minimum_time = 0.0f;
        block.first_time = std::max(
            0,
            static_cast<int>(std::floor(
                (minimum_time - parameters.causal_margin) /
                parameters.dt)));
        block.first_time = std::min(block.first_time, parameters.nt - 1);
        blocks.push_back(block);
    }
    return blocks;
}

void estimate_and_check_memory(
    const PaddedGrid &padded,
    const Parameters &parameters)
{
    const long double points =
        static_cast<long double>(parameters.nt) * padded.nx * padded.nz;
    const int full_fields =
        parameters.solver == "waveform_jacobi" ? 2 : 1;
    const long double bytes =
        points * sizeof(float) * full_fields +
        static_cast<long double>(padded.nx) * padded.nz * sizeof(float) * 6;
    const long double gib = bytes / (1024.0L * 1024.0L * 1024.0L);
    INFO((
        "estimated main memory for space-time solver: %.3Lf GiB (%d full fields)",
        gib,
        full_fields));
    if (gib > parameters.max_memory_gb)
        throw std::runtime_error(
            "estimated space-time storage exceeds max_memory_gb");
}

void initialize_ray_field(
    const PaddedGrid &padded,
    const std::vector<float> &travel_time,
    int nt,
    float dt,
    float frequency,
    float source_velocity,
    float ray_scale,
    std::vector<float> &u)
{
    const size_t total =
        static_cast<size_t>(nt) * padded.nx * padded.nz;
    u.assign(total, 0.0f);
    const float spacing = std::min(padded.dx, padded.dz);
    const int half_window = std::max(
        8,
        static_cast<int>(std::ceil(3.5f / (frequency * dt))));

    for (int ix = 0; ix < padded.nx; ++ix)
    {
        for (int iz = 0; iz < padded.nz; ++iz)
        {
            const size_t id2 = index2(iz, ix, padded.nz);
            const float tt = travel_time[id2];
            const float distance =
                std::max(source_velocity * tt, spacing);
            const float amplitude =
                ray_scale * std::sqrt(spacing / distance);
            const int center = static_cast<int>(std::lround(
                (tt + 1.5f / frequency) / dt));
            const int begin = std::max(1, center - half_window);
            const int end = std::min(nt, center + half_window + 1);
            for (int it = begin; it < end; ++it)
            {
                u[index3(it, iz, ix, nt, padded.nz, padded.nx)] =
                    padded.damping[id2] * amplitude *
                    ricker(it * dt - tt, frequency);
            }
        }
    }

    // The PDE solve uses zero initial displacement.
    std::fill(
        u.begin(),
        u.begin() + static_cast<size_t>(padded.nx) * padded.nz,
        0.0f);
}

void initialize_zero_field(
    const PaddedGrid &padded,
    int nt,
    std::vector<float> &u)
{
    u.assign(
        static_cast<size_t>(nt) * padded.nx * padded.nz,
        0.0f);
}

float laplacian_fd8(
    const float *snapshot,
    int iz,
    int ix,
    const PaddedGrid &padded)
{
    const float d2x =
        (-1.0f / 560.0f) *
            (snapshot[index2(iz, ix - 4, padded.nz)] +
             snapshot[index2(iz, ix + 4, padded.nz)]) +
        (8.0f / 315.0f) *
            (snapshot[index2(iz, ix - 3, padded.nz)] +
             snapshot[index2(iz, ix + 3, padded.nz)]) -
        (1.0f / 5.0f) *
            (snapshot[index2(iz, ix - 2, padded.nz)] +
             snapshot[index2(iz, ix + 2, padded.nz)]) +
        (8.0f / 5.0f) *
            (snapshot[index2(iz, ix - 1, padded.nz)] +
             snapshot[index2(iz, ix + 1, padded.nz)]) -
        (205.0f / 72.0f) * snapshot[index2(iz, ix, padded.nz)];

    const float d2z =
        (-1.0f / 560.0f) *
            (snapshot[index2(iz - 4, ix, padded.nz)] +
             snapshot[index2(iz + 4, ix, padded.nz)]) +
        (8.0f / 315.0f) *
            (snapshot[index2(iz - 3, ix, padded.nz)] +
             snapshot[index2(iz + 3, ix, padded.nz)]) -
        (1.0f / 5.0f) *
            (snapshot[index2(iz - 2, ix, padded.nz)] +
             snapshot[index2(iz + 2, ix, padded.nz)]) +
        (8.0f / 5.0f) *
            (snapshot[index2(iz - 1, ix, padded.nz)] +
             snapshot[index2(iz + 1, ix, padded.nz)]) -
        (205.0f / 72.0f) * snapshot[index2(iz, ix, padded.nz)];

    return d2x / (padded.dx * padded.dx) +
           d2z / (padded.dz * padded.dz);
}

inline float recurrence_candidate(
    const std::vector<float> &u,
    int it,
    int iz,
    int ix,
    int source_iz,
    int source_ix,
    const std::vector<float> &source,
    const PaddedGrid &padded)
{
    const size_t plane =
        static_cast<size_t>(padded.nx) * padded.nz;
    const float *current =
        u.data() + static_cast<size_t>(it) * plane;
    const float *previous =
        it == 0
            ? nullptr
            : u.data() + static_cast<size_t>(it - 1) * plane;
    const size_t id2 = index2(iz, ix, padded.nz);
    const float current_value = current[id2];
    const float previous_value = previous == nullptr ? 0.0f : previous[id2];
    float value =
        2.0f * current_value -
        padded.damping[id2] * previous_value +
        padded.v2dt2[id2] *
            laplacian_fd8(current, iz, ix, padded);
    if (iz == source_iz && ix == source_ix)
        value += source[it];
    return padded.damping[id2] * value;
}

long long solve_block_waveform_inplace(
    std::vector<float> &u,
    const Block &block,
    const PaddedGrid &padded,
    const Parameters &parameters,
    const std::vector<float> &source,
    int source_iz,
    int source_ix)
{
    const int z_begin = std::max(4, block.start);
    const int z_end = std::min(padded.nz - 4, block.end);
    if (z_begin >= z_end)
        return 0;

    long long updated = 0;
    for (int it = block.first_time; it < parameters.nt - 1; ++it)
    {
        const size_t next_offset =
            static_cast<size_t>(it + 1) * padded.nx * padded.nz;
        for (int ix = 4; ix < padded.nx - 4; ++ix)
        {
            for (int iz = z_begin; iz < z_end; ++iz)
            {
                const float candidate = recurrence_candidate(
                    u,
                    it,
                    iz,
                    ix,
                    source_iz,
                    source_ix,
                    source,
                    padded);
                const size_t id =
                    next_offset + index2(iz, ix, padded.nz);
                u[id] += parameters.omega * (candidate - u[id]);
                ++updated;
            }
        }
    }
    return updated;
}

long long waveform_jacobi_iteration(
    const std::vector<float> &old_u,
    std::vector<float> &new_u,
    const std::vector<Block> &blocks,
    const PaddedGrid &padded,
    const Parameters &parameters,
    const std::vector<float> &source,
    int source_iz,
    int source_ix)
{
    new_u = old_u;
    long long updated = 0;
    for (const Block &block : blocks)
    {
        const int z_begin = std::max(4, block.start);
        const int z_end = std::min(padded.nz - 4, block.end);
        for (int it = block.first_time; it < parameters.nt - 1; ++it)
        {
            const size_t next_offset =
                static_cast<size_t>(it + 1) * padded.nx * padded.nz;
            for (int ix = 4; ix < padded.nx - 4; ++ix)
            {
                for (int iz = z_begin; iz < z_end; ++iz)
                {
                    const float candidate = recurrence_candidate(
                        old_u,
                        it,
                        iz,
                        ix,
                        source_iz,
                        source_ix,
                        source,
                        padded);
                    const size_t id =
                        next_offset + index2(iz, ix, padded.nz);
                    new_u[id] =
                        old_u[id] + parameters.omega *
                                        (candidate - old_u[id]);
                    ++updated;
                }
            }
        }
    }
    return updated;
}

// This is the true all-at-once recurrence residual:
// r^{n+1} = F(u^n,u^{n-1}) - u^{n+1}.
double relative_recurrence_residual(
    const std::vector<float> &u,
    const PaddedGrid &padded,
    const Parameters &parameters,
    const std::vector<float> &source,
    int source_iz,
    int source_ix)
{
    long double numerator = 0.0L;
    long double denominator = 0.0L;
    for (int it = 0; it < parameters.nt - 1; ++it)
    {
        const size_t next_offset =
            static_cast<size_t>(it + 1) * padded.nx * padded.nz;
        for (int ix = 4; ix < padded.nx - 4; ++ix)
        {
            for (int iz = 4; iz < padded.nz - 4; ++iz)
            {
                const float candidate = recurrence_candidate(
                    u,
                    it,
                    iz,
                    ix,
                    source_iz,
                    source_ix,
                    source,
                    padded);
                const float actual =
                    u[next_offset + index2(iz, ix, padded.nz)];
                const long double residual =
                    static_cast<long double>(candidate) - actual;
                numerator += residual * residual;
                denominator +=
                    static_cast<long double>(candidate) * candidate;
            }
        }
    }
    return std::sqrt(
        static_cast<double>(
            numerator /
            std::max(denominator, std::numeric_limits<long double>::min())));
}

std::vector<float> extract_record(
    const std::vector<float> &u,
    const PaddedGrid &padded,
    int nt,
    int receiver_physical_z)
{
    std::vector<float> record(
        static_cast<size_t>(nt) * padded.physical_nx,
        0.0f);
    const int iz = padded.nbc + receiver_physical_z;
    for (int it = 0; it < nt; ++it)
    {
        for (int ix = 0; ix < padded.physical_nx; ++ix)
        {
            record[trace_index(it, ix, padded.physical_nx)] =
                u[index3(
                    it,
                    iz,
                    padded.nbc + ix,
                    nt,
                    padded.nz,
                    padded.nx)];
        }
    }
    return record;
}

std::vector<float> direct_reference(
    const PaddedGrid &padded,
    const Parameters &parameters,
    const std::vector<float> &source,
    int source_iz,
    int source_ix,
    int receiver_physical_z)
{
    const size_t plane =
        static_cast<size_t>(padded.nx) * padded.nz;
    std::vector<float> previous(plane, 0.0f);
    std::vector<float> current(plane, 0.0f);
    std::vector<float> following(plane, 0.0f);
    std::vector<float> record(
        static_cast<size_t>(parameters.nt) * padded.physical_nx,
        0.0f);

    const int receiver_iz = padded.nbc + receiver_physical_z;
    for (int it = 0; it < parameters.nt - 1; ++it)
    {
        std::fill(following.begin(), following.end(), 0.0f);
        for (int ix = 4; ix < padded.nx - 4; ++ix)
        {
            for (int iz = 4; iz < padded.nz - 4; ++iz)
            {
                const size_t id = index2(iz, ix, padded.nz);
                float value =
                    2.0f * current[id] -
                    padded.damping[id] * previous[id] +
                    padded.v2dt2[id] *
                        laplacian_fd8(current.data(), iz, ix, padded);
                if (iz == source_iz && ix == source_ix)
                    value += source[it];
                following[id] = padded.damping[id] * value;
            }
        }
        previous.swap(current);
        current.swap(following);
        for (int ix = 0; ix < padded.physical_nx; ++ix)
        {
            record[trace_index(it + 1, ix, padded.physical_nx)] =
                current[index2(
                    receiver_iz,
                    padded.nbc + ix,
                    padded.nz)];
        }
    }
    return record;
}

double relative_error(
    const std::vector<float> &a,
    const std::vector<float> &b)
{
    if (a.size() != b.size())
        throw std::runtime_error("record sizes differ");
    long double numerator = 0.0L;
    long double denominator = 0.0L;
    for (size_t i = 0; i < a.size(); ++i)
    {
        const long double difference =
            static_cast<long double>(a[i]) - b[i];
        numerator += difference * difference;
        denominator += static_cast<long double>(b[i]) * b[i];
    }
    return std::sqrt(
        static_cast<double>(
            numerator /
            std::max(denominator, std::numeric_limits<long double>::min())));
}

void write_record(
    const char *path,
    const std::vector<float> &record,
    const sewave::Grid2D &grid,
    const Parameters &parameters,
    int shot_index)
{
    (void)shot_index;
    sewave::Data3D output;
    output.nt = parameters.nt;
    output.nr = grid.nx;
    output.ns = 1;
    output.dt = parameters.dt;
    output.dr = grid.dx;
    output.ds = parameters.ds;
    output.t0 = 0.0f;
    output.r0 = grid.x0;
    output.s0 = parameters.sx;
    output.d.assign(
        static_cast<size_t>(output.nt) * output.nr,
        0.0f);
    for (int ix = 0; ix < output.nr; ++ix)
    {
        for (int it = 0; it < output.nt; ++it)
        {
            output.d[static_cast<size_t>(ix) * output.nt + it] =
                record[trace_index(it, ix, output.nr)];
        }
    }
    sewave::write_rsf3d(path, output);
}

void write_snapshot(
    const char *path,
    const std::vector<float> &u,
    const sewave::Grid2D &grid,
    const PaddedGrid &padded,
    int nt,
    int snapshot_it)
{
    const int it = clamp_int(snapshot_it, 0, nt - 1);
    std::vector<float> snapshot(
        static_cast<size_t>(grid.nx) * grid.nz,
        0.0f);
    for (int ix = 0; ix < grid.nx; ++ix)
    {
        for (int iz = 0; iz < grid.nz; ++iz)
        {
            snapshot[index2(iz, ix, grid.nz)] =
                u[index3(
                    it,
                    padded.nbc + iz,
                    padded.nbc + ix,
                    nt,
                    padded.nz,
                    padded.nx)];
        }
    }
    sewave::write_rsf2d(path, grid, snapshot, "space-time iterative snapshot");
}

int coordinate_index(float coordinate, float origin, float spacing)
{
    return static_cast<int>(std::lround((coordinate - origin) / spacing));
}

void validate(
    const sewave::Grid2D &grid,
    const Parameters &parameters,
    int source_z,
    int source_x,
    int receiver_z)
{
    if (parameters.nt < 3 || parameters.dt <= 0.0f ||
        parameters.fdom <= 0.0f)
        throw std::runtime_error("invalid nt, dt, or fdom");
    if (parameters.blocks < 1)
        throw std::runtime_error("blocks must be positive");
    if (parameters.cycles < 1)
        throw std::runtime_error("cycles must be positive");
    if (!(parameters.omega > 0.0f && parameters.omega <= 1.5f))
        throw std::runtime_error("require 0 < omega <= 1.5");
    if (parameters.nbc < 8)
        throw std::runtime_error("nbc must be at least 8 for FD8");
    if (parameters.solver != "block_gs" &&
        parameters.solver != "waveform_jacobi")
        throw std::runtime_error(
            "solver must be block_gs or waveform_jacobi");
    if (parameters.initial != "ray" && parameters.initial != "zero")
        throw std::runtime_error("initial must be ray or zero");
    if (source_z < 0 || source_z >= grid.nz ||
        source_x < 0 || source_x >= grid.nx ||
        receiver_z < 0 || receiver_z >= grid.nz)
        throw std::runtime_error("source or receiver lies outside model");

    const float maximum_velocity =
        *std::max_element(grid.v.begin(), grid.v.end());
    const double cfl =
        maximum_velocity * parameters.dt *
        std::sqrt(
            1.0 / (grid.dx * grid.dx) +
            1.0 / (grid.dz * grid.dz));
    if (cfl > 0.45)
        INFO(("WARNING: conservative FD8 CFL limit exceeded: %g", cfl));
}

} // namespace

int main(int argc, char **argv)
{
    if (argc <= 1 ||
        (argc == 2 &&
         (!std::strcmp(argv[1], "-h") ||
          !std::strcmp(argv[1], "--help") ||
          !std::strcmp(argv[1], "help"))))
    {
        print_usage(argv[0]);
        return 0;
    }

    se_par_init(argc, argv);
    try
    {
        if (!se_have_par("velocity"))
            throw std::runtime_error("Need velocity= velocity RSF");
        if (!se_have_par("output"))
            throw std::runtime_error("Need output= output RSF");

        const char *velocity_path = se_get_par_str("velocity");
        const char *output_path = se_get_par_str("output");
        const sewave::Grid2D grid = sewave::read_rsf2d(velocity_path);

        Parameters p;
        p.nt = se_have_par("nt") ? se_get_par_int("nt") : p.nt;
        p.dt = se_have_par("dt") ? se_get_par_float("dt") : p.dt;
        p.fdom = se_have_par("fdom") ? se_get_par_float("fdom") : p.fdom;
        p.sx = se_have_par("sx")
                   ? se_get_par_float("sx")
                   : grid.x0 + 0.5f * (grid.nx - 1) * grid.dx;
        p.sz = se_have_par("sz") ? se_get_par_float("sz") : grid.z0;
        p.rz = se_have_par("rz") ? se_get_par_float("rz") : p.sz;
        p.ns = se_have_par("ns") ? se_get_par_int("ns") : p.ns;
        p.ds = se_have_par("ds") ? se_get_par_float("ds") : p.ds;
        p.blocks = se_have_par("blocks")
                       ? se_get_par_int("blocks")
                       : p.blocks;
        p.cycles = se_have_par("cycles")
                       ? se_get_par_int("cycles")
                       : p.cycles;
        p.omega = se_have_par("omega")
                      ? se_get_par_float("omega")
                      : p.omega;
        p.tol = se_have_par("tol") ? se_get_par_float("tol") : p.tol;
        p.double_sweep = se_have_par("double_sweep")
                             ? se_get_par_int("double_sweep")
                             : p.double_sweep;
        p.solver = se_have_par("solver")
                       ? se_get_par_str("solver")
                       : p.solver;
        p.initial = se_have_par("initial")
                        ? se_get_par_str("initial")
                        : p.initial;
        p.nbc = se_have_par("nbc") ? se_get_par_int("nbc") : p.nbc;
        p.absorb_alpha = se_have_par("alpha")
                             ? se_get_par_float("alpha")
                             : p.absorb_alpha;
        p.ray_scale = se_have_par("ray_scale")
                          ? se_get_par_float("ray_scale")
                          : p.ray_scale;
        p.causal_margin = se_have_par("causal_margin")
                              ? se_get_par_float("causal_margin")
                              : p.causal_margin;
        p.check_interval = se_have_par("check_interval")
                               ? se_get_par_int("check_interval")
                               : p.check_interval;
        p.compare_reference = se_have_par("compare_reference")
                                  ? se_get_par_int("compare_reference")
                                  : p.compare_reference;
        p.max_memory_gb = se_have_par("max_memory_gb")
                              ? se_get_par_float("max_memory_gb")
                              : p.max_memory_gb;

        if (p.ns != 1)
            throw std::runtime_error(
                "This research prototype currently supports ns=1 only");

        const int source_x = coordinate_index(p.sx, grid.x0, grid.dx);
        const int source_z = coordinate_index(p.sz, grid.z0, grid.dz);
        const int receiver_z = coordinate_index(p.rz, grid.z0, grid.dz);
        validate(grid, p, source_z, source_x, receiver_z);

        const PaddedGrid padded = make_padded_grid(grid, p);
        estimate_and_check_memory(padded, p);

        const auto total_begin = std::chrono::steady_clock::now();
        const std::vector<float> physical_tt =
            fmm_travel_time(grid, p.sx, p.sz);
        const std::vector<float> padded_tt =
            make_padded_travel_time(grid, padded, physical_tt);
        const std::vector<Block> blocks = make_blocks(padded, padded_tt, p);

        for (const Block &block : blocks)
        {
            INFO((
                "block %d/%d: z=[%d,%d), first update sample=%d",
                block.index + 1,
                p.blocks,
                block.start,
                block.end,
                block.first_time));
        }

        std::vector<float> source(static_cast<size_t>(p.nt));
        for (int it = 0; it < p.nt; ++it)
            source[it] = ricker(it * p.dt, p.fdom);

        const int source_padded_x = padded.nbc + source_x;
        const int source_padded_z = padded.nbc + source_z;
        const float source_velocity =
            grid.v[index2(source_z, source_x, grid.nz)];

        IterationStats stats;
        std::vector<float> u;
        const auto initial_begin = std::chrono::steady_clock::now();
        if (p.initial == "ray")
        {
            initialize_ray_field(
                padded,
                padded_tt,
                p.nt,
                p.dt,
                p.fdom,
                source_velocity,
                p.ray_scale,
                u);
        }
        else
        {
            initialize_zero_field(padded, p.nt, u);
        }
        stats.initial_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - initial_begin)
                .count();

        const std::vector<float> initial_record =
            extract_record(u, padded, p.nt, receiver_z);
        if (se_have_par("ray_output"))
            write_record(
                se_get_par_str("ray_output"),
                initial_record,
                grid,
                p,
                0);

        FILE *residual_file = nullptr;
        if (se_have_par("residual_log"))
        {
            residual_file = std::fopen(se_get_par_str("residual_log"), "w");
            if (residual_file == nullptr)
                throw std::runtime_error("cannot open residual_log");
            std::fprintf(
                residual_file,
                "cycle,relative_recurrence_residual,updated_points,elapsed_seconds\n");
        }

        std::vector<float> new_u;
        if (p.solver == "waveform_jacobi")
            new_u.resize(u.size());

        const auto iteration_begin = std::chrono::steady_clock::now();
        for (int cycle = 0; cycle < p.cycles; ++cycle)
        {
            const auto cycle_begin = std::chrono::steady_clock::now();
            long long cycle_updates = 0;

            if (p.solver == "block_gs")
            {
                for (const Block &block : blocks)
                {
                    cycle_updates += solve_block_waveform_inplace(
                        u,
                        block,
                        padded,
                        p,
                        source,
                        source_padded_z,
                        source_padded_x);
                }
                if (p.double_sweep)
                {
                    for (auto iterator = blocks.rbegin();
                         iterator != blocks.rend();
                         ++iterator)
                    {
                        cycle_updates += solve_block_waveform_inplace(
                            u,
                            *iterator,
                            padded,
                            p,
                            source,
                            source_padded_z,
                            source_padded_x);
                    }
                }
            }
            else
            {
                cycle_updates = waveform_jacobi_iteration(
                    u,
                    new_u,
                    blocks,
                    padded,
                    p,
                    source,
                    source_padded_z,
                    source_padded_x);
                u.swap(new_u);
            }

            stats.updated_points += cycle_updates;
            stats.cycles_done = cycle + 1;
            const double cycle_seconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - cycle_begin)
                    .count();

            const bool need_check =
                p.check_interval > 0 &&
                ((cycle + 1) % p.check_interval == 0 ||
                 cycle + 1 == p.cycles);
            if (need_check)
            {
                const auto residual_begin =
                    std::chrono::steady_clock::now();
                stats.relative_residual = relative_recurrence_residual(
                    u,
                    padded,
                    p,
                    source,
                    source_padded_z,
                    source_padded_x);
                stats.residual_seconds +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - residual_begin)
                        .count();
                INFO((
                    "cycle %d/%d: recurrence residual=%e updates=%lld time=%.6f s",
                    cycle + 1,
                    p.cycles,
                    stats.relative_residual,
                    cycle_updates,
                    cycle_seconds));
                if (residual_file != nullptr)
                {
                    std::fprintf(
                        residual_file,
                        "%d,%.12e,%lld,%.9f\n",
                        cycle + 1,
                        stats.relative_residual,
                        cycle_updates,
                        cycle_seconds);
                    std::fflush(residual_file);
                }
                if (stats.relative_residual <= p.tol)
                    break;
            }
            else
            {
                INFO((
                    "cycle %d/%d: updates=%lld time=%.6f s",
                    cycle + 1,
                    p.cycles,
                    cycle_updates,
                    cycle_seconds));
            }
        }
        stats.iteration_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - iteration_begin)
                .count();
        if (residual_file != nullptr)
            std::fclose(residual_file);

        const std::vector<float> record =
            extract_record(u, padded, p.nt, receiver_z);
        write_record(output_path, record, grid, p, 0);

        if (se_have_par("snapshot"))
        {
            const int snapshot_it = se_have_par("snapshot_it")
                                        ? se_get_par_int("snapshot_it")
                                        : p.nt - 1;
            write_snapshot(
                se_get_par_str("snapshot"),
                u,
                grid,
                padded,
                p.nt,
                snapshot_it);
        }

        if (p.compare_reference || se_have_par("reference_output"))
        {
            const auto reference_begin = std::chrono::steady_clock::now();
            const std::vector<float> reference = direct_reference(
                padded,
                p,
                source,
                source_padded_z,
                source_padded_x,
                receiver_z);
            stats.reference_seconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - reference_begin)
                    .count();
            stats.record_error = relative_error(record, reference);
            INFO((
                "reference: time=%.6f s, relative record error=%e",
                stats.reference_seconds,
                stats.record_error));
            if (se_have_par("reference_output"))
                write_record(
                    se_get_par_str("reference_output"),
                    reference,
                    grid,
                    p,
                    0);
        }

        const double total_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - total_begin)
                .count();
        INFO((
            "space-time block summary: solver=%s initial=%s cycles=%d "
            "initial=%.6f s iteration=%.6f s residual_checks=%.6f s "
            "reference=%.6f s total=%.6f s updated_points=%lld "
            "final_residual=%e record_error=%e",
            p.solver.c_str(),
            p.initial.c_str(),
            stats.cycles_done,
            stats.initial_seconds,
            stats.iteration_seconds,
            stats.residual_seconds,
            stats.reference_seconds,
            total_seconds,
            stats.updated_points,
            stats.relative_residual,
            stats.record_error));

        se_par_destroy();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::fprintf(
            stderr,
            "space_time_block_waveform_gs: %s\n",
            error.what());
        se_par_destroy();
        return 1;
    }
}
