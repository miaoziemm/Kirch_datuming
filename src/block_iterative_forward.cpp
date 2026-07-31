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
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr float PI = 3.14159265358979323846f;

inline size_t index2(int iz, int ix, int nz)
{
    return static_cast<size_t>(ix) * static_cast<size_t>(nz) +
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
    int core_start = 0;
    int core_end = 0;
    int local_start = 0;
    int local_end = 0;
};

struct AbsorbPoint
{
    size_t id = 0;
    size_t id1 = 0;
    size_t id2 = 0;
    size_t id3 = 0;
    size_t id4 = 0;
    float weight = 0.0f;
    float t1 = 0.0f;
    float t2 = 0.0f;
    float t3 = 0.0f;
};

struct Domain
{
    Block block;
    int physical_nx = 0;
    int nx = 0;
    int nz = 0;
    int nbc = 0;
    int top_pad = 0;
    int bottom_pad = 0;
    bool physical_top = false;
    bool physical_bottom = false;
    float dx = 1.0f;
    float dz = 1.0f;
    float dt = 0.001f;
    int absorb_length = 30;
    float absorb_alpha = 1.0f;
    std::vector<float> velocity;
    std::vector<float> v2dt2;
    std::vector<AbsorbPoint> absorb_points;
};

struct State
{
    std::map<int, std::vector<float>> traces;
    std::vector<float> receiver;
};

struct SolveStats
{
    long long updated_points = 0;
    int active_steps = 0;
};

struct ShotStats
{
    double fmm_seconds = 0.0;
    double initial_seconds = 0.0;
    double sweep_seconds = 0.0;
    double ray_scale = 1.0;
    double final_interface_residual = 0.0;
    long long updated_points = 0;
    int local_solves = 0;
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
    int overlap = 22;
    int cycles = 1;
    float beta = 1.0f;
    float gate = 0.003f;
    int nbc = 40;
    int absorb_length = 30;
    float absorb_alpha = 1.0f;
    float ray_scale = 1.0f;
    int progress_interval = 0;
};

void print_usage(const char *program)
{
    const char *name =
        (program != nullptr && program[0] != '\0')
            ? program
            : "block_iterative_forward";
    std::printf(
        "Usage:\n"
        "  %s velocity=vel.rsf output=data.rsf [parameters]\n\n"
        "Wave and geometry parameters:\n"
        "  nt=INT dt=FLOAT fdom=FLOAT  Time samples, step, and Ricker dominant frequency.\n"
        "  sx=FLOAT sz=FLOAT rz=FLOAT  First source x, source depth, and receiver depth.\n"
        "  ns=INT ds=FLOAT             Number of shots and shot spacing.\n\n"
        "Block-iteration parameters:\n"
        "  blocks=INT                  Number of equal-depth core blocks. Default: 4.\n"
        "  overlap=INT                 Per-side overlap in grid rows. Default: 22.\n"
        "  cycles=INT                  Down/up Robin sweep cycles. Default: 1.\n"
        "  beta=FLOAT                  Characteristic/Robin coefficient. Default: 1.0.\n"
        "  gate=FLOAT                  Relative causal activation gate. Default: 0.003.\n"
        "  ray_scale=FLOAT             FMM geometrical-spreading amplitude scale. Default: 1.0.\n\n"
        "Boundary and diagnostic parameters:\n"
        "  nbc=INT L=INT alpha=FLOAT   Side/physical-boundary padding and absorber settings.\n"
        "  ray_output=FILE             Optional initial FMM ray record in RSF format.\n"
        "  traveltime=FILE             Optional FMM travel-time field (one shot only).\n"
        "  progress_interval=INT       Print local time-step progress every INT samples.\n\n"
        "The program is serial. It uses the existing efmm implementation for the\n"
        "first-arrival field and the repository RSF readers/writers for I/O.\n",
        name);
}

std::vector<Block> make_blocks(int nz, int count, int overlap)
{
    std::vector<int> edges(static_cast<size_t>(count) + 1);
    for (int ib = 0; ib <= count; ++ib)
    {
        const double position =
            static_cast<double>(ib) * static_cast<double>(nz) /
            static_cast<double>(count);
        edges[ib] = static_cast<int>(std::floor(position + 0.5));
    }
    edges.front() = 0;
    edges.back() = nz;

    std::vector<Block> blocks;
    blocks.reserve(count);
    for (int ib = 0; ib < count; ++ib)
    {
        Block block;
        block.index = ib;
        block.core_start = edges[ib];
        block.core_end = edges[ib + 1];
        block.local_start = std::max(0, block.core_start - overlap);
        block.local_end = std::min(nz, block.core_end + overlap);
        blocks.push_back(block);
    }
    return blocks;
}

std::vector<int> robin_query_rows(
    const std::vector<Block> &blocks,
    int block_index)
{
    std::vector<int> rows;
    auto append_three = [&rows](int row)
    {
        rows.push_back(row - 1);
        rows.push_back(row);
        rows.push_back(row + 1);
    };

    if (block_index > 0)
    {
        append_three(blocks[block_index].local_start);
        append_three(blocks[block_index - 1].local_end - 1);
    }
    if (block_index + 1 < static_cast<int>(blocks.size()))
    {
        append_three(blocks[block_index + 1].local_start);
        append_three(blocks[block_index].local_end - 1);
    }

    const Block &block = blocks[block_index];
    rows.erase(
        std::remove_if(
            rows.begin(),
            rows.end(),
            [&block](int row)
            {
                return row < block.local_start || row >= block.local_end;
            }),
        rows.end());
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    return rows;
}

Domain make_domain(
    const sewave::Grid2D &grid,
    const Block &block,
    const Parameters &parameters)
{
    Domain domain;
    domain.block = block;
    domain.physical_nx = grid.nx;
    domain.nbc = parameters.nbc;
    domain.physical_top = block.local_start == 0;
    domain.physical_bottom = block.local_end == grid.nz;
    domain.top_pad = domain.physical_top ? parameters.nbc : 0;
    domain.bottom_pad = domain.physical_bottom ? parameters.nbc : 0;
    domain.nx = grid.nx + 2 * parameters.nbc;
    domain.nz =
        block.local_end - block.local_start +
        domain.top_pad + domain.bottom_pad;
    domain.dx = grid.dx;
    domain.dz = grid.dz;
    domain.dt = parameters.dt;
    domain.absorb_length = parameters.absorb_length;
    domain.absorb_alpha = parameters.absorb_alpha;
    domain.velocity.resize(
        static_cast<size_t>(domain.nx) * static_cast<size_t>(domain.nz));
    domain.v2dt2.resize(domain.velocity.size());

    for (int ix = 0; ix < domain.nx; ++ix)
    {
        const int global_x =
            clamp_int(ix - parameters.nbc, 0, grid.nx - 1);
        for (int iz = 0; iz < domain.nz; ++iz)
        {
            const int global_z = clamp_int(
                block.local_start + iz - domain.top_pad,
                0,
                grid.nz - 1);
            const float velocity =
                grid.v[index2(global_z, global_x, grid.nz)];
            const size_t id = index2(iz, ix, domain.nz);
            domain.velocity[id] = velocity;
            domain.v2dt2[id] =
                velocity * velocity * parameters.dt * parameters.dt;
        }
    }

    const int length = std::min(
        domain.absorb_length,
        std::min(domain.nx, domain.nz) / 2);
    if (length >= 5)
    {
        auto append = [&](int ix, int iz, int step_x, int step_z, float weight)
        {
            AbsorbPoint point;
            point.id = index2(iz, ix, domain.nz);
            point.id1 = index2(
                iz + step_z,
                ix + step_x,
                domain.nz);
            point.id2 = index2(
                iz + 2 * step_z,
                ix + 2 * step_x,
                domain.nz);
            point.id3 = index2(
                iz + 3 * step_z,
                ix + 3 * step_x,
                domain.nz);
            point.id4 = index2(
                iz + 4 * step_z,
                ix + 4 * step_x,
                domain.nz);
            point.weight = weight;
            const float sigma =
                domain.absorb_alpha *
                domain.velocity[point.id] *
                domain.dt /
                std::max(domain.dx, domain.dz);
            point.t1 = (2.0f - sigma) * (1.0f - sigma) / 2.0f;
            point.t2 = sigma * (2.0f - sigma);
            point.t3 = sigma * (sigma - 1.0f) / 2.0f;
            domain.absorb_points.push_back(point);
        };

        domain.absorb_points.reserve(
            static_cast<size_t>(2 * (length - 4)) *
            static_cast<size_t>(domain.nx + domain.nz));
        for (int ix = 0; ix < length - 4; ++ix)
        {
            const float weight =
                1.0f - ix / static_cast<float>(length);
            for (int iz = 0; iz < domain.nz; ++iz)
                append(ix, iz, 1, 0, weight);
        }
        for (int ix = domain.nx - length + 4;
             ix < domain.nx;
             ++ix)
        {
            const float weight =
                1.0f -
                (domain.nx - ix - 1) /
                    static_cast<float>(length);
            for (int iz = 0; iz < domain.nz; ++iz)
                append(ix, iz, -1, 0, weight);
        }
        if (domain.physical_top)
        {
            for (int ix = 0; ix < domain.nx; ++ix)
            {
                for (int iz = 0; iz < length - 4; ++iz)
                {
                    append(
                        ix,
                        iz,
                        0,
                        1,
                        1.0f -
                            iz / static_cast<float>(length));
                }
            }
        }
        if (domain.physical_bottom)
        {
            for (int ix = 0; ix < domain.nx; ++ix)
            {
                for (int iz = domain.nz - length + 4;
                     iz < domain.nz;
                     ++iz)
                {
                    append(
                        ix,
                        iz,
                        0,
                        -1,
                        1.0f -
                            (domain.nz - iz - 1) /
                                static_cast<float>(length));
                }
            }
        }
    }
    return domain;
}

void advance_wavefield(
    const Domain &domain,
    const std::vector<float> &previous,
    const std::vector<float> &current,
    std::vector<float> &following)
{
    const float dx2 = domain.dx * domain.dx;
    const float dz2 = domain.dz * domain.dz;

    // The old implementation first evaluated FD2 over almost the complete
    // domain, overwrote it with FD8 in the interior, and then made a third
    // pass to advance the wavefield.  Initialize the zero-Laplacian boundary
    // values once and directly overwrite only the points that need FD8/FD2.
    for (size_t id = 0; id < following.size(); ++id)
        following[id] = 2.0f * current[id] - previous[id];

    for (int ix = 4; ix < domain.nx - 4; ++ix)
    {
        for (int iz = 4; iz < domain.nz - 4; ++iz)
        {
            const size_t id = index2(iz, ix, domain.nz);
            const float d2x =
                (-1.0f / 560.0f) *
                    (current[index2(iz, ix - 4, domain.nz)] +
                     current[index2(iz, ix + 4, domain.nz)]) +
                (8.0f / 315.0f) *
                    (current[index2(iz, ix - 3, domain.nz)] +
                     current[index2(iz, ix + 3, domain.nz)]) -
                (1.0f / 5.0f) *
                    (current[index2(iz, ix - 2, domain.nz)] +
                     current[index2(iz, ix + 2, domain.nz)]) +
                (8.0f / 5.0f) *
                    (current[index2(iz, ix - 1, domain.nz)] +
                     current[index2(iz, ix + 1, domain.nz)]) -
                (205.0f / 72.0f) * current[id];
            const float d2z =
                (-1.0f / 560.0f) *
                    (current[index2(iz - 4, ix, domain.nz)] +
                     current[index2(iz + 4, ix, domain.nz)]) +
                (8.0f / 315.0f) *
                    (current[index2(iz - 3, ix, domain.nz)] +
                     current[index2(iz + 3, ix, domain.nz)]) -
                (1.0f / 5.0f) *
                    (current[index2(iz - 2, ix, domain.nz)] +
                     current[index2(iz + 2, ix, domain.nz)]) +
                (8.0f / 5.0f) *
                    (current[index2(iz - 1, ix, domain.nz)] +
                     current[index2(iz + 1, ix, domain.nz)]) -
                (205.0f / 72.0f) * current[id];
            const float laplacian = d2x / dx2 + d2z / dz2;
            following[id] += domain.v2dt2[id] * laplacian;
        }
    }

    // Only the three rows adjacent to an artificial boundary use FD2 so
    // that a one-row Robin condition can enter the FD8 interior.
    auto advance_fd2_row = [&](int iz)
    {
        for (int ix = 4; ix < domain.nx - 4; ++ix)
        {
            const size_t id = index2(iz, ix, domain.nz);
            const float laplacian =
                (current[index2(iz, ix - 1, domain.nz)] -
                 2.0f * current[id] +
                 current[index2(iz, ix + 1, domain.nz)]) /
                    dx2 +
                (current[index2(iz - 1, ix, domain.nz)] -
                 2.0f * current[id] +
                 current[index2(iz + 1, ix, domain.nz)]) /
                    dz2;
            following[id] += domain.v2dt2[id] * laplacian;
        }
    };
    if (!domain.physical_top)
    {
        for (int iz = 1; iz < 4; ++iz)
            advance_fd2_row(iz);
    }
    if (!domain.physical_bottom)
    {
        for (int iz = domain.nz - 4; iz < domain.nz - 1; ++iz)
            advance_fd2_row(iz);
    }
}

void apply_absorber(
    const Domain &domain,
    std::vector<float> &following,
    const std::vector<float> &current,
    const std::vector<float> &previous)
{
    for (const AbsorbPoint &point : domain.absorb_points)
    {
        const float outgoing =
            2.0f *
                (point.t1 * current[point.id] +
                 point.t2 * current[point.id1] +
                 point.t3 * current[point.id2]) -
            (point.t1 * point.t1 * previous[point.id] +
             2.0f * point.t1 * point.t2 * previous[point.id1] +
             (2.0f * point.t1 * point.t3 +
              point.t2 * point.t2) *
                 previous[point.id2] +
             2.0f * point.t2 * point.t3 * previous[point.id3] +
             point.t3 * point.t3 * previous[point.id4]);
        following[point.id] =
            point.weight * outgoing +
            (1.0f - point.weight) * following[point.id];
    }
}

void apply_mur_top(
    const Domain &domain,
    std::vector<float> &following,
    const std::vector<float> &current)
{
    for (int ix = 1; ix < domain.nx - 1; ++ix)
    {
        const size_t top = index2(0, ix, domain.nz);
        const size_t inner = index2(1, ix, domain.nz);
        const float coefficient =
            (domain.velocity[top] * domain.dt - domain.dz) /
            (domain.velocity[top] * domain.dt + domain.dz);
        following[top] =
            current[inner] +
            coefficient * (following[inner] - current[top]);
    }
}

void apply_mur_bottom(
    const Domain &domain,
    std::vector<float> &following,
    const std::vector<float> &current)
{
    const int bottom_z = domain.nz - 1;
    for (int ix = 1; ix < domain.nx - 1; ++ix)
    {
        const size_t bottom = index2(bottom_z, ix, domain.nz);
        const size_t inner = index2(bottom_z - 1, ix, domain.nz);
        const float coefficient =
            (domain.velocity[bottom] * domain.dt - domain.dz) /
            (domain.velocity[bottom] * domain.dt + domain.dz);
        following[bottom] =
            current[inner] +
            coefficient * (following[inner] - current[bottom]);
    }
}

int first_effective_sample(
    const std::vector<float> *source,
    const std::vector<float> *top_robin,
    const std::vector<float> *bottom_robin,
    int nt,
    int nx,
    float gate)
{
    int first = nt;
    auto inspect = [&](const std::vector<float> *values, bool time_by_x)
    {
        if (values == nullptr || values->empty())
            return;
        float scale = 0.0f;
        for (float value : *values)
            scale = std::max(scale, std::fabs(value));
        if (!(scale > 0.0f))
            return;
        const float threshold =
            std::max(gate * scale, std::numeric_limits<float>::min());
        for (int it = 0; it < nt; ++it)
        {
            float amplitude = 0.0f;
            if (time_by_x)
            {
                for (int ix = 0; ix < nx; ++ix)
                    amplitude = std::max(
                        amplitude,
                        std::fabs((*values)[trace_index(it, ix, nx)]));
            }
            else
            {
                amplitude = std::fabs((*values)[it]);
            }
            if (amplitude >= threshold)
            {
                first = std::min(first, it);
                break;
            }
        }
    };
    inspect(source, false);
    inspect(top_robin, true);
    inspect(bottom_robin, true);
    return first;
}

State solve_local(
    const Domain &domain,
    const std::vector<int> &query_rows,
    int nt,
    const std::vector<float> *source,
    int source_global_z,
    int source_x,
    const std::vector<float> *top_robin,
    const std::vector<float> *bottom_robin,
    int receiver_global_z,
    bool record_receiver,
    float beta,
    float gate,
    int progress_interval,
    SolveStats &stats)
{
    State state;
    for (int row : query_rows)
        state.traces[row].assign(
            static_cast<size_t>(nt) * static_cast<size_t>(domain.nx),
            0.0f);
    if (record_receiver)
        state.receiver.assign(
            static_cast<size_t>(nt) *
                static_cast<size_t>(domain.physical_nx),
            0.0f);

    const int start = first_effective_sample(
        source,
        top_robin,
        bottom_robin,
        nt,
        domain.nx,
        gate);
    if (start >= nt)
        return state;

    std::vector<float> previous(domain.velocity.size(), 0.0f);
    std::vector<float> current(domain.velocity.size(), 0.0f);
    std::vector<float> following(domain.velocity.size(), 0.0f);
    const int source_local_z =
        domain.top_pad + source_global_z - domain.block.local_start;
    const int source_local_x = domain.nbc + source_x;
    const int receiver_local_z =
        domain.top_pad + receiver_global_z - domain.block.local_start;

    for (int it = start; it < nt; ++it)
    {
        if (source != nullptr &&
            source_local_z >= 0 && source_local_z < domain.nz &&
            source_local_x >= 0 && source_local_x < domain.nx)
        {
            current[index2(
                source_local_z,
                source_local_x,
                domain.nz)] += (*source)[it];
        }

        advance_wavefield(domain, previous, current, following);

        apply_absorber(domain, following, current, previous);

        if (!domain.physical_top)
        {
            if (top_robin == nullptr)
            {
                apply_mur_top(domain, following, current);
            }
            else
            {
                for (int ix = 1; ix < domain.nx - 1; ++ix)
                {
                    const size_t boundary = index2(0, ix, domain.nz);
                    const size_t inner = index2(1, ix, domain.nz);
                    following[boundary] =
                        current[boundary] +
                        domain.dt *
                            ((*top_robin)[trace_index(it, ix, domain.nx)] +
                             beta * domain.velocity[boundary] *
                                 (current[inner] - current[boundary]) /
                                 domain.dz);
                }
                following[index2(0, 0, domain.nz)] =
                    following[index2(0, 1, domain.nz)];
                following[index2(0, domain.nx - 1, domain.nz)] =
                    following[index2(0, domain.nx - 2, domain.nz)];
            }
        }

        if (!domain.physical_bottom)
        {
            const int bottom_z = domain.nz - 1;
            if (bottom_robin == nullptr)
            {
                apply_mur_bottom(domain, following, current);
            }
            else
            {
                for (int ix = 1; ix < domain.nx - 1; ++ix)
                {
                    const size_t boundary =
                        index2(bottom_z, ix, domain.nz);
                    const size_t inner =
                        index2(bottom_z - 1, ix, domain.nz);
                    following[boundary] =
                        current[boundary] +
                        domain.dt *
                            ((*bottom_robin)[
                                 trace_index(it, ix, domain.nx)] -
                             beta * domain.velocity[boundary] *
                                 (current[boundary] - current[inner]) /
                                 domain.dz);
                }
                following[index2(bottom_z, 0, domain.nz)] =
                    following[index2(bottom_z, 1, domain.nz)];
                following[index2(
                    bottom_z,
                    domain.nx - 1,
                    domain.nz)] =
                    following[index2(
                        bottom_z,
                        domain.nx - 2,
                        domain.nz)];
            }
        }

        previous.swap(current);
        current.swap(following);

        for (auto &entry : state.traces)
        {
            const int local_z =
                domain.top_pad +
                entry.first -
                domain.block.local_start;
            for (int ix = 0; ix < domain.nx; ++ix)
            {
                entry.second[trace_index(it, ix, domain.nx)] =
                    current[index2(local_z, ix, domain.nz)];
            }
        }
        if (record_receiver)
        {
            for (int ix = 0; ix < domain.physical_nx; ++ix)
            {
                state.receiver[trace_index(
                    it,
                    ix,
                    domain.physical_nx)] =
                    current[index2(
                        receiver_local_z,
                        ix + domain.nbc,
                        domain.nz)];
            }
        }
        if (progress_interval > 0 &&
            (it == start ||
             (it + 1) % progress_interval == 0 ||
             it == nt - 1))
        {
            INFO((
                "local block %d: step %d/%d",
                domain.block.index + 1,
                it + 1,
                nt));
        }
    }

    stats.active_steps = nt - start;
    stats.updated_points =
        static_cast<long long>(stats.active_steps) *
        static_cast<long long>(domain.nx) *
        static_cast<long long>(domain.nz);
    return state;
}

void add_state(State &target, const State &increment)
{
    for (const auto &entry : increment.traces)
    {
        std::vector<float> &values = target.traces.at(entry.first);
        for (size_t index = 0; index < values.size(); ++index)
            values[index] += entry.second[index];
    }
    if (!increment.receiver.empty())
    {
        if (target.receiver.empty())
            target.receiver.assign(increment.receiver.size(), 0.0f);
        for (size_t index = 0; index < target.receiver.size(); ++index)
            target.receiver[index] += increment.receiver[index];
    }
}

std::vector<float> extended_velocity_row(
    const sewave::Grid2D &grid,
    int row,
    int nbc)
{
    std::vector<float> velocity(
        static_cast<size_t>(grid.nx + 2 * nbc));
    for (int ix = 0; ix < static_cast<int>(velocity.size()); ++ix)
    {
        const int global_x = clamp_int(ix - nbc, 0, grid.nx - 1);
        velocity[ix] = grid.v[index2(row, global_x, grid.nz)];
    }
    return velocity;
}

std::vector<float> robin_operator(
    const State &state,
    int row,
    const std::vector<float> &velocity,
    int nt,
    int nx,
    float dt,
    float dz,
    float sign,
    float beta)
{
    const std::vector<float> &center = state.traces.at(row);
    const std::vector<float> &neighbor =
        state.traces.at(sign < 0.0f ? row + 1 : row - 1);
    std::vector<float> output(
        static_cast<size_t>(nt) * static_cast<size_t>(nx),
        0.0f);

    for (int it = 1; it < nt; ++it)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            const size_t now = trace_index(it, ix, nx);
            const size_t before = trace_index(it - 1, ix, nx);
            float derivative = 0.0f;
            if (sign < 0.0f)
                derivative = (neighbor[before] - center[before]) / dz;
            else
                derivative = (center[before] - neighbor[before]) / dz;
            output[now] =
                (center[now] - center[before]) / dt +
                sign * beta * velocity[ix] * derivative;
        }
    }
    return output;
}

std::vector<float> subtract(
    const std::vector<float> &left,
    const std::vector<float> &right)
{
    if (left.size() != right.size())
        throw std::runtime_error("interface vectors have different sizes");
    std::vector<float> output(left.size());
    for (size_t index = 0; index < output.size(); ++index)
        output[index] = left[index] - right[index];
    return output;
}

std::vector<float> fmm_travel_time(
    const sewave::Grid2D &grid,
    float source_x,
    float source_z)
{
    efmm_t context;
    std::memset(&context, 0, sizeof(context));

    // efmm's first array index is contiguous. RSF/sewave stores depth as
    // the contiguous n1 axis, so depth is passed as efmm's x-like axis.
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
    std::vector<float> travel_time(
        context.tt,
        context.tt +
            static_cast<size_t>(grid.nz) *
                static_cast<size_t>(grid.nx));
    efmm_free(&context);
    return travel_time;
}

State ray_state(
    const sewave::Grid2D &grid,
    const Domain &domain,
    const std::vector<int> &query_rows,
    const std::vector<float> &travel_time,
    int nt,
    float dt,
    float frequency,
    float source_velocity,
    float scale)
{
    State state;
    const float spacing = std::min(grid.dx, grid.dz);
    for (int row : query_rows)
    {
        std::vector<float> &trace = state.traces[row];
        trace.assign(
            static_cast<size_t>(nt) * static_cast<size_t>(domain.nx),
            0.0f);
        for (int ix = 0; ix < domain.nx; ++ix)
        {
            const int global_x =
                clamp_int(ix - domain.nbc, 0, grid.nx - 1);
            const int outside_cells =
                std::abs(ix - domain.nbc - global_x);
            const float local_velocity =
                grid.v[index2(row, global_x, grid.nz)];
            const float time =
                travel_time[index2(row, global_x, grid.nz)] +
                outside_cells * grid.dx /
                    std::max(local_velocity, 1.0f);
            const float distance =
                std::max(source_velocity * time, spacing);
            const float amplitude =
                scale * std::sqrt(spacing / distance);
            for (int it = 0; it < nt; ++it)
            {
                trace[trace_index(it, ix, domain.nx)] =
                    amplitude *
                    ricker(it * dt - time, frequency);
            }
        }
    }
    return state;
}

std::vector<float> ray_receiver(
    const sewave::Grid2D &grid,
    const std::vector<float> &travel_time,
    int receiver_row,
    int nt,
    float dt,
    float frequency,
    float source_velocity,
    float scale)
{
    std::vector<float> output(
        static_cast<size_t>(nt) * static_cast<size_t>(grid.nx),
        0.0f);
    const float spacing = std::min(grid.dx, grid.dz);
    for (int ix = 0; ix < grid.nx; ++ix)
    {
        const float time =
            travel_time[index2(receiver_row, ix, grid.nz)];
        const float distance =
            std::max(source_velocity * time, spacing);
        const float amplitude =
            scale * std::sqrt(spacing / distance);
        for (int it = 0; it < nt; ++it)
        {
            output[trace_index(it, ix, grid.nx)] =
                amplitude * ricker(it * dt - time, frequency);
        }
    }
    return output;
}

std::vector<float> run_shot(
    const sewave::Grid2D &grid,
    const std::vector<Block> &blocks,
    const std::vector<Domain> &domains,
    const std::vector<std::vector<int>> &queries,
    const Parameters &parameters,
    float source_x,
    int source_z,
    int source_x_index,
    int receiver_z,
    std::vector<float> &initial_ray_record,
    std::vector<float> &travel_time,
    ShotStats &shot_stats)
{
    const auto fmm_begin = std::chrono::steady_clock::now();
    travel_time = fmm_travel_time(
        grid,
        source_x,
        grid.z0 + source_z * grid.dz);
    shot_stats.fmm_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - fmm_begin)
            .count();

    std::vector<float> source(static_cast<size_t>(parameters.nt));
    for (int it = 0; it < parameters.nt; ++it)
        source[it] = ricker(it * parameters.dt, parameters.fdom);

    const float source_velocity =
        grid.v[index2(source_z, source_x_index, grid.nz)];
    const auto initial_begin = std::chrono::steady_clock::now();
    std::vector<State> states(static_cast<size_t>(parameters.blocks));
    SolveStats source_solve;
    states[0] = solve_local(
        domains[0],
        queries[0],
        parameters.nt,
        &source,
        source_z,
        source_x_index,
        nullptr,
        nullptr,
        receiver_z,
        true,
        parameters.beta,
        parameters.gate,
        parameters.progress_interval,
        source_solve);
    shot_stats.updated_points += source_solve.updated_points;
    ++shot_stats.local_solves;

    shot_stats.ray_scale = parameters.ray_scale;

    for (int ib = 1; ib < parameters.blocks; ++ib)
    {
        states[ib] = ray_state(
            grid,
            domains[ib],
            queries[ib],
            travel_time,
            parameters.nt,
            parameters.dt,
            parameters.fdom,
            source_velocity,
            static_cast<float>(shot_stats.ray_scale));
    }
    initial_ray_record = ray_receiver(
        grid,
        travel_time,
        receiver_z,
        parameters.nt,
        parameters.dt,
        parameters.fdom,
        source_velocity,
        static_cast<float>(shot_stats.ray_scale));
    shot_stats.initial_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - initial_begin)
            .count();

    const auto sweep_begin = std::chrono::steady_clock::now();
    double final_residual_numerator = 0.0;
    double final_residual_denominator = 0.0;
    for (int cycle = 0; cycle < parameters.cycles; ++cycle)
    {
        final_residual_numerator = 0.0;
        final_residual_denominator = 0.0;

        // Multiplicative downward sweep: each block reads only the already
        // updated block above it. No lower-block reflection can enter here.
        for (int ib = 1; ib < parameters.blocks; ++ib)
        {
            const int row = blocks[ib].local_start;
            const std::vector<float> velocity =
                extended_velocity_row(grid, row, parameters.nbc);
            const std::vector<float> neighbor = robin_operator(
                states[ib - 1],
                row,
                velocity,
                parameters.nt,
                domains[ib].nx,
                parameters.dt,
                grid.dz,
                -1.0f,
                parameters.beta);
            const std::vector<float> own = robin_operator(
                states[ib],
                row,
                velocity,
                parameters.nt,
                domains[ib].nx,
                parameters.dt,
                grid.dz,
                -1.0f,
                parameters.beta);
            const std::vector<float> forcing = subtract(neighbor, own);
            for (size_t index = 0; index < forcing.size(); ++index)
            {
                final_residual_numerator +=
                    static_cast<double>(forcing[index]) *
                    static_cast<double>(forcing[index]);
                final_residual_denominator +=
                    static_cast<double>(neighbor[index]) *
                    static_cast<double>(neighbor[index]);
            }

            SolveStats solve_stats;
            State increment = solve_local(
                domains[ib],
                queries[ib],
                parameters.nt,
                nullptr,
                -100000,
                source_x_index,
                &forcing,
                nullptr,
                receiver_z,
                false,
                parameters.beta,
                parameters.gate,
                parameters.progress_interval,
                solve_stats);
            add_state(states[ib], increment);
            shot_stats.updated_points += solve_stats.updated_points;
            ++shot_stats.local_solves;
        }

        // Multiplicative upward sweep returns reflections/scattering created
        // in deeper blocks to the surface block.
        for (int ib = parameters.blocks - 2; ib >= 0; --ib)
        {
            const int row = blocks[ib].local_end - 1;
            const std::vector<float> velocity =
                extended_velocity_row(grid, row, parameters.nbc);
            const std::vector<float> neighbor = robin_operator(
                states[ib + 1],
                row,
                velocity,
                parameters.nt,
                domains[ib].nx,
                parameters.dt,
                grid.dz,
                +1.0f,
                parameters.beta);
            const std::vector<float> own = robin_operator(
                states[ib],
                row,
                velocity,
                parameters.nt,
                domains[ib].nx,
                parameters.dt,
                grid.dz,
                +1.0f,
                parameters.beta);
            const std::vector<float> forcing = subtract(neighbor, own);
            for (size_t index = 0; index < forcing.size(); ++index)
            {
                final_residual_numerator +=
                    static_cast<double>(forcing[index]) *
                    static_cast<double>(forcing[index]);
                final_residual_denominator +=
                    static_cast<double>(neighbor[index]) *
                    static_cast<double>(neighbor[index]);
            }

            SolveStats solve_stats;
            State increment = solve_local(
                domains[ib],
                queries[ib],
                parameters.nt,
                nullptr,
                -100000,
                source_x_index,
                nullptr,
                &forcing,
                receiver_z,
                ib == 0,
                parameters.beta,
                parameters.gate,
                parameters.progress_interval,
                solve_stats);
            add_state(states[ib], increment);
            shot_stats.updated_points += solve_stats.updated_points;
            ++shot_stats.local_solves;
        }
        shot_stats.final_interface_residual =
            std::sqrt(
                final_residual_numerator /
                std::max(
                    final_residual_denominator,
                    std::numeric_limits<double>::min()));
        INFO((
            "shot x=%g cycle %d/%d: Robin residual=%e",
            source_x,
            cycle + 1,
            parameters.cycles,
            shot_stats.final_interface_residual));
    }
    shot_stats.sweep_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - sweep_begin)
            .count();
    return states[0].receiver;
}

int coordinate_index(float coordinate, float origin, float spacing)
{
    return static_cast<int>(
        std::lround((coordinate - origin) / spacing));
}

void validate(
    const sewave::Grid2D &grid,
    const Parameters &parameters,
    int source_z,
    int receiver_z,
    const std::vector<Block> &blocks)
{
    if (grid.nz < 16 || grid.nx < 16)
        throw std::runtime_error("velocity grid is too small for FD8");
    if (parameters.nt <= 0 || parameters.dt <= 0.0f ||
        parameters.fdom <= 0.0f)
        throw std::runtime_error("nt, dt, and fdom must be positive");
    if (parameters.blocks < 2)
        throw std::runtime_error("blocks must be at least 2");
    if (parameters.blocks > grid.nz / 4)
        throw std::runtime_error("too many blocks for the depth grid");
    if (parameters.overlap < 4)
        throw std::runtime_error("overlap must be at least 4 for FD8");
    if (parameters.cycles < 1)
        throw std::runtime_error("cycles must be at least 1");
    if (parameters.beta <= 0.0f)
        throw std::runtime_error("beta must be positive");
    if (parameters.ray_scale <= 0.0f)
        throw std::runtime_error("ray_scale must be positive");
    if (parameters.gate < 0.0f || parameters.gate >= 1.0f)
        throw std::runtime_error("gate must satisfy 0 <= gate < 1");
    if (parameters.nbc < parameters.absorb_length ||
        parameters.absorb_length < 5)
        throw std::runtime_error("require nbc >= L >= 5");
    if (source_z < 0 || source_z >= grid.nz ||
        receiver_z < 0 || receiver_z >= grid.nz)
        throw std::runtime_error("source or receiver depth is outside the model");
    if (source_z >= blocks.front().core_end)
        throw std::runtime_error(
            "source must lie in the first core block");
    if (receiver_z >= blocks.front().core_end)
        throw std::runtime_error(
            "receiver row must lie in the first core block");

    const float maximum_velocity =
        *std::max_element(grid.v.begin(), grid.v.end());
    const double cfl =
        maximum_velocity * parameters.dt *
        std::sqrt(
            1.0 / (grid.dx * grid.dx) +
            1.0 / (grid.dz * grid.dz));
    if (cfl > 0.45)
    {
        INFO((
            "WARNING: FD8 CFL=%g exceeds the conservative 0.45 limit",
            cfl));
    }
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
            throw std::runtime_error("Need velocity= RSF velocity model");
        if (!se_have_par("output"))
            throw std::runtime_error("Need output= output seismic RSF");

        const char *velocity_path = se_get_par_str("velocity");
        const char *output_path = se_get_par_str("output");
        sewave::Grid2D grid = sewave::read_rsf2d(velocity_path);

        Parameters parameters;
        parameters.nt =
            se_have_par("nt") ? se_get_par_int("nt") : 1000;
        parameters.dt =
            se_have_par("dt") ? se_get_par_float("dt") : 0.001f;
        parameters.fdom =
            se_have_par("fdom") ? se_get_par_float("fdom") : 20.0f;
        parameters.sx =
            se_have_par("sx")
                ? se_get_par_float("sx")
                : grid.x0 + 0.5f * (grid.nx - 1) * grid.dx;
        parameters.sz =
            se_have_par("sz") ? se_get_par_float("sz") : grid.z0;
        parameters.rz =
            se_have_par("rz") ? se_get_par_float("rz") : parameters.sz;
        parameters.ns =
            se_have_par("ns") ? se_get_par_int("ns") : 1;
        parameters.ds =
            se_have_par("ds") ? se_get_par_float("ds") : 0.0f;
        parameters.blocks =
            se_have_par("blocks") ? se_get_par_int("blocks") : 4;
        parameters.overlap =
            se_have_par("overlap") ? se_get_par_int("overlap") : 22;
        parameters.cycles =
            se_have_par("cycles") ? se_get_par_int("cycles") : 1;
        parameters.beta =
            se_have_par("beta") ? se_get_par_float("beta") : 1.0f;
        parameters.gate =
            se_have_par("gate") ? se_get_par_float("gate") : 0.003f;
        parameters.nbc =
            se_have_par("nbc") ? se_get_par_int("nbc") : 40;
        parameters.absorb_length =
            se_have_par("L") ? se_get_par_int("L") : 30;
        parameters.absorb_alpha =
            se_have_par("alpha") ? se_get_par_float("alpha") : 1.0f;
        parameters.ray_scale =
            se_have_par("ray_scale")
                ? se_get_par_float("ray_scale")
                : 1.0f;
        parameters.progress_interval =
            se_have_par("progress_interval")
                ? se_get_par_int("progress_interval")
                : 0;

        const int source_z =
            coordinate_index(parameters.sz, grid.z0, grid.dz);
        const int receiver_z =
            coordinate_index(parameters.rz, grid.z0, grid.dz);
        const std::vector<Block> blocks = make_blocks(
            grid.nz,
            parameters.blocks,
            parameters.overlap);
        validate(grid, parameters, source_z, receiver_z, blocks);

        std::vector<Domain> domains;
        std::vector<std::vector<int>> queries;
        domains.reserve(parameters.blocks);
        queries.reserve(parameters.blocks);
        for (int ib = 0; ib < parameters.blocks; ++ib)
        {
            domains.push_back(
                make_domain(grid, blocks[ib], parameters));
            queries.push_back(robin_query_rows(blocks, ib));
            INFO((
                "block %d/%d: core=[%d,%d) local=[%d,%d) overlap=%d",
                ib + 1,
                parameters.blocks,
                blocks[ib].core_start,
                blocks[ib].core_end,
                blocks[ib].local_start,
                blocks[ib].local_end,
                parameters.overlap));
        }

        sewave::Data3D output;
        output.nt = parameters.nt;
        output.nr = grid.nx;
        output.ns = parameters.ns;
        output.dt = parameters.dt;
        output.dr = grid.dx;
        output.ds = parameters.ds;
        output.t0 = 0.0f;
        output.r0 = grid.x0;
        output.s0 = parameters.sx;
        output.d.assign(
            static_cast<size_t>(output.nt) *
                static_cast<size_t>(output.nr) *
                static_cast<size_t>(output.ns),
            0.0f);
        sewave::Data3D ray_output = output;

        std::vector<float> first_travel_time;
        const auto total_begin = std::chrono::steady_clock::now();
        double fmm_total = 0.0;
        double initial_total = 0.0;
        double sweep_total = 0.0;
        long long update_total = 0;
        int solve_total = 0;

        for (int is = 0; is < parameters.ns; ++is)
        {
            const float shot_x = parameters.sx + is * parameters.ds;
            const int shot_x_index =
                coordinate_index(shot_x, grid.x0, grid.dx);
            if (shot_x_index < 0 || shot_x_index >= grid.nx)
                throw std::runtime_error("shot is outside the velocity model");

            std::vector<float> ray_record;
            std::vector<float> travel_time;
            ShotStats shot_stats;
            const std::vector<float> record = run_shot(
                grid,
                blocks,
                domains,
                queries,
                parameters,
                shot_x,
                source_z,
                shot_x_index,
                receiver_z,
                ray_record,
                travel_time,
                shot_stats);
            if (is == 0)
                first_travel_time = travel_time;

            for (int ix = 0; ix < grid.nx; ++ix)
            {
                for (int it = 0; it < parameters.nt; ++it)
                {
                    output.d[
                        (static_cast<size_t>(is) * grid.nx + ix) *
                            parameters.nt +
                        it] =
                        record[trace_index(it, ix, grid.nx)];
                    ray_output.d[
                        (static_cast<size_t>(is) * grid.nx + ix) *
                            parameters.nt +
                        it] =
                        ray_record[trace_index(it, ix, grid.nx)];
                }
            }

            fmm_total += shot_stats.fmm_seconds;
            initial_total += shot_stats.initial_seconds;
            sweep_total += shot_stats.sweep_seconds;
            update_total += shot_stats.updated_points;
            solve_total += shot_stats.local_solves;
            INFO((
                "shot %d/%d: ray_scale=%g fmm=%.6f s initial=%.6f s "
                "sweeps=%.6f s solves=%d updated_points=%lld "
                "final_robin_residual=%e",
                is + 1,
                parameters.ns,
                shot_stats.ray_scale,
                shot_stats.fmm_seconds,
                shot_stats.initial_seconds,
                shot_stats.sweep_seconds,
                shot_stats.local_solves,
                shot_stats.updated_points,
                shot_stats.final_interface_residual));
        }

        sewave::write_rsf3d(output_path, output);
        if (se_have_par("ray_output"))
            sewave::write_rsf3d(
                se_get_par_str("ray_output"),
                ray_output);
        if (se_have_par("traveltime"))
        {
            if (parameters.ns != 1)
                throw std::runtime_error(
                    "traveltime= currently requires ns=1");
            sewave::write_rsf2d(
                se_get_par_str("traveltime"),
                grid,
                first_travel_time,
                "FMM first-arrival time");
        }

        const double total_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - total_begin)
                .count();
        INFO((
            "block forward summary: total=%.6f s fmm=%.6f s "
            "initial=%.6f s sweeps=%.6f s local_solves=%d "
            "updated_points=%lld",
            total_seconds,
            fmm_total,
            initial_total,
            sweep_total,
            solve_total,
            update_total));
        se_par_destroy();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::fprintf(stderr, "block_iterative_forward: %s\n", error.what());
        se_par_destroy();
        return 1;
    }
}