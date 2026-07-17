#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_par_sep.h>
#include "../SEWAVE/sewave2d.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <sstream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
void truncate_nt(sewave::Data3D& d, int max_nt) {
    if (max_nt <= 0 || max_nt >= d.nt) return;
    std::vector<float> out((size_t)d.ns * d.nr * max_nt, 0.0f);
    for (int is = 0; is < d.ns; ++is)
        for (int ir = 0; ir < d.nr; ++ir)
            std::copy_n(&d.d[((size_t)is * d.nr + ir) * d.nt], max_nt,
                        &out[((size_t)is * d.nr + ir) * max_nt]);
    d.nt = max_nt;
    d.d.swap(out);
}



std::vector<int> parse_int_list(const char* text) {
    std::vector<int> values;
    if (!text) return values;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) values.push_back(std::max(1, std::atoi(item.c_str())));
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

void check_stability(const sewave::Grid2D& vel, float dt, int type_laplace, const char* label) {
    if (dt <= 0.0f || vel.v.empty() || vel.dx <= 0.0f || vel.dz <= 0.0f) return;
    const float vmax = *std::max_element(vel.v.begin(), vel.v.end());
    const double inv_h = std::sqrt(1.0 / (vel.dx * vel.dx) + 1.0 / (vel.dz * vel.dz));
    const double cfl = vmax * dt * inv_h;
    double dt_limit = 0.45 / (vmax * inv_h);
    const char* method = "FD8";
    if (type_laplace == 1) {
        const double kmax = std::sqrt((M_PI / vel.dx) * (M_PI / vel.dx) + (M_PI / vel.dz) * (M_PI / vel.dz));
        dt_limit = 2.0 / (vmax * kmax);
        method = "PS";
    }
    if (dt > dt_limit) {
        INFO(("WARNING: %s stability check may be unstable for %s: vmax=%.6g dx=%.6g dz=%.6g dt=%.9g CFL=%.6g; suggested dt <= %.9g",
              label, method, vmax, vel.dx, vel.dz, dt, cfl, dt_limit));
    } else {
        INFO(("%s stability check OK for %s: vmax=%.6g dx=%.6g dz=%.6g dt=%.9g CFL=%.6g; limit %.9g",
              label, method, vmax, vel.dx, vel.dz, dt, cfl, dt_limit));
    }

}

void resample_time(sewave::Data3D& d, float resample_dt) {
    if (resample_dt <= 0.0f || d.nt <= 1 || std::fabs(resample_dt - d.dt) <= 1e-12f) return;
    const int old_nt = d.nt;
    const float old_dt = d.dt;
    const int new_nt = std::max(1, (int)std::floor(((old_nt - 1) * old_dt) / resample_dt + 1.0e-6f) + 1);
    std::vector<float> out((size_t)d.ns * d.nr * new_nt, 0.0f);
    for (int is = 0; is < d.ns; ++is) {
        for (int ir = 0; ir < d.nr; ++ir) {
            const float* in = &d.d[((size_t)is * d.nr + ir) * old_nt];
            float* trace = &out[((size_t)is * d.nr + ir) * new_nt];
            for (int it = 0; it < new_nt; ++it) {
                const float pos = it * resample_dt / old_dt;
                const int i0 = std::min(old_nt - 1, (int)std::floor(pos));
                const int i1 = std::min(old_nt - 1, i0 + 1);
                const float a = pos - i0;
                trace[it] = (1.0f - a) * in[i0] + a * in[i1];
            }
        }
    }
    d.nt = new_nt;
    d.dt = resample_dt;
    d.d.swap(out);
}
}

int main(int argc, char** argv) {
    se_par_init(argc, argv);
    try {
#ifdef _OPENMP
        INFO(("OpenMP enabled, max_threads=%d", omp_get_max_threads()));
#else
        INFO(("OpenMP disabled at compile time"));
#endif
        const char* mode = se_have_par("mode") ? se_get_par_str("mode") : "forward";
        if (!se_have_par("velocity")) ERROR(("Need velocity= RSF velocity model"));
        sewave::Grid2D vel = sewave::read_rsf2d(se_get_par_str("velocity"));
        auto load_q = [&]() -> sewave::Grid2D {
            if (se_have_par("qfile")) return sewave::read_rsf2d(se_get_par_str("qfile"));
            if (se_have_par("Qfile")) return sewave::read_rsf2d(se_get_par_str("Qfile"));
            ERROR(("visco=1 requires qfile= (or Qfile=) RSF Q model matching velocity="));
            return sewave::Grid2D{};
        };
        if (mode[0]=='f') {
            if (!se_have_par("output")) ERROR(("Need output= output seismic RSF"));
            sewave::ForwardParams p;
            p.nt = se_have_par("nt") ? se_get_par_int("nt") : 1000;
            p.dt = se_have_par("dt") ? se_get_par_float("dt") : 0.001f;
            p.fdom = se_have_par("fdom") ? se_get_par_float("fdom") : 20.0f;
            p.sx = se_have_par("sx") ? se_get_par_float("sx") : vel.x0 + 0.5f*(vel.nx-1)*vel.dx;
            p.sz = se_have_par("sz") ? se_get_par_float("sz") : vel.z0;
            p.rz = se_have_par("rz") ? se_get_par_float("rz") : p.sz;
            p.nr = se_have_par("nr") ? se_get_par_int("nr") : vel.nx;
            p.r0 = se_have_par("r0") ? se_get_par_float("r0") : vel.x0;
            p.dr = se_have_par("dr") ? se_get_par_float("dr") : vel.dx;
            p.ns = se_have_par("ns") ? se_get_par_int("ns") : 1;
            p.ds = se_have_par("ds") ? se_get_par_float("ds") : 0.0f;
            p.nbc = se_have_par("nbc") ? se_get_par_int("nbc") : 100;
            p.alpha = se_have_par("alpha") ? se_get_par_float("alpha") : 1.0f;
            p.type_compute_laplace = se_have_par("type_compute_Laplace") ? se_get_par_int("type_compute_Laplace") : (se_have_par("type_compute_laplace") ? se_get_par_int("type_compute_laplace") : 1);
            p.L = se_have_par("L") ? se_get_par_int("L") : 30;
            p.type = se_have_par("type") ? se_get_par_int("type") : 1;
            p.f0 = se_have_par("f0") ? se_get_par_float("f0") : 20.0f;
            p.Omega0 = se_have_par("Omega0") ? se_get_par_float("Omega0") : 0.0f;
            p.fc = se_have_par("fc") ? se_get_par_float("fc") : 120.0f;
            p.order = se_have_par("order") ? se_get_par_int("order") : 2;
            p.flag_smooth = se_have_par("flag_smooth") ? se_get_par_int("flag_smooth") != 0 : true;
            p.flag_homo = se_have_par("flag_homo") ? se_get_par_int("flag_homo") != 0 : false;
            p.visco = se_have_par("visco") ? se_get_par_int("visco") != 0 : false;
            p.q = se_have_par("q") ? se_get_par_float("q") : 1000.0f;
            p.progress_interval = se_have_par("progress_interval") ? se_get_par_int("progress_interval") : (se_have_par("log_interval") ? se_get_par_int("log_interval") : std::max(1, p.nt/10));
            sewave::Grid2D qmodel;
            sewave::Grid2D* qptr = nullptr;
            if (p.visco || p.type == 1) { qmodel = load_q(); qptr = &qmodel; }
            check_stability(vel, p.dt, p.type_compute_laplace, "forward");
            sewave::Data3D d = sewave::forward(vel, p, qptr);
            sewave::write_rsf3d(se_get_par_str("output"), d);
        } else if (mode[0]=='i') {
            if (!se_have_par("seismic_data")) ERROR(("Need seismic_data= input seismic RSF"));
            if (!se_have_par("migration")) ERROR(("Need migration= output image RSF"));
            sewave::Data3D d = sewave::read_rsf3d(se_get_par_str("seismic_data"));
            int max_nt = se_have_par("max_nt") ? se_get_par_int("max_nt") : -1;
            float resample_dt = se_have_par("resample_dt") ? se_get_par_float("resample_dt") : -1.0f;
            const int input_nt = d.nt;
            const float input_dt = d.dt;
            truncate_nt(d, max_nt);
            resample_time(d, resample_dt);
            if (d.nt != input_nt || d.dt != input_dt) {
                INFO(("seismic_data time preprocessing: nt %d -> %d, dt %.9g -> %.9g", input_nt, d.nt, input_dt, d.dt));
            }
            sewave::ImageParams p;
            p.cmp = se_have_par("cmp") ? se_get_par_int("cmp") != 0 : true;
            p.compensate = se_have_par("compensate") ? se_get_par_int("compensate") != 0 : false;
            p.amp_compensation_sign = se_have_par("amp_compensation_sign") ? se_get_par_int("amp_compensation_sign") : (p.compensate ? -1 : 1);
            if (p.amp_compensation_sign != -1 && p.amp_compensation_sign != 1) ERROR(("amp_compensation_sign must be -1 (Q compensation) or 1 (Q attenuation)"));
            p.visco = se_have_par("visco") ? se_get_par_int("visco") != 0 : false;
            p.q = se_have_par("q") ? se_get_par_float("q") : 1000.0f;
            p.nbc = se_have_par("nbc") ? se_get_par_int("nbc") : 100;
            p.alpha = se_have_par("alpha") ? se_get_par_float("alpha") : 1.0f;
            p.type_compute_laplace = se_have_par("type_compute_Laplace") ? se_get_par_int("type_compute_Laplace") : (se_have_par("type_compute_laplace") ? se_get_par_int("type_compute_laplace") : 1);
            p.L = se_have_par("L") ? se_get_par_int("L") : 30;
            p.type = se_have_par("type") ? se_get_par_int("type") : 1;
            p.f0 = se_have_par("f0") ? se_get_par_float("f0") : 20.0f;
            p.Omega0 = se_have_par("Omega0") ? se_get_par_float("Omega0") : 0.0f;
            p.fc = se_have_par("fc") ? se_get_par_float("fc") : 120.0f;
            p.order = se_have_par("order") ? se_get_par_int("order") : 2;
            p.flag_smooth = se_have_par("flag_smooth") ? se_get_par_int("flag_smooth") != 0 : true;
            p.flag_homo = se_have_par("flag_homo") ? se_get_par_int("flag_homo") != 0 : false;
            p.sz = se_have_par("sz") ? se_get_par_float("sz") : vel.z0;
            p.rz = se_have_par("rz") ? se_get_par_float("rz") : p.sz;
            p.fdom = se_have_par("fdom") ? se_get_par_float("fdom") : 20.0f;
            p.shot_begin = se_have_par("shot_begin") ? se_get_par_int("shot_begin") : (se_have_par("first_shot") ? se_get_par_int("first_shot") : 0);
            p.shot_end = se_have_par("shot_end") ? se_get_par_int("shot_end") : (se_have_par("last_shot") ? se_get_par_int("last_shot") : -1);
            p.progress_interval = se_have_par("progress_interval") ? se_get_par_int("progress_interval") : (se_have_par("log_interval") ? se_get_par_int("log_interval") : std::max(1, d.nt/10));
            p.max_parallel_shots = se_have_par("max_parallel_shots") ? se_get_par_int("max_parallel_shots") : (se_have_par("parallel_shots") ? se_get_par_int("parallel_shots") : (se_have_par("omp_shots") ? se_get_par_int("omp_shots") : 1));
            p.in_memory_snapshots = se_have_par("in_memory_snapshots") ? se_get_par_int("in_memory_snapshots") != 0 : false;
            if (se_have_par("debug_wavefield_prefix")) p.debug_wavefield_prefix = se_get_par_str("debug_wavefield_prefix");
            else if (se_have_par("wavefield_snapshot_prefix")) p.debug_wavefield_prefix = se_get_par_str("wavefield_snapshot_prefix");
            if (se_have_par("debug_snapshot_steps")) {
                int nsteps = 0;
                int* steps = se_get_pararray_int("debug_snapshot_steps", &nsteps);
                if (steps && nsteps > 0) p.debug_snapshot_steps.assign(steps, steps + nsteps);
                else p.debug_snapshot_steps = parse_int_list(se_get_par_str("debug_snapshot_steps"));
            } else if (se_have_par("wavefield_snapshot_steps")) {
                int nsteps = 0;
                int* steps = se_get_pararray_int("wavefield_snapshot_steps", &nsteps);
                if (steps && nsteps > 0) p.debug_snapshot_steps.assign(steps, steps + nsteps);
                else p.debug_snapshot_steps = parse_int_list(se_get_par_str("wavefield_snapshot_steps"));
            }
            for (int& step : p.debug_snapshot_steps) step = std::max(1, step);
            std::sort(p.debug_snapshot_steps.begin(), p.debug_snapshot_steps.end());
            p.debug_snapshot_steps.erase(std::unique(p.debug_snapshot_steps.begin(), p.debug_snapshot_steps.end()), p.debug_snapshot_steps.end());
            if (!p.debug_wavefield_prefix.empty() && p.debug_snapshot_steps.empty()) INFO(("debug_wavefield_prefix is set but debug_snapshot_steps is empty; no wavefield snapshots will be written"));
            sewave::Grid2D qmodel;
            sewave::Grid2D* qptr = nullptr;
            if (p.visco || p.type == 1) { qmodel = load_q(); qptr = &qmodel; }
            check_stability(vel, d.dt, p.type_compute_laplace, "image");
            auto img = sewave::rtm_image(vel, d, p, qptr);
            sewave::write_rsf2d(se_get_par_str("migration"), vel, img, "RTM image");
        } else {
            ERROR(("Unknown mode=%s (use forward or image)", mode));
        }
    } catch (const std::exception& e) {
        ERROR(("sewave2d failed: %s", e.what()));
    }
    se_par_destroy();
    return 0;
}
