#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_par_sep.h>
#include "../SEWAVE/sewave2d.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <sstream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
void print_usage(const char* prog) {
    const char* name = (prog && prog[0]) ? prog : "sewave2d";
    std::printf(
        "Usage:\n"
        "  %s mode=forward velocity=vel.rsf output=data.rsf [forward parameters]\n"
        "  %s mode=image velocity=vel.rsf seismic_data=data.rsf migration=image.rsf [imaging parameters]\n"
        "\n"
        "Required parameters:\n"
        "  mode=forward|image          Run synthetic modeling or RTM imaging. Default: forward.\n"
        "  velocity=FILE              2D RSF velocity model.\n"
        "  output=FILE                Forward-mode output seismic RSF.\n"
        "  seismic_data=FILE          Image-mode input seismic RSF.\n"
        "  migration=FILE             Image-mode output RTM image RSF.\n"
        "\n"
        "Common wave-equation parameters:\n"
        "  dt=FLOAT                   Time step for forward modeling. Default: 0.001.\n"
        "  nt=INT                     Number of time samples for forward modeling. Default: 1000.\n"
        "  fdom=FLOAT                 Ricker dominant frequency. Default: 20.\n"
        "  sz=FLOAT                   Source depth. Default: velocity z origin.\n"
        "  rz=FLOAT                   Receiver depth. Default: sz.\n"
        "  nbc=INT                    Absorbing boundary cells. Default: 100.\n"
        "  L=INT                      Absorbing boundary taper length. Default: 30.\n"
        "  alpha=FLOAT                Absorbing boundary strength. Default: 1.\n"
        "  type_compute_Laplace=0|1   Laplacian method: 0=FD8, 1=pseudospectral. Default: 1.\n"
        "  visco=0|1                  Wave mode: 0=acoustic, 1=visco-acoustic/Q path. Default: 0.\n"
        "  qfile=FILE or Qfile=FILE   RSF Q model required when visco=1.\n"
        "  f0=FLOAT                   Reference frequency for visco mode. Default: fdom.\n"
        "  Omega0=FLOAT               Reference angular frequency override. Default: 2*pi*f0.\n"
        "  fc=FLOAT order=INT         Wavenumber filter cutoff/order. Defaults: 120, 2.\n"
        "  flag_smooth=0|1            Smooth extended velocity/Q model. Default: 1.\n"
        "  flag_homo=0|1              Use homogeneous model at source velocity. Default: 0.\n"
        "  progress_interval=INT      Print propagation progress every INT steps.\n"
        "\n"
        "Forward parameters:\n"
        "  sx=FLOAT                   First shot x coordinate. Default: model center.\n"
        "  ns=INT                     Number of shots. Default: 1.\n"
        "  ds=FLOAT                   Shot x spacing. Default: 0.\n"
        "  nr=INT r0=FLOAT dr=FLOAT   Receiver count/origin/spacing; receiver output uses model x grid.\n"
        "\n"
        "Image parameters:\n"
        "  cmp=0|1                    Treat receiver axis as offset from shot when 1. Default: 1.\n"
        "  shot_begin=INT             First shot index to image, zero-based. Alias: first_shot. Default: 0.\n"
        "  shot_end=INT               Last shot index to image, zero-based. Alias: last_shot. Default: last input shot.\n"
        "  max_parallel_shots=INT     Maximum shots imaged concurrently. Aliases: parallel_shots, omp_shots. Default: 1.\n"
        "  in_memory_snapshots=0|1    Store source snapshots in memory instead of temp_img files. Default: 0.\n"
        "  snapshot_interval=INT      Store and correlate one RTM wavefield snapshot every INT time steps. Alias: snapshot_stride. Default: 1.\n"
        "  max_nt=INT                 Truncate input data time samples before imaging.\n"
        "  resample_dt=FLOAT          Resample input data in time before imaging.\n"
        "  compensate=0|1             Enable RTM Q compensation: 0=attenuation/no compensation, 1=Q compensation. Default: 0.\n"
        "  debug_wavefield_prefix=STR Write selected debug wavefields/gathers with this prefix.\n"
        "  debug_snapshot_steps=LIST  Comma-separated or array time-step list for debug snapshots.\n"
        "\n"
        "Examples:\n"
        "  %s mode=forward velocity=vel.rsf qfile=q.rsf output=data.rsf nt=2000 dt=0.001 sx=0 ns=10 ds=25\n"
        "  %s mode=image velocity=vel.rsf qfile=q.rsf seismic_data=data.rsf migration=img.rsf shot_begin=0 shot_end=9 snapshot_interval=5\n",
        name, name, name, name);
}

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
const char* laplace_name(int type_laplace) {
    return type_laplace == 1 ? "pseudospectral" : "FD8";
}

double effective_omega0(float f0, float Omega0) {
    return Omega0 > 0.0f ? Omega0 : 2.0 * M_PI * f0;
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
    if (argc <= 1 || (argc == 2 && (!std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help") || !std::strcmp(argv[1], "help")))) {
        print_usage(argv[0]);
        return 0;
    }
    se_par_init(argc, argv);
    try {
#ifdef _OPENMP
        INFO(("OpenMP enabled, max_threads=%d", omp_get_max_threads()));
#else
        INFO(("OpenMP disabled at compile time"));
#endif
        const char* mode = se_have_par("mode") ? se_get_par_str("mode") : "forward";
        if (!se_have_par("velocity")) ERROR(("Need velocity= RSF velocity model"));
        const char* velocity_path = se_get_par_str("velocity");
        sewave::Grid2D vel = sewave::read_rsf2d(velocity_path);
        auto q_path = [&]() -> const char* {
            if (se_have_par("qfile")) return se_get_par_str("qfile");
            if (se_have_par("Qfile")) return se_get_par_str("Qfile");
            return nullptr;
        };
        auto load_q = [&]() -> sewave::Grid2D {
            const char* path = q_path();
            if (path) return sewave::read_rsf2d(path);
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
            p.visco = se_have_par("visco") ? se_get_par_int("visco") != 0 : false;
            p.type = p.visco ? 1 : 0;
            p.f0 = se_have_par("f0") ? se_get_par_float("f0") : p.fdom;
            p.Omega0 = se_have_par("Omega0") ? se_get_par_float("Omega0") : 0.0f;
            p.fc = se_have_par("fc") ? se_get_par_float("fc") : 120.0f;
            p.order = se_have_par("order") ? se_get_par_int("order") : 2;
            p.flag_smooth = se_have_par("flag_smooth") ? se_get_par_int("flag_smooth") != 0 : true;
            p.flag_homo = se_have_par("flag_homo") ? se_get_par_int("flag_homo") != 0 : false;
            p.q = 1000.0f;
            p.progress_interval = se_have_par("progress_interval") ? se_get_par_int("progress_interval") : (se_have_par("log_interval") ? se_get_par_int("log_interval") : std::max(1, p.nt/10));
            sewave::Grid2D qmodel;
            sewave::Grid2D* qptr = nullptr;
            if (p.visco) { qmodel = load_q(); qptr = &qmodel; }
            const char* output_path = se_get_par_str("output");
            INFO(("forward wave mode: %s (visco=%d)", p.visco ? "visco-acoustic/Q" : "acoustic", p.visco ? 1 : 0));
            INFO(("forward velocity model: %s", velocity_path));
            INFO(("forward Q model: %s", p.visco ? q_path() : "not used (acoustic mode)"));
            INFO(("forward seismic output: %s", output_path));
            INFO(("forward Q compensation/imaging compensation: not applicable in forward mode"));
            INFO(("forward wave-equation parameters: nt=%d dt=%g fdom=%g type_compute_Laplace=%d (%s)",
                  p.nt, p.dt, p.fdom, p.type_compute_laplace, laplace_name(p.type_compute_laplace)));
            INFO(("forward geometry: sx0=%g sz=%g rz=%g ns=%d ds=%g nr=%d r0=%g dr=%g",
                  p.sx, p.sz, p.rz, p.ns, p.ds, p.nr, p.r0, p.dr));
            INFO(("forward boundary/model options: nbc=%d L=%d alpha=%g flag_smooth=%d flag_homo=%d",
                  p.nbc, p.L, p.alpha, p.flag_smooth ? 1 : 0, p.flag_homo ? 1 : 0));
            INFO(("forward visco/filter parameters: f0=%g Omega0_input=%g Omega0_effective=%g fc=%g order=%d",
                  p.f0, p.Omega0, effective_omega0(p.f0, p.Omega0), p.fc, p.order));
            INFO(("forward progress_interval=%d", p.progress_interval));
            INFO(("forward shots: total=%d, calculating shot index range 0-%d", p.ns, std::max(0, p.ns - 1)));
            check_stability(vel, p.dt, p.type_compute_laplace, "forward");
            sewave::Data3D d = sewave::forward(vel, p, qptr);
            sewave::write_rsf3d(output_path, d);
        } else if (mode[0]=='i') {
            if (!se_have_par("seismic_data")) ERROR(("Need seismic_data= input seismic RSF"));
            if (!se_have_par("migration")) ERROR(("Need migration= output image RSF"));
            const char* seismic_data_path = se_get_par_str("seismic_data");
            const char* migration_path = se_get_par_str("migration");
            sewave::Data3D d = sewave::read_rsf3d(seismic_data_path);
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
            p.amp_compensation_sign = p.compensate ? -1 : 1;
            p.q = 1000.0f;
            p.nbc = se_have_par("nbc") ? se_get_par_int("nbc") : 100;
            p.alpha = se_have_par("alpha") ? se_get_par_float("alpha") : 1.0f;
            p.type_compute_laplace = se_have_par("type_compute_Laplace") ? se_get_par_int("type_compute_Laplace") : (se_have_par("type_compute_laplace") ? se_get_par_int("type_compute_laplace") : 1);
            p.L = se_have_par("L") ? se_get_par_int("L") : 30;
            p.visco = se_have_par("visco") ? se_get_par_int("visco") != 0 : false;
            p.type = p.visco ? 1 : 0;
            p.fdom = se_have_par("fdom") ? se_get_par_float("fdom") : 20.0f;
            p.f0 = se_have_par("f0") ? se_get_par_float("f0") : p.fdom;
            p.Omega0 = se_have_par("Omega0") ? se_get_par_float("Omega0") : 0.0f;
            p.fc = se_have_par("fc") ? se_get_par_float("fc") : 120.0f;
            p.order = se_have_par("order") ? se_get_par_int("order") : 2;
            p.flag_smooth = se_have_par("flag_smooth") ? se_get_par_int("flag_smooth") != 0 : true;
            p.flag_homo = se_have_par("flag_homo") ? se_get_par_int("flag_homo") != 0 : false;
            p.sz = se_have_par("sz") ? se_get_par_float("sz") : vel.z0;
            p.rz = se_have_par("rz") ? se_get_par_float("rz") : p.sz;
            p.shot_begin = se_have_par("shot_begin") ? se_get_par_int("shot_begin") : (se_have_par("first_shot") ? se_get_par_int("first_shot") : 0);
            p.shot_end = se_have_par("shot_end") ? se_get_par_int("shot_end") : (se_have_par("last_shot") ? se_get_par_int("last_shot") : -1);
            p.progress_interval = se_have_par("progress_interval") ? se_get_par_int("progress_interval") : (se_have_par("log_interval") ? se_get_par_int("log_interval") : std::max(1, d.nt/10));
            p.max_parallel_shots = se_have_par("max_parallel_shots") ? se_get_par_int("max_parallel_shots") : (se_have_par("parallel_shots") ? se_get_par_int("parallel_shots") : (se_have_par("omp_shots") ? se_get_par_int("omp_shots") : 1));
            p.in_memory_snapshots = se_have_par("in_memory_snapshots") ? se_get_par_int("in_memory_snapshots") != 0 : false;
            p.snapshot_interval = se_have_par("snapshot_interval") ? se_get_par_int("snapshot_interval") : (se_have_par("snapshot_stride") ? se_get_par_int("snapshot_stride") : 1);
            if (p.snapshot_interval < 1) ERROR(("snapshot_interval must be >= 1"));
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
            if (p.visco) { qmodel = load_q(); qptr = &qmodel; }
            INFO(("image wave mode: %s (visco=%d)", p.visco ? "visco-acoustic/Q" : "acoustic", p.visco ? 1 : 0));
            INFO(("image velocity model: %s", velocity_path));
            INFO(("image Q model: %s", p.visco ? q_path() : "not used (acoustic mode)"));
            INFO(("image seismic data: %s", seismic_data_path));
            INFO(("image migration output: %s", migration_path));
            INFO(("image compensation: compensate=%d (%s), internal amp_compensation_sign=%d",
                  p.compensate ? 1 : 0,
                  p.compensate ? "Q compensation" : "Q attenuation/no compensation",
                  p.amp_compensation_sign));
            INFO(("image wave-equation parameters: nt=%d dt=%g fdom=%g type_compute_Laplace=%d (%s)",
                  d.nt, d.dt, p.fdom, p.type_compute_laplace, laplace_name(p.type_compute_laplace)));
            INFO(("image geometry/options: cmp=%d sz=%g rz=%g shot_begin=%d shot_end=%d",
                  p.cmp ? 1 : 0, p.sz, p.rz, p.shot_begin, p.shot_end));
            INFO(("image snapshot/parallel options: snapshot_interval=%d in_memory_snapshots=%d max_parallel_shots=%d",
                  p.snapshot_interval, p.in_memory_snapshots ? 1 : 0, p.max_parallel_shots));
            INFO(("image boundary/model options: nbc=%d L=%d alpha=%g flag_smooth=%d flag_homo=%d",
                  p.nbc, p.L, p.alpha, p.flag_smooth ? 1 : 0, p.flag_homo ? 1 : 0));
            INFO(("image visco/filter parameters: f0=%g Omega0_input=%g Omega0_effective=%g fc=%g order=%d",
                  p.f0, p.Omega0, effective_omega0(p.f0, p.Omega0), p.fc, p.order));
            INFO(("image progress_interval=%d", p.progress_interval));
            {
                int shot0 = std::max(0, p.shot_begin);
                int shot1 = (p.shot_end < 0) ? d.ns - 1 : std::min(d.ns - 1, p.shot_end);
                INFO(("image shots: total=%d, calculating shot index range %d-%d (%d shots)",
                      d.ns, shot0, shot1, shot0 <= shot1 ? shot1 - shot0 + 1 : 0));
            }
            check_stability(vel, d.dt, p.type_compute_laplace, "image");
            auto img = sewave::rtm_image(vel, d, p, qptr);
            sewave::write_rsf2d(migration_path, vel, img, "RTM image");
        } else {
            ERROR(("Unknown mode=%s (use forward or image)", mode));
        }
    } catch (const std::exception& e) {
        ERROR(("sewave2d failed: %s", e.what()));
    }
    se_par_destroy();
    return 0;
}
