#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_par_sep.h>
#include "../SEWAVE/sewave2d.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>

namespace {
int xidx(const sewave::Grid2D& g, float x) { return (int)std::lround((x - g.x0) / g.dx); }
}

int main(int argc, char** argv) {
    se_par_init(argc, argv);
    try {
        const char* vfile = se_have_par("velocity") ? se_get_par_str("velocity") : "synth_vel.rsf";
        const char* dfile = se_have_par("seismic_data") ? se_get_par_str("seismic_data") : "synth_data.rsf";
        sewave::Grid2D vel;
        vel.nz = se_have_par("nz") ? se_get_par_int("nz") : 101;
        vel.nx = se_have_par("nx") ? se_get_par_int("nx") : 201;
        vel.dz = se_have_par("dz") ? se_get_par_float("dz") : 10.0f;
        vel.dx = se_have_par("dx") ? se_get_par_float("dx") : 10.0f;
        vel.z0 = se_have_par("z0") ? se_get_par_float("z0") : 0.0f;
        vel.x0 = se_have_par("x0") ? se_get_par_float("x0") : 0.0f;
        float v0 = se_have_par("v0") ? se_get_par_float("v0") : 2000.0f;
        float v1 = se_have_par("v1") ? se_get_par_float("v1") : 2600.0f;
        float zint = se_have_par("z_interface") ? se_get_par_float("z_interface") : vel.z0 + 0.5f * (vel.nz - 1) * vel.dz;
        vel.v.resize((size_t)vel.nz * vel.nx);
        for (int ix = 0; ix < vel.nx; ++ix) {
            for (int iz = 0; iz < vel.nz; ++iz) {
                float z = vel.z0 + iz * vel.dz;
                vel.v[(size_t)ix * vel.nz + iz] = (z < zint) ? v0 : v1;
            }
        }
        sewave::write_rsf2d(vfile, vel, vel.v, "Velocity");

        sewave::ForwardParams p;
        p.nt = se_have_par("nt") ? se_get_par_int("nt") : 1000;
        p.dt = se_have_par("dt") ? se_get_par_float("dt") : 0.001f;
        p.fdom = se_have_par("fdom") ? se_get_par_float("fdom") : 20.0f;
        p.sx = se_have_par("sx") ? se_get_par_float("sx") : vel.x0 + 0.5f * (vel.nx - 1) * vel.dx;
        p.sz = se_have_par("sz") ? se_get_par_float("sz") : vel.z0;
        p.rz = se_have_par("rz") ? se_get_par_float("rz") : p.sz;
        p.ns = se_have_par("ns") ? se_get_par_int("ns") : 1;
        p.ds = se_have_par("ds") ? se_get_par_float("ds") : 0.0f;
        p.nbc = se_have_par("nbc") ? se_get_par_int("nbc") : 60;
        p.L = se_have_par("L") ? se_get_par_int("L") : 30;
        p.alpha = se_have_par("alpha") ? se_get_par_float("alpha") : 1.0f;
        p.type_compute_laplace = se_have_par("type_compute_Laplace") ? se_get_par_int("type_compute_Laplace") : 0;
        p.type = 0;
        p.visco = false;
        p.flag_smooth = se_have_par("flag_smooth") ? se_get_par_int("flag_smooth") != 0 : false;
        p.progress_interval = se_have_par("progress_interval") ? se_get_par_int("progress_interval") : std::max(1, p.nt / 10);
        sewave::Data3D absdata = sewave::forward(vel, p, nullptr);
        bool subtract_direct = se_have_par("subtract_direct") ? se_get_par_int("subtract_direct") != 0 : true;
        if (subtract_direct) {
            sewave::Grid2D bg = vel;
            std::fill(bg.v.begin(), bg.v.end(), v0);
            sewave::Data3D direct = sewave::forward(bg, p, nullptr);
            for (size_t i = 0; i < absdata.d.size(); ++i) absdata.d[i] -= direct.d[i];
        }

        int cmp = se_have_par("cmp") ? se_get_par_int("cmp") : 0;
        if (!cmp) {
            sewave::write_rsf3d(dfile, absdata);
        } else {
            int nr = se_have_par("nr") ? se_get_par_int("nr") : vel.nx;
            float dh = se_have_par("dh") ? se_get_par_float("dh") : vel.dx;
            float h0 = se_have_par("h0") ? se_get_par_float("h0") : -0.5f * (nr - 1) * dh;
            sewave::Data3D cmpdata;
            cmpdata.nt = absdata.nt; cmpdata.dt = absdata.dt; cmpdata.t0 = absdata.t0;
            cmpdata.nr = nr; cmpdata.r0 = h0; cmpdata.dr = dh;
            cmpdata.ns = absdata.ns; cmpdata.s0 = absdata.s0; cmpdata.ds = absdata.ds;
            cmpdata.d.assign((size_t)cmpdata.ns * cmpdata.nr * cmpdata.nt, 0.0f);
            for (int is = 0; is < cmpdata.ns; ++is) {
                float sx = cmpdata.s0 + is * cmpdata.ds;
                for (int ir = 0; ir < cmpdata.nr; ++ir) {
                    float rx = sx + cmpdata.r0 + ir * cmpdata.dr;
                    int ix = xidx(vel, rx);
                    if (ix < 0 || ix >= vel.nx) continue;
                    std::memcpy(&cmpdata.d[((size_t)is * cmpdata.nr + ir) * cmpdata.nt],
                                &absdata.d[((size_t)is * absdata.nr + ix) * absdata.nt],
                                (size_t)cmpdata.nt * sizeof(float));
                }
            }
            sewave::write_rsf3d(dfile, cmpdata);
        }
        INFO(("Wrote synthetic velocity=%s seismic_data=%s cmp=%d subtract_direct=%d", vfile, dfile, cmp, subtract_direct ? 1 : 0));
    } catch (const std::exception& e) {
        ERROR(("sewave2d_synth failed: %s", e.what()));
        se_par_destroy();
        return 1;
    }
    se_par_destroy();
    return 0;
}
