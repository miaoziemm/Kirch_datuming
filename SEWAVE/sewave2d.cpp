#include "sewave2d.h"
#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <fftw3.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sewave
{
    namespace
    {
        constexpr float PI = 3.14159265358979323846f;
        inline int idx(int iz, int ix, int nz) { return ix * nz + iz; }
        float ricker(float t, float f)
        {
            float tau = t - 1.5f / f;
            float a = PI * f * tau;
            return (1.0f - 2.0f * a * a) * std::exp(-a * a);
        }
        int clampi(int v, int lo, int hi) { return std::max(lo, std::min(hi, v)); }
        int xindex(const Grid2D &g, float x) { return clampi((int)std::lround((x - g.x0) / g.dx), 0, g.nx - 1); }
        int xindex_unclamped(const Grid2D &g, float x) { return (int)std::lround((x - g.x0) / g.dx); }
        int zindex(const Grid2D &g, float z) { return clampi((int)std::lround((z - g.z0) / g.dz), 0, g.nz - 1); }
        int zindex_unclamped(const Grid2D &g, float z) { return (int)std::lround((z - g.z0) / g.dz); }
        void check_same_grid(const Grid2D &a, const Grid2D &b, const char *name)
        {
            if (a.nz != b.nz || a.nx != b.nx || a.dz != b.dz || a.dx != b.dx || a.z0 != b.z0 || a.x0 != b.x0)
                throw std::runtime_error(std::string(name) + " grid must match velocity grid");
        }

        void tripd(std::vector<float> &d, std::vector<float> &e, std::vector<float> &b, int n)
        {
            for (int k = 1; k < n; ++k)
            {
                float temp = e[k - 1];
                e[k - 1] = temp / d[k - 1];
                d[k] -= temp * e[k - 1];
            }
            for (int k = 1; k < n; ++k)
                b[k] -= e[k - 1] * b[k - 1];
            b[n - 1] /= d[n - 1];
            for (int k = n - 1; k > 0; --k)
                b[k - 1] = b[k - 1] / d[k - 1] - e[k - 1] * b[k];
        }
        void smooth2d_flat(int n1, int n2, float r1, float r2, std::vector<float> &v)
        {
            int nmax = std::max(n1, n2);
            std::vector<float> w((size_t)n1 * n2, 1.0f), d(nmax), e(nmax), f(nmax);
            r1 = r1 * r1 * 0.25f;
            r2 = r2 * r2 * 0.25f;
            for (int iz = 0; iz < n1; ++iz)
            {
                for (int ix = 0; ix < n2 - 1; ++ix)
                {
                    d[ix] = 1.0f + r2 * (w[idx(iz, ix, n1)] + w[idx(iz, ix + 1, n1)]);
                    e[ix] = -r2 * w[idx(iz, ix + 1, n1)];
                    f[ix] = v[idx(iz, ix, n1)];
                }
                d[0] -= r2 * w[idx(iz, 0, n1)];
                d[n2 - 1] = 1.0f + r2 * w[idx(iz, n2 - 1, n1)];
                f[n2 - 1] = v[idx(iz, n2 - 1, n1)];
                tripd(d, e, f, n2);
                for (int ix = 0; ix < n2; ++ix)
                    v[idx(iz, ix, n1)] = f[ix];
            }
            for (int ix = 0; ix < n2; ++ix)
            {
                for (int iz = 0; iz < n1 - 2; ++iz)
                {
                    d[iz] = 1.0f + r1 * (w[idx(iz + 1, ix, n1)] + w[idx(iz + 2, ix, n1)]);
                    e[iz] = -r1 * w[idx(iz + 2, ix, n1)];
                    f[iz] = v[idx(iz + 1, ix, n1)];
                }
                f[0] += r1 * w[idx(1, ix, n1)] * v[idx(0, ix, n1)];
                d[n1 - 2] = 1.0f + r1 * w[idx(n1 - 1, ix, n1)];
                f[n1 - 2] = v[idx(n1 - 1, ix, n1)];
                tripd(d, e, f, n1 - 1);
                for (int iz = 0; iz < n1 - 1; ++iz)
                    v[idx(iz + 1, ix, n1)] = f[iz];
            }
        }
        struct ExtModel
        {
            int nx = 0, nz = 0, nbc = 0, L = 30;
            float dx = 1, dz = 1, alpha = 1, Omega0 = 0, fc = 120;
            int order = 2, type = 0;
            std::vector<float> vp, Q;
        };

        template <class P>
        ExtModel make_ext(const Grid2D &vel, const Grid2D *q, const P &p)
        {
            if (p.visco && !q)
                throw std::runtime_error("visco=1 requires qfile= RSF Q model");
            if (q)
                check_same_grid(vel, *q, "Q model");
            ExtModel e;
            e.nbc = p.nbc;
            e.L = p.L;
            e.nx = vel.nx + 2 * p.nbc;
            e.nz = vel.nz + 2 * p.nbc;
            e.dx = vel.dx;
            e.dz = vel.dz;
            e.alpha = p.alpha;
            e.Omega0 = (p.Omega0 > 0 ? p.Omega0 : 2 * PI * p.f0);
            e.fc = p.fc;
            e.order = p.order;
            e.type = p.visco ? 1 : 0;
            e.vp.assign((size_t)e.nx * e.nz, 0);
            e.Q.assign((size_t)e.nx * e.nz, std::max(1.0f, p.q));
            for (int ix = 0; ix < e.nx; ++ix)
            {
                int sx = clampi(ix - p.nbc, 0, vel.nx - 1);
                for (int iz = 0; iz < e.nz; ++iz)
                {
                    int sz = clampi(iz - p.nbc, 0, vel.nz - 1);
                    e.vp[idx(iz, ix, e.nz)] = vel.v[idx(sz, sx, vel.nz)];
                    e.Q[idx(iz, ix, e.nz)] = q ? q->v[idx(sz, sx, vel.nz)] : p.q;
                }
            }
            if (p.flag_smooth)
            {
                smooth2d_flat(e.nz, e.nx, 5, 5, e.vp);
                smooth2d_flat(e.nz, e.nx, 5, 5, e.Q);
            }
            return e;
        }

        void apply_homo(ExtModel &e, int sx, int sz)
        {
            int is = idx(clampi(sz, 0, e.nz - 1), clampi(sx, 0, e.nx - 1), e.nz);
            std::fill(e.vp.begin(), e.vp.end(), e.vp[is]);
            std::fill(e.Q.begin(), e.Q.end(), e.Q[is]);
        }

        void wavenumbers(const ExtModel &e, std::vector<float> &kx, std::vector<float> &kz)
        {
            kx.resize(e.nx);
            kz.resize(e.nz);
            for (int i = 0; i < e.nx; ++i)
                kx[i] = 2 * PI * (i <= e.nx / 2 ? i : i - e.nx) / (e.nx * e.dx);
            for (int i = 0; i < e.nz; ++i)
                kz[i] = 2 * PI * (i <= e.nz / 2 ? i : i - e.nz) / (e.nz * e.dz);
        }
        void butter(const ExtModel &e, const std::vector<float> &kx, const std::vector<float> &kz, std::vector<float> &H)
        {
            H.assign((size_t)e.nx * e.nz, 1.0f);
            float vmax = *std::max_element(e.vp.begin(), e.vp.end());
            float kc = e.fc * 2 * PI / std::max(1.0f, vmax);
            for (int ix = 0; ix < e.nx; ++ix)
                for (int iz = 0; iz < e.nz; ++iz)
                {
                    float k = std::sqrt(kx[ix] * kx[ix] + kz[iz] * kz[iz]);
                    H[idx(iz, ix, e.nz)] = 1.0f / (1.0f + std::pow(k / std::max(kc, 1e-6f), 2 * e.order));
                }
        }

        void acoustic_fd8_step(
            const ExtModel &e,
            const std::vector<float> &cur,
            const std::vector<float> &old,
            const std::vector<float> &velocity_dt2,
            std::vector<float> &next)
        {
            const int nx = e.nx;
            const int nz = e.nz;
            const int n = nx * nz;
            const float dx2 = e.dx * e.dx;
            const float dz2 = e.dz * e.dz;

            // The FD8 stencil is not evaluated in the outer four samples.
            // Preserve the original zero-Laplacian update there; absorb()
            // subsequently replaces the absorbing-boundary samples.
            for (int i = 0; i < n; ++i)
                next[i] = 2.0f * cur[i] - old[i];

            for (int ix = 4; ix < nx - 4; ++ix)
            {
                for (int iz = 4; iz < nz - 4; ++iz)
                {
                    const int i = idx(iz, ix, nz);
                    const float lap =
                        ((-1.f / 560 * (cur[idx(iz, ix - 4, nz)] + cur[idx(iz, ix + 4, nz)]) +
                          8.f / 315 * (cur[idx(iz, ix - 3, nz)] + cur[idx(iz, ix + 3, nz)]) -
                          1.f / 5 * (cur[idx(iz, ix - 2, nz)] + cur[idx(iz, ix + 2, nz)]) +
                          8.f / 5 * (cur[idx(iz, ix - 1, nz)] + cur[idx(iz, ix + 1, nz)]) -
                          205.f / 72 * cur[i]) /
                         dx2) +
                        ((-1.f / 560 * (cur[idx(iz - 4, ix, nz)] + cur[idx(iz + 4, ix, nz)]) +
                          8.f / 315 * (cur[idx(iz - 3, ix, nz)] + cur[idx(iz + 3, ix, nz)]) -
                          1.f / 5 * (cur[idx(iz - 2, ix, nz)] + cur[idx(iz + 2, ix, nz)]) +
                          8.f / 5 * (cur[idx(iz - 1, ix, nz)] + cur[idx(iz + 1, ix, nz)]) -
                          205.f / 72 * cur[i]) /
                         dz2);
                    next[i] = 2.0f * cur[i] - old[i] + velocity_dt2[i] * lap;
                }
            }
        }

        void fft_ops(const ExtModel &e, const std::vector<float> &cur, const std::vector<float> &old, float dt, int type_lap, bool filt, const std::vector<float> &H, std::vector<float> &lap, std::vector<float> &amp, std::vector<float> &KP)
        {
            int nx = e.nx, nz = e.nz, nc = nz / 2 + 1, n = nx * nz;
            if (type_lap != 0 && type_lap != 1)
                throw std::runtime_error("type_compute_Laplace must be 0 (FD8) or 1 (PS)");
            lap.assign(n, 0);
            amp.assign(n, 0);
            KP.assign(n, 0);
            if (type_lap == 0)
            {
                for (int ix = 4; ix < nx - 4; ++ix)
                    for (int iz = 4; iz < nz - 4; ++iz)
                    {
                        int i = idx(iz, ix, nz);
                        float dx2 = e.dx * e.dx, dz2 = e.dz * e.dz;
                        lap[i] = ((-1.f / 560 * (cur[idx(iz, ix - 4, nz)] + cur[idx(iz, ix + 4, nz)]) + 8.f / 315 * (cur[idx(iz, ix - 3, nz)] + cur[idx(iz, ix + 3, nz)]) - 1.f / 5 * (cur[idx(iz, ix - 2, nz)] + cur[idx(iz, ix + 2, nz)]) + 8.f / 5 * (cur[idx(iz, ix - 1, nz)] + cur[idx(iz, ix + 1, nz)]) - 205.f / 72 * cur[i]) / dx2) + ((-1.f / 560 * (cur[idx(iz - 4, ix, nz)] + cur[idx(iz + 4, ix, nz)]) + 8.f / 315 * (cur[idx(iz - 3, ix, nz)] + cur[idx(iz + 3, ix, nz)]) - 1.f / 5 * (cur[idx(iz - 2, ix, nz)] + cur[idx(iz + 2, ix, nz)]) + 8.f / 5 * (cur[idx(iz - 1, ix, nz)] + cur[idx(iz + 1, ix, nz)]) - 205.f / 72 * cur[i]) / dz2);
                    }
            }
            std::vector<float> p1t(n), in(cur), in1(n), out(n), outa(n), outk(n);
            for (int i = 0; i < n; ++i)
            {
                p1t[i] = (cur[i] - old[i]) / dt;
                in1[i] = p1t[i];
            }
            std::vector<fftwf_complex> P((size_t)nx * nc), P1((size_t)nx * nc);
            fftwf_plan fp, fp1, bl, ba, bk;
#pragma omp critical(sewave_fftw_plan)
            {
                fp = fftwf_plan_dft_r2c_2d(nx, nz, in.data(), P.data(), FFTW_ESTIMATE);
                fp1 = fftwf_plan_dft_r2c_2d(nx, nz, in1.data(), P1.data(), FFTW_ESTIMATE);
                bl = fftwf_plan_dft_c2r_2d(nx, nz, P.data(), out.data(), FFTW_ESTIMATE);
                ba = fftwf_plan_dft_c2r_2d(nx, nz, P1.data(), outa.data(), FFTW_ESTIMATE);
                bk = fftwf_plan_dft_c2r_2d(nx, nz, P.data(), outk.data(), FFTW_ESTIMATE);
            }
            float scale = 1.0f / n;
            if (type_lap == 1)
            {
                fftwf_execute(fp);
                for (int ix = 0; ix < nx; ++ix)
                {
                    float kx = 2 * PI * (ix <= nx / 2 ? ix : ix - nx) / (nx * e.dx);
                    for (int iz = 0; iz < nc; ++iz)
                    {
                        float kz = 2 * PI * iz / (nz * e.dz), k2 = kx * kx + kz * kz;
                        size_t si = (size_t)ix * nc + iz;
                        P[si][0] *= -k2;
                        P[si][1] *= -k2;
                    }
                }
                fftwf_execute(bl);
                for (int i = 0; i < n; ++i)
                    lap[i] = out[i] * scale;
            }
            fftwf_execute(fp);
            fftwf_execute(fp1);
            for (int ix = 0; ix < nx; ++ix)
            {
                float kx = 2 * PI * (ix <= nx / 2 ? ix : ix - nx) / (nx * e.dx);
                for (int iz = 0; iz < nc; ++iz)
                {
                    float kz = 2 * PI * iz / (nz * e.dz), k = std::sqrt(kx * kx + kz * kz);
                    size_t si = (size_t)ix * nc + iz;
                    float h = (!H.empty() && filt) ? H[idx(iz, ix, nz)] : 1.0f;
                    P[si][0] *= k * h;
                    P[si][1] *= k * h;
                    P1[si][0] *= k * h;
                    P1[si][1] *= k * h;
                }
            }
            fftwf_execute(ba);
            fftwf_execute(bk);
            for (int i = 0; i < n; ++i)
            {
                amp[i] = outa[i] * scale;
                KP[i] = outk[i] * scale;
            }
#pragma omp critical(sewave_fftw_plan)
            {
                fftwf_destroy_plan(fp);
                fftwf_destroy_plan(fp1);
                fftwf_destroy_plan(bl);
                fftwf_destroy_plan(ba);
                fftwf_destroy_plan(bk);
            }
        }

        void absorb(const ExtModel &e, std::vector<float> &next, const std::vector<float> &cur, const std::vector<float> &old, float dt)
        {
            int L = std::min(e.L, std::min(e.nx, e.nz) / 2);
            if (L < 5)
                return;
            auto upd = [&](int ix, int iz, int dx, int dz, float w)
            { int i=idx(iz,ix,e.nz); float s=e.alpha*e.vp[i]*dt/std::max(e.dx,e.dz); float t1=(2-s)*(1-s)/2,t2=s*(2-s),t3=s*(s-1)/2; int i1=idx(iz+dz,ix+dx,e.nz),i2=idx(iz+2*dz,ix+2*dx,e.nz),i3=idx(iz+3*dz,ix+3*dx,e.nz),i4=idx(iz+4*dz,ix+4*dx,e.nz); float ow=2*(t1*cur[i]+t2*cur[i1]+t3*cur[i2])-(t1*t1*old[i]+2*t1*t2*old[i1]+(2*t1*t3+t2*t2)*old[i2]+2*t2*t3*old[i3]+t3*t3*old[i4]); next[i]=w*ow+(1-w)*next[i]; };
            for (int ix = 0; ix < e.nx; ++ix)
                for (int iz = 0; iz < e.nz; ++iz)
                {
                    if (ix < L - 4)
                        upd(ix, iz, 1, 0, 1.f - ix / (float)L);
                    if (ix > e.nx - L + 3)
                        upd(ix, iz, -1, 0, 1.f - (e.nx - ix - 1) / (float)L);
                    if (iz < L - 4)
                        upd(ix, iz, 0, 1, 1.f - iz / (float)L);
                    if (iz > e.nz - L + 3)
                        upd(ix, iz, 0, -1, 1.f - (e.nz - iz - 1) / (float)L);
                }
        }

        void write_debug_gather(const std::vector<float> &gather, int nx, int nt, float x0, float dx, float dt, const std::string &path)
        {
            sep_t *out = sep_open(path.c_str(), SEP_WRITE, 0);
            if (!out)
                throw std::runtime_error("cannot create debug receiver gather " + path);
            out->headers->ndim = 2;
            out->headers->n[0] = nt;
            out->headers->n[1] = nx;
            out->headers->o[0] = 0.0f;
            out->headers->o[1] = x0;
            out->headers->d[0] = dt;
            out->headers->d[1] = dx;
            sep_set_header(out, "label1", "Time");
            sep_set_header(out, "label2", "ReceiverX");
            sep_set_header(out, "label", "Receiver gather used for RTM");
            se_fsio_write_float(out->data->io, const_cast<float *>(gather.data()), gather.size());
            sep_write_headers(out);
            sep_close(out);
        }

        void report_gather_energy(const std::vector<float> &gather, int nx, int nt, int shot, const std::string &prefix, float x0, float dx, float dt)
        {
            float max_all = 0.0f, max_reverse_head = 0.0f;
            int nhead = std::min(nt, 4000);
            for (size_t i = 0; i < gather.size(); ++i)
                max_all = std::max(max_all, std::fabs(gather[i]));
            for (int ix = 0; ix < nx; ++ix)
                for (int it = nt - nhead; it < nt; ++it)
                    max_reverse_head = std::max(max_reverse_head, std::fabs(gather[(size_t)ix * nt + it]));
            INFO(("image shot %d: receiver gather maxabs=%g, maxabs in first %d reverse steps (last data samples)=%g", shot, max_all, nhead, max_reverse_head));
            if (!prefix.empty())
            {
                char name[512];
                std::snprintf(name, sizeof(name), "%s_shot%d_gather.rsf", prefix.c_str(), shot);
                write_debug_gather(gather, nx, nt, x0, dx, dt, name);
                INFO(("Wrote debug receiver gather %s", name));
            }
        }

        void write_debug_snapshot(const ExtModel &e, const Grid2D &g, const std::vector<float> &wave, const std::string &path)
        {
            std::vector<float> crop((size_t)g.nz * g.nx, 0.0f);
            for (int ix = 0; ix < g.nx; ++ix)
                for (int iz = 0; iz < g.nz; ++iz)
                    crop[idx(iz, ix, g.nz)] = wave[idx(iz + e.nbc, ix + e.nbc, e.nz)];
            sep_t *out = sep_open(path.c_str(), SEP_WRITE, 0);
            if (!out)
                throw std::runtime_error("cannot create debug wavefield snapshot " + path);
            out->headers->ndim = 2;
            out->headers->n[0] = g.nz;
            out->headers->n[1] = g.nx;
            out->headers->o[0] = g.z0;
            out->headers->o[1] = g.x0;
            out->headers->d[0] = g.dz;
            out->headers->d[1] = g.dx;
            sep_set_header(out, "label1", "Depth");
            sep_set_header(out, "label2", "Lateral");
            sep_set_header(out, "label", "Wavefield");
            se_fsio_write_float(out->data->io, crop.data(), crop.size());
            sep_write_headers(out);
            sep_close(out);
        }

        void maybe_write_debug_snapshot(const ExtModel &e, const Grid2D *g, const std::vector<float> &wave, int step, const std::vector<int> *steps, const std::string *prefix, const char *kind, int shot)
        {
            if (!g || !steps || !prefix || prefix->empty() || !kind)
                return;
            if (std::find(steps->begin(), steps->end(), step) == steps->end())
                return;
            char name[512];
            std::snprintf(name, sizeof(name), "%s_shot%d_%s_it%d.rsf", prefix->c_str(), shot, kind, step);
            write_debug_snapshot(e, *g, wave, name);
            INFO(("Wrote debug wavefield snapshot %s", name));
        }
        void ensure_directory(const char *dir, const char *purpose)
        {
            if (::mkdir(dir, 0777) != 0 && errno != EEXIST)
                throw std::runtime_error(std::string("cannot create ") + dir + " directory for " + purpose);
        }

        std::string shot_image_path(int shot)
        {
            char name[256];
            std::snprintf(name, sizeof(name), "temp_shot_img/shot_%06d.rsf", shot);
            return name;
        }

        void write_shot_image_rsf(const std::string &path, const Grid2D &g, const std::vector<float> &data, int shot, float shot_x, float shot_z)
        {
            sep_t *s = sep_open(path.c_str(), SEP_WRITE, 0);
            if (!s)
                throw std::runtime_error("cannot create per-shot image " + path);
            s->headers->ndim = 2;
            s->headers->n[0] = g.nz;
            s->headers->n[1] = g.nx;
            s->headers->o[0] = g.z0;
            s->headers->o[1] = g.x0;
            s->headers->d[0] = g.dz;
            s->headers->d[1] = g.dx;
            sep_set_header(s, "label1", "Depth");
            sep_set_header(s, "label2", "Lateral");
            sep_set_header(s, "label", "Per-shot RTM Image");
            sep_set_header_int(s, "shot_index", shot);
            sep_set_header_float(s, "shot_x", shot_x);
            sep_set_header_float(s, "shot_z", shot_z);
            sep_set_header_float(s, "sx", shot_x);
            sep_set_header_float(s, "sz", shot_z);
            se_fsio_write_float(s->data->io, const_cast<float *>(data.data()), data.size());
            sep_write_headers(s);
            sep_close(s);
        }

        FILE *open_snapshot_temp(int shot, std::string &path)
        {
            ensure_directory("temp_img", "source snapshots");
            char name[256];
            std::snprintf(name, sizeof(name), "temp_img/sewave_snap_%ld_shot_%d_XXXXXX", (long)::getpid(), shot);
            int fd = ::mkstemp(name);
            if (fd < 0)
                throw std::runtime_error("cannot create source snapshot temporary file in temp_img");
            FILE *fp = ::fdopen(fd, "w+b");
            if (!fp)
            {
                ::close(fd);
                ::unlink(name);
                throw std::runtime_error("cannot open source snapshot temporary file stream");
            }
            path = name;
            return fp;
        }

        int count_stored_snapshots(int nt, int interval, int offset)
        {
            if (interval <= 1)
                return nt;
            if (offset < 0)
                offset = 0;
            if (offset >= interval)
                offset %= interval;
            if (offset >= nt)
                return 0;
            return (nt - 1 - offset) / interval + 1;
        }

        bool store_snapshot_step(int it, int interval, int offset)
        {
            if (interval <= 1)
                return true;
            int r = (it - offset) % interval;
            return r == 0;
        }

        int stored_snapshot_index(int it, int interval, int offset)
        {
            if (interval <= 1)
                return it;
            return (it - offset) / interval;
        }

        void log_progress(const char *label, int it, int nt)
        {
#pragma omp critical(sewave_progress_log)
            {
                INFO(("%s: step %d/%d", label, it, nt));
            }
        }

        void propagate(const ExtModel &e, int nt, float dt, int type_lap, bool filter, int amp_sign,
                       int srcx, int srcz, int record_z, bool source_add, int progress_interval,
                       const char *progress_label, const std::vector<float> *inject,
                       std::vector<float> *record, std::vector<float> *snaps, std::vector<float> *image,
                       const std::vector<float> *srcsnaps, FILE *snap_write = nullptr, FILE *snap_read = nullptr,
                       int snapshot_interval = 1, int snapshot_offset = 0, int src_snapshot_interval = 1, int src_snapshot_offset = 0,
                       const Grid2D *debug_grid = nullptr, const std::vector<int> *debug_steps = nullptr,
                       const std::string *debug_prefix = nullptr, const char *debug_kind = nullptr, int debug_shot = -1)
        {
            int n = e.nx * e.nz;
            std::vector<float> old(n, 0), cur(n, 0), nxt(n, 0), lap, amp, KP, H, snapbuf;
            const bool acoustic_fd8 = (e.type == 0 && type_lap == 0);
            std::vector<float> velocity_dt2;
            if (acoustic_fd8)
            {
                velocity_dt2.resize(n);
                for (int i = 0; i < n; ++i)
                    velocity_dt2[i] = e.vp[i] * e.vp[i] * dt * dt;
            }
            if (image && snap_read)
                snapbuf.resize(n);
            std::vector<float> kx, kz;
            if (!acoustic_fd8)
            {
                wavenumbers(e, kx, kz);
                if (filter)
                    butter(e, kx, kz, H);
            }
            const bool receiver_gather = inject && inject->size() != (size_t)nt;
            for (int it = 0; it < nt; ++it)
            {
                if (progress_interval > 0 && (it == 0 || (it + 1) % progress_interval == 0 || it == nt - 1))
                    log_progress(progress_label, it + 1, nt);
                if (inject && !receiver_gather)
                {
                    if (source_add)
                        cur[idx(srcz, srcx, e.nz)] += (*inject)[it];
                    else
                        cur[idx(srcz, srcx, e.nz)] = (*inject)[it];
                }
                else if (receiver_gather)
                {
                    int rit = nt - it - 1;
                    for (int ix = e.nbc; ix < e.nx - e.nbc; ++ix)
                        cur[idx(srcz, ix, e.nz)] = (*inject)[(size_t)(ix - e.nbc) * nt + rit];
                }
                if (acoustic_fd8)
                {
                    acoustic_fd8_step(e, cur, old, velocity_dt2, nxt);
                }
                else
                {
                    fft_ops(e, cur, old, dt, type_lap, filter, H, lap, amp, KP);
                }
                const int general_update_x_begin = acoustic_fd8 ? e.nx : 0;
                for (int ix = general_update_x_begin; ix < e.nx; ++ix)
                    for (int iz = 0; iz < e.nz; ++iz)
                    {
                        int i = idx(iz, ix, e.nz);
                        if (e.type == 1)
                        {
                            float c0 = 411.7f + 4.36f, c1 = -51.64f, c2 = 4.366f / 2.0f, Q = std::max(1.0f, e.Q[i]), v = e.vp[i];
                            float C1 = (1 - 4 * c2 / (PI * Q) + 2 * std::log(e.Omega0) / (PI * Q)) / (v * v);
                            float C2 = 2 * c1 / (PI * Q * v);
                            float C3 = 2 * (c0 - c2) / (PI * Q * v * v);
                            float C4 = 1 / (v * Q);
                            nxt[i] = 2 * cur[i] - old[i] + (dt * dt / C1) * (C4 * amp_sign * (-amp[i]) + lap[i] - C2 * KP[i] - C3 * cur[i]);
                        }
                        else
                        {
                            nxt[i] = 2 * cur[i] - old[i] + e.vp[i] * e.vp[i] * dt * dt * lap[i];
                        }
                    }
                absorb(e, nxt, cur, old, dt);
                old.swap(cur);
                cur.swap(nxt);
                maybe_write_debug_snapshot(e, debug_grid, cur, it + 1, debug_steps, debug_prefix, debug_kind, debug_shot);
                if (record)
                {
                    for (int ix = e.nbc; ix < e.nx - e.nbc; ++ix)
                        (*record)[(size_t)(ix - e.nbc) * nt + it] = cur[idx(record_z, ix, e.nz)];
                }
                if (snaps && store_snapshot_step(it, snapshot_interval, snapshot_offset))
                {
                    int si = stored_snapshot_index(it, snapshot_interval, snapshot_offset);
                    std::memcpy(&(*snaps)[(size_t)si * n], cur.data(), n * sizeof(float));
                }
                if (snap_write && store_snapshot_step(it, snapshot_interval, snapshot_offset))
                {
                    if (std::fwrite(cur.data(), sizeof(float), n, snap_write) != (size_t)n)
                        throw std::runtime_error("failed writing source snapshot temporary file");
                }
                if (image && (srcsnaps || snap_read))
                {
                    int rit = nt - it - 1;
                    const float *srcptr = nullptr;
                    int rsi = stored_snapshot_index(rit, src_snapshot_interval, src_snapshot_offset);
                    if (srcsnaps)
                        srcptr = &(*srcsnaps)[(size_t)rsi * n];
                    else
                    {
                        if (fseeko(snap_read, (off_t)rsi * n * sizeof(float), SEEK_SET) != 0 || std::fread(snapbuf.data(), sizeof(float), n, snap_read) != (size_t)n)
                            throw std::runtime_error("failed reading source snapshot temporary file");
                        srcptr = snapbuf.data();
                    }
                    for (int ix = 0; ix < e.nx - 2 * e.nbc; ++ix)
                        for (int iz = 0; iz < e.nz - 2 * e.nbc; ++iz)
                        {
                            (*image)[idx(iz, ix, e.nz - 2 * e.nbc)] += cur[idx(iz + e.nbc, ix + e.nbc, e.nz)] * srcptr[idx(iz + e.nbc, ix + e.nbc, e.nz)];
                        }
                }
            }
        }

        void correlate_snapshots(const ExtModel &e, int nt, int interval, int src_offset, int rec_offset, const std::vector<float> *srcsnaps, FILE *srcfile, FILE *recfile, std::vector<float> &image)
        {
            int n = e.nx * e.nz;
            std::vector<float> srcbuf, recbuf(n);
            if (srcfile)
                srcbuf.resize(n);
            for (int it = 0; it < nt; ++it)
            {
                if (!store_snapshot_step(it, interval, rec_offset))
                    continue;
                int sit = nt - it - 1;
                int ri = stored_snapshot_index(it, interval, rec_offset);
                if (fseeko(recfile, (off_t)ri * n * sizeof(float), SEEK_SET) != 0 || std::fread(recbuf.data(), sizeof(float), n, recfile) != (size_t)n)
                    throw std::runtime_error("failed reading receiver snapshot temporary file");
                const float *srcptr = nullptr;
                int si = stored_snapshot_index(sit, interval, src_offset);
                if (srcsnaps)
                    srcptr = &(*srcsnaps)[(size_t)si * n];
                else
                {
                    if (fseeko(srcfile, (off_t)si * n * sizeof(float), SEEK_SET) != 0 || std::fread(srcbuf.data(), sizeof(float), n, srcfile) != (size_t)n)
                        throw std::runtime_error("failed reading source snapshot temporary file");
                    srcptr = srcbuf.data();
                }
                for (int ix = 0; ix < e.nx - 2 * e.nbc; ++ix)
                    for (int iz = 0; iz < e.nz - 2 * e.nbc; ++iz)
                    {
                        image[idx(iz, ix, e.nz - 2 * e.nbc)] += recbuf[idx(iz + e.nbc, ix + e.nbc, e.nz)] * srcptr[idx(iz + e.nbc, ix + e.nbc, e.nz)];
                    }
            }
        }
    } // namespace

    Grid2D read_rsf2d(const std::string &path)
    {
        sep_t *s = sep_open(path.c_str(), SEP_READ, 0);
        if (!s)
            throw std::runtime_error("cannot open " + path);
        if (s->headers->ndim < 2)
            throw std::runtime_error("need 2D RSF " + path);
        Grid2D g;
        g.nz = s->headers->n[0];
        g.nx = s->headers->n[1];
        g.z0 = s->headers->o[0];
        g.dz = s->headers->d[0];
        g.x0 = s->headers->o[1];
        g.dx = s->headers->d[1];
        g.v.resize((size_t)g.nz * g.nx);
        se_fsio_seek(s->data->io, 0);
        se_fsio_read_float(s->data->io, g.v.data(), g.v.size());
        sep_close(s);
        return g;
    }
    Data3D read_rsf3d(const std::string &path)
    {
        sep_t *s = sep_open(path.c_str(), SEP_READ, 0);
        if (!s)
            throw std::runtime_error("cannot open " + path);
        if (s->headers->ndim < 3)
            throw std::runtime_error("need 3D RSF " + path);
        Data3D d;
        d.nt = s->headers->n[0];
        d.nr = s->headers->n[1];
        d.ns = s->headers->n[2];
        d.t0 = s->headers->o[0];
        d.dt = s->headers->d[0];
        d.r0 = s->headers->o[1];
        d.dr = s->headers->d[1];
        d.s0 = s->headers->o[2];
        d.ds = s->headers->d[2];
        d.d.resize((size_t)d.nt * d.nr * d.ns);
        se_fsio_seek(s->data->io, 0);
        se_fsio_read_float(s->data->io, d.d.data(), d.d.size());
        sep_close(s);
        return d;
    }
    void write_rsf2d(const std::string &path, const Grid2D &g, const std::vector<float> &data, const char *label)
    {
        sep_t *s = sep_open(path.c_str(), SEP_WRITE, 0);
        s->headers->ndim = 2;
        s->headers->n[0] = g.nz;
        s->headers->n[1] = g.nx;
        s->headers->o[0] = g.z0;
        s->headers->o[1] = g.x0;
        s->headers->d[0] = g.dz;
        s->headers->d[1] = g.dx;
        sep_set_header(s, "label1", "Depth");
        sep_set_header(s, "label2", "Lateral");
        sep_set_header(s, "label", label);
        se_fsio_write_float(s->data->io, const_cast<float *>(data.data()), data.size());
        sep_write_headers(s);
        sep_close(s);
    }
    void write_rsf3d(const std::string &path, const Data3D &d)
    {
        sep_t *s = sep_open(path.c_str(), SEP_WRITE, 0);
        s->headers->ndim = 3;
        s->headers->n[0] = d.nt;
        s->headers->n[1] = d.nr;
        s->headers->n[2] = d.ns;
        s->headers->o[0] = d.t0;
        s->headers->o[1] = d.r0;
        s->headers->o[2] = d.s0;
        s->headers->d[0] = d.dt;
        s->headers->d[1] = d.dr;
        s->headers->d[2] = d.ds;
        sep_set_header(s, "label1", "Time");
        sep_set_header(s, "label2", "Offset/Receiver");
        sep_set_header(s, "label3", "Shot");
        se_fsio_write_float(s->data->io, const_cast<float *>(d.d.data()), d.d.size());
        sep_write_headers(s);
        sep_close(s);
    }

    Data3D forward(const Grid2D &vel, const ForwardParams &p, const Grid2D *qmodel)
    {
        ExtModel ebase = make_ext(vel, qmodel, p);
        Data3D out;
        out.nt = p.nt;
        out.dt = p.dt;
        out.nr = vel.nx;
        out.r0 = vel.x0;
        out.dr = vel.dx;
        out.ns = p.ns;
        out.s0 = p.sx;
        out.ds = p.ds;
        out.d.assign((size_t)out.nt * out.nr * out.ns, 0);
#pragma omp parallel for schedule(dynamic)
        for (int is = 0; is < out.ns; ++is)
        {
            std::vector<float> src(out.nt), rec((size_t)out.nr * out.nt);
            char label[64];
            std::snprintf(label, sizeof(label), "forward shot %d/%d", is + 1, out.ns);
            for (int it = 0; it < out.nt; ++it)
                src[it] = ricker(it * out.dt, p.fdom);
            ExtModel e = ebase;
            int sx = xindex(vel, p.sx + is * p.ds) + e.nbc, sz = zindex(vel, p.sz) + e.nbc, rz = zindex(vel, p.rz) + e.nbc;
            if (p.flag_homo)
                apply_homo(e, sx, sz);
            propagate(e, out.nt, out.dt, p.type_compute_laplace, false, 1, sx, sz, rz, true, p.progress_interval, label, &src, &rec, nullptr, nullptr, nullptr);
            std::memcpy(&out.d[(size_t)is * out.nr * out.nt], rec.data(), rec.size() * sizeof(float));
        }
        return out;
    }

    std::vector<float> rtm_image(const Grid2D &vel, const Data3D &data, const ImageParams &p, const Grid2D *qmodel)
    {
        ExtModel ebase = make_ext(vel, qmodel, p);
        int shot0 = std::max(0, p.shot_begin);
        int shot1 = (p.shot_end < 0) ? data.ns - 1 : std::min(data.ns - 1, p.shot_end);
        if (shot0 > shot1)
            throw std::runtime_error("empty shot range for imaging");
        int shot_threads = std::max(1, p.max_parallel_shots);
#ifdef _OPENMP
        shot_threads = std::min(shot_threads, omp_get_max_threads());
#endif
        int snapshot_interval = std::max(1, p.snapshot_interval);
        int rec_snapshot_offset = 0;
        int src_snapshot_offset = (data.nt - 1) % snapshot_interval;
        int stored_snaps = count_stored_snapshots(data.nt, snapshot_interval, rec_snapshot_offset);
        size_t snap_bytes = (size_t)stored_snaps * ebase.nx * ebase.nz * sizeof(float);
        int amp_sign = p.compensate ? -1 : 1;
        INFO(("RTM shot parallelism=%d, snapshot_interval=%d, stored snapshots per wavefield=%d/%d, source snapshots per active shot=%.3f GiB (%s)",
              shot_threads, snapshot_interval, stored_snaps, data.nt, snap_bytes / (1024.0 * 1024.0 * 1024.0), p.in_memory_snapshots ? "memory" : "disk temporary file"));
        ensure_directory("temp_shot_img", "per-shot images");
        std::vector<int> shot_imaged(data.ns, 0);
#pragma omp parallel for schedule(dynamic) num_threads(shot_threads)
        for (int is = shot0; is <= shot1; ++is)
        {
            ExtModel e = ebase;
            float sx = data.s0 + is * data.ds;
            int six0 = xindex_unclamped(vel, sx);
            int siz0 = zindex_unclamped(vel, p.sz);
            if (six0 < 0 || six0 >= vel.nx || siz0 < 0 || siz0 >= vel.nz)
            {
#pragma omp critical(sewave_shot_skip_log)
                {
                    INFO(("image shot %d skipped: source position outside velocity model (shot_x=%g shot_z=%g, ix=%d iz=%d, valid ix=0-%d iz=0-%d)",
                          is, sx, p.sz, six0, siz0, vel.nx - 1, vel.nz - 1));
                }
                continue;
            }
            std::vector<float> local((size_t)vel.nx * vel.nz, 0.0f);
            std::vector<float> src(data.nt), shot((size_t)(e.nx - 2 * e.nbc) * data.nt, 0.0f), snaps;
            if (p.in_memory_snapshots)
                snaps.resize((size_t)stored_snaps * e.nx * e.nz);
            std::string snap_path, rec_snap_path;
            FILE *snap_tmp = p.in_memory_snapshots ? nullptr : open_snapshot_temp(is, snap_path);
            FILE *rec_snap_tmp = open_snapshot_temp(is, rec_snap_path);
            char src_label[80], rec_label[80];
            std::snprintf(src_label, sizeof(src_label), "image source shot %d", is);
            std::snprintf(rec_label, sizeof(rec_label), "image receiver shot %d", is);
            for (int it = 0; it < data.nt; ++it)
                src[it] = ricker(it * data.dt, p.fdom);
            int six = six0 + e.nbc, siz = siz0 + e.nbc;
            if (p.flag_homo)
                apply_homo(e, six, siz);
            propagate(e, data.nt, data.dt, p.type_compute_laplace, true, amp_sign, six, siz, siz, false, p.progress_interval, src_label, &src, nullptr, p.in_memory_snapshots ? &snaps : nullptr, nullptr, nullptr, snap_tmp, nullptr,
                      snapshot_interval, src_snapshot_offset, snapshot_interval, src_snapshot_offset,
                      &vel, &p.debug_snapshot_steps, &p.debug_wavefield_prefix, "src", is);
            int skipped_receivers = 0;
            for (int ir = 0; ir < data.nr; ++ir)
            {
                float rx = p.cmp ? sx + data.r0 + ir * data.dr : data.r0 + ir * data.dr;
                int ix = xindex_unclamped(vel, rx);
                if (ix < 0 || ix >= vel.nx)
                {
                    ++skipped_receivers;
                    continue;
                }
                for (int it = 0; it < data.nt; ++it)
                    shot[(size_t)ix * data.nt + it] += data.d[((size_t)is * data.nr + ir) * data.nt + it];
            }
            if (skipped_receivers > 0)
                INFO(("image shot %d: skipped %d receivers outside velocity model", is, skipped_receivers));
            if (!p.debug_wavefield_prefix.empty())
                report_gather_energy(shot, e.nx - 2 * e.nbc, data.nt, is, p.debug_wavefield_prefix, vel.x0, vel.dx, data.dt);
            propagate(e, data.nt, data.dt, p.type_compute_laplace, true, amp_sign, 0, zindex(vel, p.rz) + e.nbc, zindex(vel, p.rz) + e.nbc, false, p.progress_interval, rec_label, &shot, nullptr, nullptr, nullptr, nullptr, rec_snap_tmp, nullptr,
                      snapshot_interval, rec_snapshot_offset, snapshot_interval, src_snapshot_offset,
                      &vel, &p.debug_snapshot_steps, &p.debug_wavefield_prefix, "rec", is);
            correlate_snapshots(e, data.nt, snapshot_interval, src_snapshot_offset, rec_snapshot_offset, p.in_memory_snapshots ? &snaps : nullptr, snap_tmp, rec_snap_tmp, local);
            if (rec_snap_tmp)
            {
                std::fclose(rec_snap_tmp);
                if (!rec_snap_path.empty())
                    ::unlink(rec_snap_path.c_str());
            }
            if (snap_tmp)
            {
                std::fclose(snap_tmp);
                if (!snap_path.empty())
                    ::unlink(snap_path.c_str());
            }
            std::string shot_img_path = shot_image_path(is);
#pragma omp critical(sewave_shot_image_write)
            {
                write_shot_image_rsf(shot_img_path, vel, local, is, sx, p.sz);
                shot_imaged[is] = 1;
                INFO(("Wrote per-shot image %s with shot_x=%g shot_z=%g", shot_img_path.c_str(), sx, p.sz));
            }
        }
        std::vector<float> image((size_t)vel.nx * vel.nz, 0.0f);
        int stacked_shots = 0;
        for (int is = shot0; is <= shot1; ++is)
        {
            if (!shot_imaged[is])
                continue;
            std::string shot_img_path = shot_image_path(is);
            Grid2D shot_img = read_rsf2d(shot_img_path);
            check_same_grid(vel, shot_img, shot_img_path.c_str());
            for (size_t i = 0; i < image.size(); ++i)
                image[i] += shot_img.v[i];
            ++stacked_shots;
        }
        if (stacked_shots == 0)
            INFO(("image shot stack: no shots were inside the velocity model; returning zero image"));
        else if (stacked_shots < shot1 - shot0 + 1)
            INFO(("image shot stack: stacked %d/%d shots; skipped %d shots outside the velocity model",
                  stacked_shots, shot1 - shot0 + 1, shot1 - shot0 + 1 - stacked_shots));
        return image;
    }
} // namespace sewave
