#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Boundary { int iz; float z; };

static std::vector<float> parse_float_list(const char *text)
{
    std::vector<float> vals;
    if (text == NULL) return vals;
    std::string s(text);
    std::replace(s.begin(), s.end(), ',', ' ');
    std::stringstream ss(s);
    float v;
    while (ss >> v) vals.push_back(v);
    return vals;
}

static int clamp_index(int i, int n)
{
    return std::max(0, std::min(n - 1, i));
}

static float velocity_at(const std::vector<float> &vel, int nz, int nx, int iz, int ix)
{
    iz = clamp_index(iz, nz);
    ix = clamp_index(ix, nx);
    const float v = vel[(size_t)ix * (size_t)nz + (size_t)iz];
    if (!(v > 0.0f) || !std::isfinite(v))
        ERROR(("Invalid velocity sample at iz=%d ix=%d", iz, ix));
    return v;
}

static void build_tau_fmm(const std::vector<float> &vel, int nz, int nx,
                          float dz, float dx,
                          int iz_top, int iz_out, std::vector<float> &tau)
{
    if (iz_out < iz_top)
        ERROR(("FMM tau requires iz_out >= iz_top"));
    const int nzl = iz_out - iz_top + 1;
    tau.assign((size_t)nx * (size_t)nx, 0.0f);

    std::vector<float> subvel((size_t)nx * (size_t)nzl);
    for (int iz = 0; iz < nzl; ++iz)
        for (int ix = 0; ix < nx; ++ix)
            subvel[(size_t)ix + (size_t)iz * (size_t)nx] = velocity_at(vel, nz, nx, iz_top + iz, ix);

    std::vector<float> tt((size_t)nx * (size_t)nzl);
    for (int ix_src = 0; ix_src < nx; ++ix_src)
    {
        efmm_t ctx;
        if (efmm_init(&ctx, nx, nzl, dx, dz, ix_src * dx, 0.0f) != 0)
            ERROR(("efmm_init failed for source ix=%d", ix_src));
        if (efmm_set_vel(&ctx, subvel.data()) != 0 || efmm_solver(&ctx) != 0)
        {
            efmm_free(&ctx);
            ERROR(("efmm_solver failed for source ix=%d", ix_src));
        }
        std::copy(ctx.tt, ctx.tt + (size_t)nx * (size_t)nzl, tt.begin());
        efmm_free(&ctx);

        for (int ix_out = 0; ix_out < nx; ++ix_out)
            tau[(size_t)ix_out * (size_t)nx + (size_t)ix_src] = tt[(size_t)ix_out + (size_t)(nzl - 1) * (size_t)nx];
    }
}

static void write_tau_file(const std::string &hdr, const std::vector<float> &tau,
                           int nx, float ox, float dx, int layer, int iz_top, int iz_bot,
                           float z_top, float z_bot)
{
    const std::string data = hdr + "@";
    FILE *fp = std::fopen(data.c_str(), "wb");
    if (!fp) ERROR(("Cannot open tau data file %s", data.c_str()));
    if (std::fwrite(tau.data(), sizeof(float), tau.size(), fp) != tau.size())
        ERROR(("Failed to write tau data %s", data.c_str()));
    std::fclose(fp);
    std::ofstream h(hdr.c_str());
    h << "n1=" << nx << " o1=" << ox << " d1=" << dx << "\n";
    h << "n2=" << nx << " o2=" << ox << " d2=" << dx << "\n";
    h << "esize=4 data_format=native_float in=\"" << data << "\"\n";
    h << "layer=" << layer << " iz_top=" << iz_top << " iz_bot=" << iz_bot
      << " z_top=" << z_top << " z_bot=" << z_bot << "\n";
}

static std::string freq_tag(float freq)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%g", freq);
    std::string s(buf);
    for (char &c : s)
        if (c == '.' || c == '-' || c == '+') c = '_';
    return s;
}

static void write_complex_boundary_file(const std::string &hdr,
                                        const std::vector<std::complex<float>> &field,
                                        int nx, float ox, float dx, float z,
                                        float freq, int layer)
{
    const std::string real_data = hdr + ".real";
    const std::string imag_data = hdr + ".imag";
    FILE *fr = std::fopen(real_data.c_str(), "wb");
    FILE *fi = std::fopen(imag_data.c_str(), "wb");
    if (!fr || !fi) ERROR(("Cannot open split boundary data files for %s", hdr.c_str()));
    std::vector<float> real_row((size_t)nx), imag_row((size_t)nx);
    for (int ix = 0; ix < nx; ++ix)
    {
        real_row[(size_t)ix] = field[ix].real();
        imag_row[(size_t)ix] = field[ix].imag();
    }
    if (std::fwrite(real_row.data(), sizeof(float), real_row.size(), fr) != real_row.size() ||
        std::fwrite(imag_row.data(), sizeof(float), imag_row.size(), fi) != imag_row.size())
        ERROR(("Failed to write split boundary data %s", hdr.c_str()));
    std::fclose(fr);
    std::fclose(fi);
    std::ofstream h(hdr.c_str());
    h << "n1=" << nx << " o1=" << ox << " d1=" << dx << "\n";
    h << "esize=4 data_format=native_float complex_storage=split_float in_real=\"" << real_data << "\" in_imag=\"" << imag_data << "\"\n";
    h << "freq=" << freq << " layer=" << layer << " z=" << z << "\n";
}

static void write_complex_block_file(const std::string &hdr,
                                     const std::vector<std::vector<std::complex<float>>> &block,
                                     int nx, float ox, float dx, float z0, float dz,
                                     float freq, int layer)
{
    const int nzl = (int)block.size();
    const std::string real_data = hdr + ".real";
    const std::string imag_data = hdr + ".imag";
    FILE *fr = std::fopen(real_data.c_str(), "wb");
    FILE *fi = std::fopen(imag_data.c_str(), "wb");
    if (!fr || !fi) ERROR(("Cannot open split block data files for %s", hdr.c_str()));
    std::vector<float> real_col((size_t)nzl), imag_col((size_t)nzl);
    for (int ix = 0; ix < nx; ++ix)
    {
        for (int iz = 0; iz < nzl; ++iz)
        {
            real_col[(size_t)iz] = block[(size_t)iz][ix].real();
            imag_col[(size_t)iz] = block[(size_t)iz][ix].imag();
        }
        if (std::fwrite(real_col.data(), sizeof(float), real_col.size(), fr) != real_col.size() ||
            std::fwrite(imag_col.data(), sizeof(float), imag_col.size(), fi) != imag_col.size())
            ERROR(("Failed to write split block data %s", hdr.c_str()));
    }
    std::fclose(fr);
    std::fclose(fi);
    std::ofstream h(hdr.c_str());
    h << "n1=" << nzl << " o1=" << z0 << " d1=" << dz << "\n";
    h << "n2=" << nx << " o2=" << ox << " d2=" << dx << "\n";
    h << "esize=4 data_format=native_float complex_storage=split_float in_real=\"" << real_data << "\" in_imag=\"" << imag_data << "\"\n";
    h << "freq=" << freq << " layer=" << layer << "\n";
}

static float taper_weight(int left, int i, int right, int tap)
{
    if (tap <= 0) return 1.0f;
    float coef = 1.0f;
    coef *= (i - left >= tap) ? 1.0f : (float)(i - left) / (float)tap;
    coef *= (right - i >= tap) ? 1.0f : (float)(right - i) / (float)tap;
    return coef;
}

static void build_kirdat_filter_F(float tau, float dt, int nsam, std::vector<float> &F)
{
    F.assign((size_t)std::max(0, nsam - 1), 0.0f);
    if (tau <= 0.0f || nsam <= 1) return;
    float prev2 = 0.0f;
    float prev1 = std::sqrt(((tau + dt) / tau) * ((tau + dt) / tau) - 1.0f);
    F[0] = prev1 - prev2;
    for (int k = 1; k < nsam - 1; ++k)
    {
        const float t = (tau + (float)(k + 1) * dt) / tau;
        const float curr = std::sqrt(t * t - 1.0f);
        F[(size_t)k] = curr - 2.0f * prev1 + prev2;
        prev2 = prev1;
        prev1 = curr;
    }
}

static std::complex<float> kirdat_frequency_filter(float tau, float omega, float dt, int nsam)
{
    if (tau <= 0.0f || !std::isfinite(tau)) return std::complex<float>(0.0f, 0.0f);
    std::vector<float> F;
    build_kirdat_filter_F(tau, dt, nsam, F);
    const float q = tau / dt;
    const int m = (int)std::ceil(q);
    const float delta = (float)m - q;
    std::complex<float> sum(0.0f, 0.0f);
    for (int k = 0; k <= nsam - 2; ++k)
    {
        const float phase = -omega * dt * (float)k;
        sum += (F[(size_t)k] / dt) * std::complex<float>(std::cos(phase), std::sin(phase));
    }
    const std::complex<float> e_m(std::cos(-omega * dt * (float)m), std::sin(-omega * dt * (float)m));
    const std::complex<float> e_m1(std::cos(-omega * dt * (float)(m - 1)), std::sin(-omega * dt * (float)(m - 1)));
    return ((1.0f - delta) * e_m + delta * e_m1) * sum;
}

static std::complex<float> kirdat_matrix_value(float tau, int ix_out, int ix_in, int nx,
                                                float dx, float datum, float omega,
                                                float dt, int nsam, int aperture, int taper,
                                                bool bf_amplitude)
{
    if (tau <= 0.0f || !std::isfinite(tau)) return std::complex<float>(0.0f, 0.0f);
    const int left = std::max(0, ix_out - aperture);
    const int right = std::min(nx - 1, ix_out + aperture);
    if (ix_in < left || ix_in > right) return std::complex<float>(0.0f, 0.0f);
    const float coef = taper_weight(left, ix_in, right, taper);
    const float lateral = (float)(ix_in - ix_out) * dx;
    const float dist = datum * datum + lateral * lateral;
    if (dist <= 0.0f) return std::complex<float>(0.0f, 0.0f);
    const float w = coef / (float)M_PI * dx * datum * tau / dist;
    std::complex<float> amp = w * kirdat_frequency_filter(tau, omega, dt, nsam);
    if (bf_amplitude)
        amp *= std::complex<float>(std::cos(omega * tau), std::sin(omega * tau));
    return amp;
}

static void apply_direct(const std::vector<float> &tau, int nx, float dx, float datum,
                         float omega, float dt, int nsam, int aperture, int taper,
                         const std::vector<std::complex<float>> &uin,
                         std::vector<std::complex<float>> &uout)
{
    std::fill(uout.begin(), uout.end(), std::complex<float>(0.0f, 0.0f));
    for (int ix_out = 0; ix_out < nx; ++ix_out)
    {
        std::complex<float> sum(0.0f, 0.0f);
        for (int ix_in = 0; ix_in < nx; ++ix_in)
        {
            const float t = tau[(size_t)ix_out * (size_t)nx + (size_t)ix_in];
            sum += kirdat_matrix_value(t, ix_out, ix_in, nx, dx, datum, omega, dt, nsam, aperture, taper, false) * uin[ix_in];
        }
        uout[ix_out] = sum;
    }
}

static float ricker_spectrum(float freq, float fpeak)
{
    if (fpeak <= 0.0f) return 1.0f;
    const float a = freq / fpeak;
    return a * a * std::exp(-a * a);
}

static void source_boundary(int nx, float ox, float dx, float z, float srcx, float srcz,
                            float freq, float vref, float fpeak,
                            std::vector<std::complex<float>> &u)
{
    const float omega = 2.0f * (float)M_PI * freq;
    const float amp = ricker_spectrum(freq, fpeak);
    for (int ix = 0; ix < nx; ++ix)
    {
        const float x = ox + ix * dx;
        const float r = std::max(std::hypot(x - srcx, z - srcz), 0.25f * std::fabs(dx));
        const float phase = omega * r / vref;
        u[ix] = amp * std::complex<float>(std::cos(phase), -std::sin(phase)) / std::sqrt(r);
    }
}

static std::vector<std::complex<float>> data_boundary_at_freq(const std::vector<float> &data, int nt, int nx, float ot, float dt, float freq)
{
    std::vector<std::complex<float>> spec(nx, std::complex<float>(0.0f, 0.0f));
    const float omega = 2.0f * (float)M_PI * freq;
    for (int ix = 0; ix < nx; ++ix)
        for (int it = 0; it < nt; ++it)
        {
            const float t = ot + it * dt;
            const std::complex<float> e(std::cos(omega * t), -std::sin(omega * t));
            spec[ix] += data[(size_t)ix * (size_t)nt + it] * e * dt;
        }
    return spec;
}

int main(int argc, char **argv)
{
    se_par_init(argc, argv);
    if (!se_have_par("model")) ERROR(("Need model= velocity SEP file"));
    if (!se_have_par("out")) ERROR(("Need out= output wavefield SEP file"));
    if (!se_have_par("tau_prefix")) ERROR(("Need tau_prefix= prefix for tau check files"));
    if (!se_have_par("nlayer")) ERROR(("Need nlayer= number of uniform layers"));
    if (!se_have_par("freqs")) ERROR(("Need freqs= comma-separated frequencies"));

    const char *model_f = se_get_par_str("model");
    const char *out_f = se_get_par_str("out");
    const std::string tau_prefix = se_get_par_str("tau_prefix");
    const std::string boundary_prefix = se_have_par("boundary_prefix") ? se_get_par_str("boundary_prefix") : std::string(out_f) + "_boundary";
    const std::string block_prefix = se_have_par("block_prefix") ? se_get_par_str("block_prefix") : std::string(out_f) + "_block";
    const std::string out_real_f = se_have_par("out_real") ? se_get_par_str("out_real") : std::string(out_f) + ".real";
    const std::string out_imag_f = se_have_par("out_imag") ? se_get_par_str("out_imag") : std::string(out_f) + ".imag";
    const int nlayer = se_get_par_int("nlayer");
    const std::string mode = se_have_par("mode") ? se_get_par_str("mode") : "source";
    const int use_bf = se_have_par("use_bf") ? se_get_par_int("use_bf") : 1;
    const int bf_p = se_have_par("bf_p") ? se_get_par_int("bf_p") : 12;
    const int bf_leaf = se_have_par("bf_n_leaf") ? se_get_par_int("bf_n_leaf") : 16;
    const int bf_panel_levels = se_have_par("bf_panel_levels") ? se_get_par_int("bf_panel_levels") : 1;
    const float bf_amp_eps = se_have_par("bf_amp_eps") ? se_get_par_float("bf_amp_eps") : 1.0e-20f;
    const float bf_phase_tol = se_have_par("bf_phase_tol") ? se_get_par_float("bf_phase_tol") : 2.0f;
    const float kirdat_dt = se_have_par("dt") ? se_get_par_float("dt") : 0.001f;
    const float kirdat_length = se_have_par("length") ? se_get_par_float("length") : 0.025f;
    const int taper = se_have_par("taper") ? se_get_par_int("taper") : 10;
    int nfreq = 0;
    float *freq_array = se_get_pararray_float("freqs", &nfreq);
    std::vector<float> freqs;
    if (freq_array != NULL && nfreq > 0)
    {
        freqs.assign(freq_array, freq_array + nfreq);
        free1float(freq_array);
    }
    else
    {
        freqs = parse_float_list(se_get_par_str("freqs"));
    }
    if (nlayer <= 0) ERROR(("nlayer must be positive"));
    if (freqs.empty()) ERROR(("No frequencies parsed from freqs="));
    if (kirdat_dt <= 0.0f || kirdat_length <= 0.0f) ERROR(("dt= and length= must be positive for the kirdat frequency-domain filter"));
    const int nsam = (int)(kirdat_length / kirdat_dt) + 2;

    sep_t *model = sep_open(model_f, SEP_READ, 0);
    sep_t *out = sep_open(out_f, SEP_WRITE, 0);
    if (model->headers->ndim < 2) ERROR(("model must be 2-D with n1=z,n2=x"));
    const int nz = model->headers->n[0];
    const int nx = model->headers->n[1];
    const float dz = (float)model->headers->d[0];
    const float dx = (float)model->headers->d[1];
    const float oz = (float)model->headers->o[0];
    const float ox = (float)model->headers->o[1];
    if (nlayer >= nz) ERROR(("nlayer=%d must be smaller than nz=%d", nlayer, nz));
    const int aperture = se_have_par("aperture") ? se_get_par_int("aperture") : nx;

    std::vector<float> vel((size_t)nz * (size_t)nx);
    se_fsio_read_float(model->data->io, vel.data(), vel.size());

    std::vector<Boundary> bnd;
    bnd.push_back({0, oz});
    for (int il = 1; il <= nlayer; ++il)
    {
        int iz = (int)std::llround((double)il * (double)(nz - 1) / (double)nlayer);
        if (iz <= bnd.back().iz) iz = bnd.back().iz + 1;
        if (iz > nz - 1) iz = nz - 1;
        bnd.push_back({iz, oz + iz * dz});
    }
    bnd.back().iz = nz - 1;
    bnd.back().z = oz + (nz - 1) * dz;

    std::vector<std::vector<float>> layer_tau((size_t)nlayer);
    for (int il = 0; il < nlayer; ++il)
    {
        build_tau_fmm(vel, nz, nx, dz, dx, bnd[il].iz, bnd[il + 1].iz, layer_tau[il]);
        char name[1024];
        std::snprintf(name, sizeof(name), "%s_layer%03d.H", tau_prefix.c_str(), il + 1);
        write_tau_file(name, layer_tau[il], nx, ox, dx, il + 1, bnd[il].iz, bnd[il + 1].iz, bnd[il].z, bnd[il + 1].z);
        INFO(("Wrote tau for layer %d: iz %d -> %d", il + 1, bnd[il].iz, bnd[il + 1].iz));
    }

    sep_t *shot = NULL;
    int nt = 0;
    float ot = 0.0f, dt = 0.0f;
    std::vector<float> shot_data;
    if (mode == "data")
    {
        if (!se_have_par("data")) ERROR(("mode=data requires data= single-shot gather"));
        shot = sep_open(se_get_par_str("data"), SEP_READ, 0);
        if (shot->headers->ndim < 2 || shot->headers->n[1] != nx)
            ERROR(("data must have n1=time and n2=model nx=%d", nx));
        nt = shot->headers->n[0];
        ot = (float)shot->headers->o[0];
        dt = (float)shot->headers->d[0];
        shot_data.resize((size_t)nt * (size_t)nx);
        se_fsio_read_float(shot->data->io, shot_data.data(), shot_data.size());
    }

    out->headers->ndim = 4;
    out->headers->n[0] = 2;
    out->headers->n[1] = nx;
    out->headers->n[2] = nz;
    out->headers->n[3] = (int)freqs.size();
    out->headers->d[0] = 1.0;
    out->headers->d[1] = dx;
    out->headers->d[2] = dz;
    out->headers->d[3] = 1.0;
    out->headers->o[0] = 0.0;
    out->headers->o[1] = ox;
    out->headers->o[2] = oz;
    out->headers->o[3] = 0.0;

    const float srcx = se_have_par("srcx") ? se_get_par_float("srcx") : ox + 0.5f * (nx - 1) * dx;
    const float srcz = se_have_par("srcz") ? se_get_par_float("srcz") : oz;
    const float vref = se_have_par("vref") ? se_get_par_float("vref") : velocity_at(vel, nz, nx, 0, nx / 2);
    const float fpeak = se_have_par("fpeak") ? se_get_par_float("fpeak") : 25.0f;

    FILE *final_real = std::fopen(out_real_f.c_str(), "wb");
    FILE *final_imag = std::fopen(out_imag_f.c_str(), "wb");
    if (!final_real || !final_imag) ERROR(("Cannot open split final wavefield files %s / %s", out_real_f.c_str(), out_imag_f.c_str()));

    std::vector<float> row((size_t)2 * (size_t)nx);
    std::vector<float> final_real_vol((size_t)nx * (size_t)nz), final_imag_vol((size_t)nx * (size_t)nz);
    auto store_final_split = [&](const std::vector<std::complex<float>> &field, int iz_global) {
        for (int ix = 0; ix < nx; ++ix)
        {
            final_real_vol[(size_t)ix * (size_t)nz + (size_t)iz_global] = field[ix].real();
            final_imag_vol[(size_t)ix * (size_t)nz + (size_t)iz_global] = field[ix].imag();
        }
    };
    auto flush_final_split = [&]() {
        if (std::fwrite(final_real_vol.data(), sizeof(float), final_real_vol.size(), final_real) != final_real_vol.size() ||
            std::fwrite(final_imag_vol.data(), sizeof(float), final_imag_vol.size(), final_imag) != final_imag_vol.size())
            ERROR(("Failed to write split final wavefield volume"));
    };
    std::vector<std::complex<float>> u_top(nx), u_bot(nx), u_z(nx);
    for (size_t jf = 0; jf < freqs.size(); ++jf)
    {
        const float freq = freqs[jf];
        const float omega = 2.0f * (float)M_PI * freq;
        std::fill(final_real_vol.begin(), final_real_vol.end(), 0.0f);
        std::fill(final_imag_vol.begin(), final_imag_vol.end(), 0.0f);
        if (mode == "data")
            u_top = data_boundary_at_freq(shot_data, nt, nx, ot, dt, freq);
        else if (mode == "source")
            source_boundary(nx, ox, dx, oz, srcx, srcz, freq, vref, fpeak, u_top);
        else
            ERROR(("Unknown mode=%s; use mode=source or mode=data", mode.c_str()));

        {
            char bname[1024];
            std::snprintf(bname, sizeof(bname), "%s_f%s_layer%03d.H", boundary_prefix.c_str(), freq_tag(freq).c_str(), 0);
            write_complex_boundary_file(bname, u_top, nx, ox, dx, bnd.front().z, freq, 0);
        }

        for (int il = 0; il < nlayer; ++il)
        {
            const int nzl_block = bnd[il + 1].iz - bnd[il].iz + 1;
            std::vector<std::vector<std::complex<float>>> block_field((size_t)nzl_block, std::vector<std::complex<float>>(nx));
            block_field[0] = u_top;

            if (il == 0)
            {
                for (int ix = 0; ix < nx; ++ix) { row[2 * ix] = u_top[ix].real(); row[2 * ix + 1] = u_top[ix].imag(); }
                se_fsio_write_float(out->data->io, row.data(), row.size());
                store_final_split(u_top, bnd[il].iz);
            }

            if (use_bf)
            {
                float **tau_ptr = alloc2float(nx, nx);
                fftwf_complex *amp = (fftwf_complex *)fftwf_malloc((size_t)nx * (size_t)nx * sizeof(fftwf_complex));
                fftwf_complex *bin = (fftwf_complex *)fftwf_malloc((size_t)nx * sizeof(fftwf_complex));
                fftwf_complex *bout = (fftwf_complex *)fftwf_malloc((size_t)nx * sizeof(fftwf_complex));
                if (!tau_ptr || !amp || !bin || !bout) ERROR(("Out of memory for butterfly work arrays"));
                for (int ix_out = 0; ix_out < nx; ++ix_out)
                    for (int ix_in = 0; ix_in < nx; ++ix_in)
                    {
                        const float t = layer_tau[il][(size_t)ix_out * nx + ix_in];
                        tau_ptr[ix_out][ix_in] = t;
                        const size_t k = (size_t)ix_out * nx + ix_in;
                        const std::complex<float> a = kirdat_matrix_value(t, ix_out, ix_in, nx, dx,
                                                                          std::fabs(bnd[il + 1].z - bnd[il].z),
                                                                          omega, kirdat_dt, nsam, aperture, taper, true);
                        amp[k][0] = a.real();
                        amp[k][1] = a.imag();
                    }
                for (int ix = 0; ix < nx; ++ix) { bin[ix][0] = u_top[ix].real(); bin[ix][1] = u_top[ix].imag(); }
                BFStrictSegmentedFactor *factor = bf1d_strict_segmented_create_phase_amp(nx, tau_ptr, amp, omega, bf_p, bf_leaf, bf_panel_levels, bf_amp_eps, bf_phase_tol);
                bf1d_strict_segmented_apply(factor, bin, bout);
                bf1d_strict_segmented_destroy(factor);
                for (int ix = 0; ix < nx; ++ix) u_bot[ix] = std::complex<float>(bout[ix][0], bout[ix][1]);
                free2float(tau_ptr); fftwf_free(amp); fftwf_free(bin); fftwf_free(bout);
            }
            else
            {
                apply_direct(layer_tau[il], nx, dx, std::fabs(bnd[il + 1].z - bnd[il].z), omega, kirdat_dt, nsam, aperture, taper, u_top, u_bot);
            }

            for (int iz = bnd[il].iz + 1; iz <= bnd[il + 1].iz; ++iz)
            {
                if (iz == bnd[il + 1].iz)
                    u_z = u_bot;
                else
                {
                    std::vector<float> tau_z;
                    build_tau_fmm(vel, nz, nx, dz, dx, bnd[il].iz, iz, tau_z);
                    apply_direct(tau_z, nx, dx, std::fabs((oz + iz * dz) - bnd[il].z), omega, kirdat_dt, nsam, aperture, taper, u_top, u_z);
                }
                block_field[(size_t)(iz - bnd[il].iz)] = u_z;
                for (int ix = 0; ix < nx; ++ix)
                {
                    row[(size_t)2 * ix] = u_z[ix].real();
                    row[(size_t)2 * ix + 1] = u_z[ix].imag();
                }
                se_fsio_write_float(out->data->io, row.data(), row.size());
                store_final_split(u_z, iz);
            }
            {
                char bname[1024];
                std::snprintf(bname, sizeof(bname), "%s_f%s_layer%03d.H", boundary_prefix.c_str(), freq_tag(freq).c_str(), il + 1);
                write_complex_boundary_file(bname, u_bot, nx, ox, dx, bnd[il + 1].z, freq, il + 1);
                char lname[1024];
                std::snprintf(lname, sizeof(lname), "%s_f%s_layer%03d.H", block_prefix.c_str(), freq_tag(freq).c_str(), il + 1);
                write_complex_block_file(lname, block_field, nx, ox, dx, bnd[il].z, dz, freq, il + 1);
                INFO(("Finished block wavefield for frequency %g Hz, layer %d/%d: iz %d -> %d",
                      freq, il + 1, nlayer, bnd[il].iz, bnd[il + 1].iz));
            }
            u_top.swap(u_bot);
        }
        flush_final_split();
        INFO(("Finished frequency %g Hz", freq));
    }

    std::fclose(final_real);
    std::fclose(final_imag);

    if (shot) sep_close(shot);
    sep_close(model);
    sep_close(out);
    se_par_destroy();
    return 0;
}
