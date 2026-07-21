#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

int main(int argc, char **argv)
{
    se_par_init(argc, argv);
    if (!se_have_par("out")) ERROR(("Need out= output SEP velocity model"));
    const char *out_f = se_get_par_str("out");
    const int nz = se_have_par("nz") ? se_get_par_int("nz") : 81;
    const int nx = se_have_par("nx") ? se_get_par_int("nx") : 96;
    const float dz = se_have_par("dz") ? se_get_par_float("dz") : 10.0f;
    const float dx = se_have_par("dx") ? se_get_par_float("dx") : 10.0f;
    const float oz = se_have_par("oz") ? se_get_par_float("oz") : 0.0f;
    const float ox = se_have_par("ox") ? se_get_par_float("ox") : 0.0f;
    if (nz <= 1 || nx <= 1 || dz <= 0.0f || dx <= 0.0f) ERROR(("Invalid grid"));

    std::vector<float> vel((size_t)nz * (size_t)nx);
    const float xmax = ox + (nx - 1) * dx;
    const float xc = ox + 0.5f * (nx - 1) * dx;
    for (int ix = 0; ix < nx; ++ix)
    {
        const float x = ox + ix * dx;
        for (int iz = 0; iz < nz; ++iz)
        {
            const float z = oz + iz * dz;
            const float sag = 180.0f * std::exp(-std::pow((x - xc) / (0.22f * (xmax - ox + dx)), 2.0f));
            const float h1 = oz + 0.28f * (nz - 1) * dz + 0.10f * (x - ox) + sag;
            const float h2 = oz + 0.58f * (nz - 1) * dz - 0.07f * (x - ox) - 0.45f * sag;
            float v = 1800.0f + 0.35f * z;
            if (z > h1) v += 350.0f;
            if (z > h2) v += 500.0f;
            vel[(size_t)ix * (size_t)nz + iz] = v;
        }
    }

    const std::string data_f = std::string(out_f) + "@";
    FILE *fp = std::fopen(data_f.c_str(), "wb");
    if (!fp) ERROR(("Cannot open %s", data_f.c_str()));
    if (std::fwrite(vel.data(), sizeof(float), vel.size(), fp) != vel.size()) ERROR(("Failed writing velocity data"));
    std::fclose(fp);

    std::ofstream h(out_f);
    h << "n1=" << nz << " o1=" << oz << " d1=" << dz << "\n";
    h << "n2=" << nx << " o2=" << ox << " d2=" << dx << "\n";
    h << "esize=4 data_format=native_float in=\"" << data_f << "\"\n";
    h << "label1=Depth label2=Distance unit1=m unit2=m\n";
    h.close();
    INFO(("Wrote layered/sag/dipping velocity model %s (%dx%d)", out_f, nz, nx));
    se_par_destroy();
    return 0;
}
