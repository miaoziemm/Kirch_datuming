#include "../include/tpow.h"


void tpow_gain_trace(
    float *trace,
    int nt,
    float tpow,
    float dt,
    float t0
)
{
    if (trace == nullptr || nt <= 0) {
        return;
    }

    if (tpow == 0.0f) {
        return;
    }

    for (int it = 0; it < nt; ++it) {
        float t = t0 + static_cast<float>(it + 1) * dt;
        trace[it] *= std::pow(t, tpow);
    }
}


void tpow_gain(
    float *data,
    int nt,
    int nx,
    int nblk,
    float tpow,
    float xpow,
    float dt,
    float t0,
    float dx,
    float x0
)
{
    if (data == nullptr || nt <= 0 || nx <= 0 || nblk <= 0) {
        return;
    }

    std::vector<float> tgain;
    std::vector<float> xgain;

    if (tpow != 0.0f) {
        tgain.resize(static_cast<std::size_t>(nt));

        for (int it = 0; it < nt; ++it) {
            float t = t0 + static_cast<float>(it + 1) * dt;
            tgain[static_cast<std::size_t>(it)] = std::pow(t, tpow);
        }
    }

    if (xpow != 0.0f) {
        xgain.resize(static_cast<std::size_t>(nx));

        for (int ix = 0; ix < nx; ++ix) {
            float x = x0 + static_cast<float>(ix + 1) * dx;
            xgain[static_cast<std::size_t>(ix)] = std::pow(x, xpow);
        }
    }

    for (int iblk = 0; iblk < nblk; ++iblk) {
        for (int ix = 0; ix < nx; ++ix) {

            const std::size_t trace_offset =
                (static_cast<std::size_t>(iblk) * static_cast<std::size_t>(nx)
                 + static_cast<std::size_t>(ix)) * static_cast<std::size_t>(nt);

            float xcoef = 1.0f;
            if (xpow != 0.0f) {
                xcoef = xgain[static_cast<std::size_t>(ix)];
            }

            for (int it = 0; it < nt; ++it) {
                float gain = xcoef;

                if (tpow != 0.0f) {
                    gain *= tgain[static_cast<std::size_t>(it)];
                }

                data[trace_offset + static_cast<std::size_t>(it)] *= gain;
            }
        }
    }
}