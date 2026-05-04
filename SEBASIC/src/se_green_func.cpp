#include "se_green_func.h"
#include "se_math_physics.h"
#include <cmath>
#include <complex>
#include <iostream>

#ifndef PI
#define PI 3.14159265358979323846
#endif

void se_green_func_2d_homogeneous(se_green_2d_params_t *params, float *real, float *imag)
{
    if (!params || !real || !imag) return;

    float v = params->v;
    float f = params->f;
    float dx = params->rec_x - params->src_x;
    float dz = params->rec_z - params->src_z;
    
    float r = std::sqrt(dx * dx + dz * dz);
    
    /* Avoid singularity at r=0 */
    if (r < 1e-6f) {
        *real = 0.0f; /* Or handle singularity appropriately */
        *imag = 0.0f;
        return;
    }

    float omega = 2.0f * (float)PI * f;
    float k = omega / v;
    float arg = k * r;

    /* H0(1)(z) = J0(z) + i Y0(z) */
    /* G = (i/4) * H0(1) = (i/4) * (J0 + i Y0) = (-Y0/4) + i (J0/4) */

    /* Using SE math functions for compatibility */
    float j0_val = (float)se_bessel_j0((double)arg);
    float y0_val = (float)se_bessel_y0((double)arg);

    *real = -y0_val / 4.0f;
    *imag = j0_val / 4.0f;
}

void se_green_func_3d_homogeneous(se_green_3d_params_t *params, float *real, float *imag)
{
    if (!params || !real || !imag) return;

    float v = params->v;
    float f = params->f;
    float dx = params->rec_x - params->src_x;
    float dy = params->rec_y - params->src_y;
    float dz = params->rec_z - params->src_z;

    float r = std::sqrt(dx * dx + dy * dy + dz * dz);

    /* Avoid singularity at r=0 */
    if (r < 1e-6f) {
        *real = 0.0f;
        *imag = 0.0f;
        return;
    }

    float omega = 2.0f * (float)PI * f;
    float k = omega / v;
    float arg = k * r;

    /* G = (1 / (4*pi*r)) * exp(i*k*r) */
    float amp = 1.0f / (4.0f * (float)PI * r);
    
    *real = amp * std::cos(arg);
    *imag = amp * std::sin(arg);
}