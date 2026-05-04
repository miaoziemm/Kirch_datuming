#ifndef SE_GREEN_FUNC_H
#define SE_GREEN_FUNC_H

#ifdef __cplusplus
extern "C" {
#endif

/* 2D Green's Function Parameters (Homogeneous) */
typedef struct se_green_2d_params_s {
    float v;            /* Velocity (m/s) */
    float f;            /* Frequency (Hz) */
    float src_x, src_z; /* Source coordinates */
    float rec_x, rec_z; /* Receiver coordinates */
} se_green_2d_params_t;

/* 3D Green's Function Parameters (Homogeneous) */
typedef struct se_green_3d_params_s {
    float v;                   /* Velocity (m/s) */
    float f;                   /* Frequency (Hz) */
    float src_x, src_y, src_z; /* Source coordinates */
    float rec_x, rec_y, rec_z; /* Receiver coordinates */
} se_green_3d_params_t;

/* Function Prototypes */

/**
 * @brief Computes 2D Green's function for homogeneous medium.
 *        G(r, w) = (i/4) * H0(1)(k*r)
 * 
 * @param params Pointer to parameter structure.
 * @param real Output real part.
 * @param imag Output imaginary part.
 */
void se_green_func_2d_homogeneous(se_green_2d_params_t *params, float *real, float *imag);

/**
 * @brief Computes 3D Green's function for homogeneous medium.
 *        G(r, w) = (1/(4*pi*r)) * exp(i*k*r)
 * 
 * @param params Pointer to parameter structure.
 * @param real Output real part.
 * @param imag Output imaginary part.
 */
void se_green_func_3d_homogeneous(se_green_3d_params_t *params, float *real, float *imag);

#ifdef __cplusplus
}
#endif

#endif /* SE_GREEN_FUNC_H */