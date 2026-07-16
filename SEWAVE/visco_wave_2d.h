#ifndef VISCO_WAVE_2D_H
#define VISCO_WAVE_2D_H
#include "../SEmain/include/SEcore.h"
#include "../SEkissfft/se_kiss_fft.h"
#include "c_io.h"
#include "c_userdifine.h"


typedef struct visco_wave_2d_s
{
    int ns;  //炮数
    float fs; //第一炮位置
    float ds; //炮间距

    float sx; //震源位置
    float sz;
    float rz; //检波器深度
    int nx;
    int nz;
   
    float dxz;
    float dt;
    int nt;
    float fdom;
    float f0;
    float Omega0;

    char *Vpfile;
    char *Qfile;
    char *recordfile;

    float **vp;
    float **Q;
    float **init_record; //初始的震动，正演时为二维的子波分布，偏移时为地震记录
    float **record; //地震记录

    int nbc,L;
    float alpha; //吸收边界的衰减系数
    float fc;   // 截止频率 
    int order;    // 滤波器阶数
    
    int type; //0: P wave; 1: visco P wave;

    int myid;
    int np;

}visco_wave_2d_t;

void visco_wave_2d_getpar(visco_wave_2d_t *p,int argc,char **argv);

void visco_wave_2d_init(visco_wave_2d_t *p, int flag_homo, int flag_smooth);

void visco_wave_2d_free(visco_wave_2d_t *p);

void get_ricker_source(visco_wave_2d_t *p,float fdom);
void get_receiver_seis(visco_wave_2d_t *p, float **pr);


void visco_wave_2d_seismic_generate_TVM_PSandCS(visco_wave_2d_t *p, int type_compute_Laplace, char *PS_or_CS);
void visco_wave_2d_propagation_TVM_PSandCS(visco_wave_2d_t *p, int flag, int type_compute_Laplace, char *PS_or_CS, int amp_compensation_sign, float ***P_all_snapshots, float **Image);

void visco_wave_2d_seismic_generate_TVM_CS(visco_wave_2d_t *p);
void visco_wave_2d_propagation_TVM_CS(visco_wave_2d_t *p, int flag, int flag_denoise, int amp_compensation_sign, float ***P_all_snapshots, float **Image);


void caffarelli_extension_2d(float **input, float **output, int nx, int nz, float dxz);
void laplace_2d(float **input, float **output, int nx, int nz, float dxz, int order);

void apply_absorbing_boundary_conditions(float **P_next, float **P_current, float **P_past, visco_wave_2d_t *p);
void calculate_wavenumbers(float *kx, float *ky, int nx, int nz, float dxz);
void butterworth_lowpass_filter(float **H, float *kx, float *ky, int nx, int nz, float **vp,  float fc, int order);
void output_snapshot(int it, int tfix, int nx, int nz, float **P_current, char *fp_name);

void smooth2d(int n1, int n2, float r1, float r2, float **v);
void tripd(float *d, float *e, float *b, int n);
void matrix_multiply(float **a, float **b, float **c, int nx, int nz);


#endif