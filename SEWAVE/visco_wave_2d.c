#include "visco_wave_2d.h"


void caffarelli_extension_2d(float **input, float **output, int nx, int nz, float dxz)
{
    int i, j, k;

    // homo模型下优化参数,第一组精度较高参数！！，但是有点慢
    // int Nz = 10;
    // float H = 100.0;
    // float s = 0.5;
    // float relax = 1.2;
    // int max_iter = 30;
    // float tol = 1e-5;


    // 获取参数值, marmousi模型下优化参数,精度高！！
    int Nz = 20;
    float H = 90.0;
    float s = 0.5;
    float relax = 1.8;
    int max_iter = 20;
    float tol = 1e-5;


    // 计算延拓方程参数
    float alpha = 1.0f - 2.0f * s;
    
    // 预计算常量，避免重复计算
    float dxz_sq = dxz * dxz;
    float inv_dxz_sq = 1.0f / dxz_sq;

    // 均匀z网格（0到H）
    float *z = (float *)malloc(Nz * sizeof(float));
    float dz = H / (Nz - 1);
    float inv_dz_sq = 1.0f / (dz * dz);
    float inv_2dz = 1.0f / (2.0f * dz);
    
    for (k = 0; k < Nz; ++k)
    {
        z[k] = k * dz;
    }

    // 预计算z相关参数
    float *z_alpha = (float *)malloc(Nz * sizeof(float));
    for (k = 0; k < Nz; ++k)
    {
        z_alpha[k] = powf(z[k], alpha);
    }
    float *z_half = (float *)malloc((Nz - 1) * sizeof(float));
    float *z_half_alpha = (float *)malloc((Nz - 1) * sizeof(float));
    for (k = 0; k < Nz - 1; ++k)
    {
        z_half[k] = 0.5f * (z[k] + z[k + 1]);
        z_half_alpha[k] = powf(z_half[k], alpha);
    }

    // 假设 input 是大小为 nx*nz 的二维数组，output 也是
    // 初始化三维解数组 U[ny][nx][Nz]
    int ny = nx, nx_ = nz; // 这里假设 input 为 ny*nx，用户需根据实际情况调整
    float ***U = (float ***)malloc(ny * sizeof(float **));
    for (i = 0; i < ny; ++i)
    {
        U[i] = (float **)malloc(nx_ * sizeof(float *));
        for (j = 0; j < nx_; ++j)
        {
            U[i][j] = (float *)calloc(Nz, sizeof(float));
        }
    }
    // z=0边界
    for (i = 0; i < ny; ++i)
        for (j = 0; j < nx_; ++j)
            U[i][j][0] = input[i][j];

    // SOR迭代
    float residual = 1e10f;
    int iter = 0;

    while (residual > tol && iter < max_iter)
    {
        residual = 0.0f;
        #pragma omp parallel for collapse(3) reduction(max:residual)
        for (k = 1; k < Nz - 1; ++k)
        {
            
            for (i = 1; i < ny - 1; ++i)
            {
                for (j = 1; j < nx_ - 1; ++j)
                {
                    // 预计算当前k层的常量
            float c_xy = z_alpha[k] * inv_dxz_sq;
            float c_z = inv_dz_sq;
            float z_half_alpha_k = z_half_alpha[k];
            float z_half_alpha_k_minus_1 = z_half_alpha[k - 1];
            float denom_const = c_z * (z_half_alpha_k + z_half_alpha_k_minus_1);
            float four_c_xy = 4.0f * c_xy;
                    // 拉普拉斯项 (x,y方向)
                    float lap_xy = (U[i + 1][j][k] + U[i - 1][j][k] +
                                    U[i][j + 1][k] + U[i][j - 1][k] -
                                    4.0f * U[i][j][k]);

                    // z方向导数 (修正的离散格式)
                    float z_term = z_half_alpha_k * (U[i][j][k + 1] - U[i][j][k]) - 
                                   z_half_alpha_k_minus_1 * (U[i][j][k] - U[i][j][k - 1]);

                    // 残差计算
                    float R = c_xy * lap_xy + c_z * z_term;

                    // SOR更新
                    float denom = four_c_xy + denom_const;
                    U[i][j][k] = U[i][j][k] + relax * R / denom;

                    // 更新残差
                    float abs_R = fabsf(R);
                    if (abs_R > residual)
                        residual = abs_R;
                }
            }
        }
        iter++;
        // 可选：打印残差
        // if (iter % 50 == 0) printf("Iter %d, Residual: %.4e\n", iter, residual);
    }

    // 预计算归一化系数
    float C_s = powf(2.0f, 2.0f * s - 1.0f) * tgammaf(s) / tgammaf(1.0f - s);
    float coeff_C_s_inv_2dz = C_s * inv_2dz;

    // 计算z=0处的法向导数（二阶精度）
    if (Nz >= 3)
    {
        for (i = 0; i < ny; ++i)
        {
            for (j = 0; j < nx_; ++j)
            {
                output[i][j] = -(4.0f * U[i][j][1] - 3.0f * U[i][j][0] - U[i][j][2]) * coeff_C_s_inv_2dz;
            }
        }
    }
    else
    {
        printf("Nz must be at least 3 to compute normal derivative\n");
    }

    // 释放内存
    for (i = 0; i < ny; ++i)
    {
        for (j = 0; j < nx_; ++j)
        {
            free(U[i][j]);
        }
        free(U[i]);
    }
    free(U);
    free(z);
    free(z_alpha);
    free(z_half);
    free(z_half_alpha);
}



void laplace_2d(float **input, float **output, int nx, int nz, float dxz, int order)
{
    int ix, iz;
    const float dxz_sq = dxz * dxz;

    // 清零
    #pragma omp parallel for collapse(2) schedule(static)
    for (iz = 0; iz < nz; ++iz)
        for (ix = 0; ix < nx; ++ix)
            output[ix][iz] = 0.0f;

    int boundary = order / 2 + 1;
    if (boundary >= nx/2 || boundary >= nz/2) {
        printf("Error: Grid size too small for order %d difference\n", order);
        return;
    }

    switch (order) {
        case 2:
            #pragma omp parallel for schedule(static) collapse(2)
            for (iz = boundary; iz < nz - boundary; ++iz) {
       
                for (ix = boundary; ix < nx - boundary; ++ix) {
                    output[ix][iz] =
                        (input[ix+1][iz] - 2.0f*input[ix][iz] + input[ix-1][iz]) / dxz_sq +
                        (input[ix][iz+1] - 2.0f*input[ix][iz] + input[ix][iz-1]) / dxz_sq;
                }
            }
            break;

        case 4: {
            const float c = 1.0f/(12.0f*dxz_sq);
            #pragma omp parallel for schedule(static) collapse(2)
            for (iz = boundary; iz < nz - boundary; ++iz) {
             
                for (ix = boundary; ix < nx - boundary; ++ix) {
                    output[ix][iz] =
                        c * (-input[ix-2][iz] + 16*input[ix-1][iz] - 30*input[ix][iz] + 16*input[ix+1][iz] - input[ix+2][iz]) +
                        c * (-input[ix][iz-2] + 16*input[ix][iz-1] - 30*input[ix][iz] + 16*input[ix][iz+1] - input[ix][iz+2]);
                }
            }
            break;
        }

        case 6: {
            const float c = 1.0f/(180.0f*dxz_sq);
            #pragma omp parallel for schedule(static) collapse(2)
            for (iz = boundary; iz < nz - boundary; ++iz) {
    
                for (ix = boundary; ix < nx - boundary; ++ix) {
                    output[ix][iz] =
                        c * (  2*input[ix-3][iz] - 27*input[ix-2][iz] + 270*input[ix-1][iz]
                             -490*input[ix][iz]   + 270*input[ix+1][iz] - 27*input[ix+2][iz]
                              + 2*input[ix+3][iz]) +
                        c * (  2*input[ix][iz-3] - 27*input[ix][iz-2] + 270*input[ix][iz-1]
                             -490*input[ix][iz]   + 270*input[ix][iz+1] - 27*input[ix][iz+2]
                              + 2*input[ix][iz+3]);
                }
            }
            break;
        }

        case 8: {
            const float c0 = (-1.0f/560.0f) / dxz_sq;
            const float c1 = ( 8.0f/315.0f) / dxz_sq;
            const float c2 = ( -1.0f/5.0f ) / dxz_sq;
            const float c3 = (  8.0f/5.0f ) / dxz_sq;
            const float c4 = (-205.0f/72.0f)/ dxz_sq;
            #pragma omp parallel for schedule(static) collapse(2)
            for (iz = boundary; iz < nz - boundary; ++iz) {
      
                for (ix = boundary; ix < nx - boundary; ++ix) {
                    output[ix][iz] =
                        ( c0*input[ix-4][iz] + c1*input[ix-3][iz] + c2*input[ix-2][iz] + c3*input[ix-1][iz]
                        + c4*input[ix][iz]   + c3*input[ix+1][iz] + c2*input[ix+2][iz] + c1*input[ix+3][iz]
                        + c0*input[ix+4][iz] ) +
                        ( c0*input[ix][iz-4] + c1*input[ix][iz-3] + c2*input[ix][iz-2] + c3*input[ix][iz-1]
                        + c4*input[ix][iz]   + c3*input[ix][iz+1] + c2*input[ix][iz+2] + c1*input[ix][iz+3]
                        + c0*input[ix][iz+4] );
                }
            }
            break;
        }

        default:
            printf("Error: Unsupported difference order %d. Supported orders: 2, 4, 6, 8\n", order);

    }
}

// void calculate_wavenumbers(float *kx, float *ky, int nx, int nz, float dxz) 
// {
//     // 计算波数k
//     int i;
//     float kmax = PI / dxz;
//     float delta_kx = kmax / (nx / 2);
//     float delta_ky = kmax / (nz / 2);

//     for (i = 0; i < nx / 2; i++) {
//         kx[i] = (i + 1) * delta_kx;
//         kx[nx / 2 + i] = -kmax + (i + 1) * delta_kx;
//     }
//     for (i = 0; i < nz / 2; i++) {
//         ky[i] = (i + 1) * delta_ky;
//         ky[nz / 2 + i] = -kmax + (i + 1) * delta_ky;
//     }
// }

void calculate_wavenumbers(float *kx, float *kz, int nx, int nz, float dxz)
{
    int i;
    float dkx = 2.0f * PI / (nx * dxz);
    float dkz = 2.0f * PI / (nz * dxz);

    for (i = 0; i < nx; ++i) {
        if (i <= nx / 2) kx[i] = i * dkx;
        else             kx[i] = (i - nx) * dkx;
    }

    for (i = 0; i < nz; ++i) {
        if (i <= nz / 2) kz[i] = i * dkz;
        else             kz[i] = (i - nz) * dkz;
    }
}


void butterworth_lowpass_filter(float **H, float *kx, float *ky, int nx, int nz, float **vp,  float fc, int order) 
{
    int i, j, ix, iz;
    // 计算最大速度vmax
    float vmax = 0;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < nz; j++) {
            if (vp[i][j] > vmax) {
                vmax = vp[i][j];
            }
        }
    }
    
    float kc = fc * 2 * PI / vmax; // 截止波数

    for (ix = 0; ix < nx; ix++) {
        for (iz = 0; iz < nz; iz++) {
            float K2 = kx[ix] * kx[ix] + ky[iz] * ky[iz];
            H[ix][iz] = 1.0 / (1.0 + pow((sqrt(K2) / kc), 2 * order)); // 构造 Butterworth低通滤波器
        }
    }

    // write_binary_file_2d("k_matrix.bin", k_matrix, p->nx, p->nz); // 输出k_matrix到二进制文件

    // 将H进行fftshift
    // fftshift_2d(H, p->nx, p->nz);

    // write_binary_file_2d("H_new.bin", H, p->nx, p->nz); // 输出滤波器到二进制文件
    // exit(0);
}

void output_snapshot(int it, int tfix, int nx, int nz, float **P_current, char *fp_name) 
{
    if (it == tfix) {
        printf("Snap Output!! nx=%d,nz=%d\n", nx, nz);
        FILE *fp = fopen(fp_name, "wb");
       
        if (fp == NULL) 
        {
            err("error in opening file snap");
        }

        for (int i = 0; i < nx; i++) 
        {
            for (int j = 0; j < nz; j++) 
            {
                fwrite(&P_current[i][j], sizeof(float), 1, fp);
            }
        }
        fclose(fp);
    }
}


void apply_absorbing_boundary_conditions(float **P_next, float **P_current, float **P_past, visco_wave_2d_t *p)
{
    int ix, iz;
    float w, s, t1, t2, t3;

    for (ix = 0; ix < p->nx; ix++)
    {
        for (iz = 0; iz < p->nz; iz++)
        {
            // 吸收边界
            // left
            if (ix >= 0 && ix < p->L - 4 && iz >= 0 && iz < p->nz)
            {
                w = 1 - 1.0 * ix / p->L;
                s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                t1 = (2 - s) * (1 - s) / 2;
                t2 = s * (2 - s);
                t3 = s * (s - 1) / 2;

                P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                    t2 * P_current[ix + 1][iz] + t3 * P_current[ix + 2][iz]) +
                                        (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                    2 * t1 * t2 * P_past[ix + 1][iz] +
                                                    (2 * t1 * t3 + t2 * t2) * P_past[ix + 2][iz] +
                                                    2 * t2 * t3 * P_past[ix + 3][iz] +
                                                    t3 * t3 * P_past[ix + 4][iz])) +
                                    (1 - w) * P_next[ix][iz];
            }

            // right
            if (ix > p->nx - p->L + 3 && ix < p->nx && iz >= 0 && iz < p->nz)
            {
                w = 1 - 1.0 * (p->nx - ix - 1) / p->L;
                s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                t1 = (2 - s) * (1 - s) / 2;
                t2 = s * (2 - s);
                t3 = s * (s - 1) / 2;

                P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                    t2 * P_current[ix - 1][iz] + t3 * P_current[ix - 2][iz]) +
                                        (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                    2 * t1 * t2 * P_past[ix - 1][iz] +
                                                    (2 * t1 * t3 + t2 * t2) * P_past[ix - 2][iz] +
                                                    2 * t2 * t3 * P_past[ix - 3][iz] +
                                                    t3 * t3 * P_past[ix - 4][iz])) +
                                    (1 - w) * P_next[ix][iz];
            }

            // top
            if (iz >= 0 && iz < p->L - 4 && ix >= 0 && ix < p->nx)
            {
                w = 1 - 1.0 * iz / p->L;
                s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                t1 = (2 - s) * (1 - s) / 2;
                t2 = s * (2 - s);
                t3 = s * (s - 1) / 2;

                P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                    t2 * P_current[ix][iz + 1] + t3 * P_current[ix][iz + 2]) +
                                        (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                    2 * t1 * t2 * P_past[ix][iz + 1] +
                                                    (2 * t1 * t3 + t2 * t2) * P_past[ix][iz + 2] +
                                                    2 * t2 * t3 * P_past[ix][iz + 3] +
                                                    t3 * t3 * P_past[ix][iz + 4])) +
                                    (1 - w) * P_next[ix][iz];
            }

            // bottom
            if (iz > p->nz - p->L + 3 && iz < p->nz && ix >= 0 && ix < p->nx)
            {
                w = 1 - 1.0 * (p->nz - iz - 1) / p->L;
                s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                t1 = (2 - s) * (1 - s) / 2;
                t2 = s * (2 - s);
                t3 = s * (s - 1) / 2;

                P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                    t2 * P_current[ix][iz - 1] + t3 * P_current[ix][iz - 2]) +
                                        (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                    2 * t1 * t2 * P_past[ix][iz - 1] +
                                                    (2 * t1 * t3 + t2 * t2) * P_past[ix][iz - 2] +
                                                    2 * t2 * t3 * P_past[ix][iz - 3] +
                                                    t3 * t3 * P_past[ix][iz - 4])) +
                                    (1 - w) * P_next[ix][iz];
            }
        }
    }
}



void smooth2d(int n1, int n2, float r1, float r2, float **v)
{

    int nmax, ix, iz;

    float **w;
    float *d, *e, *f;
    int *win;
    float rw = 0.0;

    r1 = r1 * r1 * 0.25;
    r2 = r2 * r2 * 0.25;

    /* allocate space */
    nmax = (n1 < n2) ? n2 : n1;
    win = alloc1int(4);
    //	v = alloc2float(n1,n2);
    //	v0 = alloc2float(n1,n2);
    w = alloc2float(n1, n2);
    d = alloc1float(nmax);
    e = alloc1float(nmax);
    f = alloc1float(nmax);
    rw = rw * rw * 0.25;

    win[0] = 0;
    win[1] = n1;
    win[2] = 0;
    win[3] = n2;

    /* define the window function */
    for (ix = 0; ix < n2; ++ix)
        for (iz = 0; iz < n1; ++iz)
            w[ix][iz] = 0;
    for (ix = win[2]; ix < win[3]; ++ix)
        for (iz = win[0]; iz < win[1]; ++iz)
            w[ix][iz] = 1;

    if (win[0] > 0 || win[1] < n1 || win[2] > 0 || win[3] < n2)
    {
        /*	smooth the window function */
        for (iz = 0; iz < n1; ++iz)
        {
            for (ix = 0; ix < n2; ++ix)
            {
                d[ix] = 1.0 + 2.0 * rw;
                e[ix] = -rw;
                f[ix] = w[ix][iz];
            }
            d[0] -= rw;
            d[n2 - 1] -= rw;
            tripd(d, e, f, n2);
            for (ix = 0; ix < n2; ++ix)
                w[ix][iz] = f[ix];
        }
        for (ix = 0; ix < n2; ++ix)
        {
            for (iz = 0; iz < n1; ++iz)
            {
                d[iz] = 1.0 + 2.0 * rw;
                e[iz] = -rw;
                f[iz] = w[ix][iz];
            }
            d[0] -= rw;
            d[n1 - 1] -= rw;
            tripd(d, e, f, n1);
            for (iz = 0; iz < n1; ++iz)
                w[ix][iz] = f[iz];
        }
    }

    /* solving for the smoothing velocity */
    for (iz = 0; iz < n1; ++iz)
    {
        for (ix = 0; ix < n2 - 1; ++ix)
        {
            d[ix] = 1.0 + r2 * (w[ix][iz] + w[ix + 1][iz]);
            e[ix] = -r2 * w[ix + 1][iz];
            f[ix] = v[ix][iz];
        }
        d[0] -= r2 * w[0][iz];
        d[n2 - 1] = 1.0 + r2 * w[n2 - 1][iz];
        f[n2 - 1] = v[n2 - 1][iz];
        tripd(d, e, f, n2);
        for (ix = 0; ix < n2; ++ix)
            v[ix][iz] = f[ix];
    }
    for (ix = 0; ix < n2; ++ix)
    {
        for (iz = 0; iz < n1 - 2; ++iz)
        {
            d[iz] = 1.0 + r1 * (w[ix][iz + 1] + w[ix][iz + 2]);
            e[iz] = -r1 * w[ix][iz + 2];
            f[iz] = v[ix][iz + 1];
        }
        f[0] += r1 * w[ix][1] * v[ix][0];
        d[n1 - 2] = 1.0 + r1 * w[ix][n1 - 1];
        f[n1 - 2] = v[ix][n1 - 1];
        tripd(d, e, f, n1 - 1);
        for (iz = 0; iz < n1 - 1; ++iz)
            v[ix][iz + 1] = f[iz];
    }

    free2float(w);
    free1float(d);
    free1float(e);
    free1float(f);
    free1int(win);
}

void tripd(float *d, float *e, float *b, int n)
/*****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b*/
{
    int k;
    float temp;

    /* decomposition */
    for (k = 1; k < n; ++k)
    {
        temp = e[k - 1];
        e[k - 1] = temp / d[k - 1];
        d[k] -= temp * e[k - 1];
    }

    /* substitution	*/
    for (k = 1; k < n; ++k)
        b[k] -= e[k - 1] * b[k - 1];

    b[n - 1] /= d[n - 1];
    for (k = n - 1; k > 0; --k)
        b[k - 1] = b[k - 1] / d[k - 1] - e[k - 1] * b[k];
}

void matrix_multiply(float **a, float **b, float **c, int nx, int nz)
{
    int i, j, k;
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < nz; j++)
        {
            c[i][j] = a[i][j] * b[i][j];
        }
    }
}

void visco_wave_2d_getpar(visco_wave_2d_t *p, int argc, char **argv)
{
    initargs(argc, argv);

    if (!getparstring("Vpfile", &p->Vpfile))
        err("error in getting Vpfile!");
    if (!getparstring("Qfile", &p->Qfile))
        err("error in getting Qfile!");
    if (!getparstring("Recordfile", &p->recordfile))
        p->recordfile = "Recordfile";
    if (!getparint("nx", &p->nx))
        p->nx = 100;
    if (!getparint("nz", &p->nz))
        p->nz = 100;
    if (!getparfloat("dxz", &p->dxz))
        p->dxz = 10.0;
    if (!getparint("ns", &p->ns))
        p->ns = 1;
    if (!getparfloat("fs", &p->fs))
        p->fs = p->nx / 2 * p->dxz;
    if (!getparfloat("ds", &p->ds))
        p->ds = p->nx / 2 * p->dxz;
    if (!getparfloat("sz", &p->sz))
        p->sz = p->nz / 2 * p->dxz;
    if (!getparfloat("rz", &p->rz))
        p->rz = p->sz;
    if (!getparfloat("dt", &p->dt))
        p->dt = 0.001;
    if (!getparint("nt", &p->nt))
        p->nt = 100;
    if (!getparfloat("fdom", &p->fdom))
        p->fdom = 20.0;
    if (!getparfloat("f0", &p->f0))
        p->f0 = 20.0;
    if (!getparfloat("Omega0", &p->Omega0))
        p->Omega0 = 2 * PI * p->f0;
    if (!getparint("type", &p->type))
        p->type = 1;
    if (!getparint("L", &p->L))
        p->L = 30;
    if (!getparint("nbc", &p->nbc))
        p->nbc = 100;
    if (!getparfloat("alpha", &p->alpha))
        p->alpha = 1.0;
    if (!getparfloat("fc", &p->fc))
        p->fc = 120;
    if (!getparint("order", &p->order))
        p->order = 2;
    if (!getparint("np", &p->np))
        p->np = 1;
}

void visco_wave_2d_init(visco_wave_2d_t *p, int flag_homo, int flag_smooth)
{
    int ix, iz;
    p->vp = alloc2float((p->nz + 2 * p->nbc), (p->nx + 2 * p->nbc));
    p->Q = alloc2float((p->nz + 2 * p->nbc), (p->nx + 2 * p->nbc));
    p->init_record = alloc2float((p->nt), (p->nx + 2 * p->nbc));
    p->record = alloc2float((p->nt), (p->nx + 2 * p->nbc));

    // p->sx = p->sx + p->nbc * p->dxz;
    // p->sz = p->sz + p->nbc * p->dxz;
    // p->rz = p->rz + p->nbc * p->dxz;

    p->myid = 0;

    memset(p->vp[0], 0, sizeof(float) * (p->nx + 2 * p->nbc) * (p->nz + 2 * p->nbc));
    memset(p->Q[0], 0, sizeof(float) * (p->nx + 2 * p->nbc) * (p->nz + 2 * p->nbc));
    memset(p->init_record[0], 0, sizeof(float) * (p->nt) * (p->nx + 2 * p->nbc));
    memset(p->record[0], 0, sizeof(float) * (p->nx + 2 * p->nbc) * (p->nt));

    FILE *fp;
    if ((fp = fopen(p->Vpfile, "rb")) == NULL)
        err("error in opening Vpfile!");

    for (ix = p->nbc; ix < p->nx + p->nbc; ix++)
    {
        for (iz = p->nbc; iz < p->nz + p->nbc; iz++)
        {
            fread(&p->vp[ix][iz], sizeof(float), 1, fp);
        }
    }
    fclose(fp);

    for (ix = 0; ix < p->nbc; ix++)
    {
        for (iz = 0; iz < p->nz + 2 * p->nbc; iz++)
        {
            p->vp[ix][iz] = p->vp[p->nbc][iz];
            p->vp[ix + p->nx + p->nbc][iz] = p->vp[p->nx + p->nbc - 1][iz];
        }
    }

    for (iz = 0; iz < p->nbc; iz++)
    {
        for (ix = 0; ix < p->nx + 2 * p->nbc; ix++)
        {
            p->vp[ix][iz] = p->vp[ix][p->nbc];
            p->vp[ix][iz + p->nz + p->nbc] = p->vp[ix][p->nz + p->nbc - 1];
        }
    }

    // write_binary_file_2d("V_test", p->vp, p->nx + 2 * p->nbc, p->nz + 2 * p->nbc);
    // exit(0);

    if ((fp = fopen(p->Qfile, "rb")) == NULL)
        err("error in opening Qfile!");

    for (ix = p->nbc; ix < p->nx + p->nbc; ix++)
    {
        for (iz = p->nbc; iz < p->nz + p->nbc; iz++)
        {
            fread(&p->Q[ix][iz], sizeof(float), 1, fp);
        }
    }
    fclose(fp);

    for (ix = 0; ix < p->nbc; ix++)
    {
        for (iz = 0; iz < p->nz + 2 * p->nbc; iz++)
        {
            p->Q[ix][iz] = p->Q[p->nbc][iz];
            p->Q[ix + p->nx + p->nbc][iz] = p->Q[p->nx + p->nbc - 1][iz];
        }
    }

    for (iz = 0; iz < p->nbc; iz++)
    {
        for (ix = 0; ix < p->nx + 2 * p->nbc; ix++)
        {
            p->Q[ix][iz] = p->Q[ix][p->nbc];
            p->Q[ix][iz + p->nz + p->nbc] = p->Q[ix][p->nz + p->nbc - 1];
        }
    }

    SE_MESSAGE("Origin model size : (%d, %d)\n", p->nx, p->nz);

    p->nx = p->nx + 2 * p->nbc;
    p->nz = p->nz + 2 * p->nbc;

    SE_MESSAGE("Extend model size : (%d, %d)\n", p->nx, p->nz);

    if (flag_smooth == 1)
    {
        smooth2d(p->nz,p->nx, 5, 5, p->vp);
        smooth2d(p->nz,p->nx, 5, 5, p->Q);
    }
    

    if (flag_homo == 1)
    {
        // 均匀模型
        int isx = (int)(p->sx);
        int isz = (int)(p->sz);
        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {
                // // 自主设置
                // p->vp[ix][iz] = 2500;
                // p->Q[ix][iz] = 20;

                // 设置为震源位置处的值
                p->vp[ix][iz] = p->vp[isx][isz];
                p->Q[ix][iz] = p->Q[isx][isz];
            }
        }
    }
}

void visco_wave_2d_free(visco_wave_2d_t *p)
{
    free2float(p->vp);
    free2float(p->Q);
    free2float(p->init_record);
    free2float(p->record);
}

void get_ricker_source(visco_wave_2d_t *p, float fdom)
{

    // 每一炮前必须清空激发记录，否则多炮循环时前一炮的震源会残留
    memset(p->init_record[0], 0, sizeof(float) * p->nt * p->nx);

    // 生成一组雷克子波，放在sx对应的位置
    int i;
    int isx = (int)(p->sx);
    float t;
    float delay = 1.5 / fdom;
    for (i = 0; i < p->nt; i++)
    {
        t = (i)*p->dt - delay;
        p->init_record[isx][i] = (1 - 2 * (PI * fdom * t) * (PI * fdom * t)) * exp(-(PI * fdom * t) * (PI * fdom * t));
    }
}

void get_receiver_seis(visco_wave_2d_t *p, float **pr)
{

    // 每一炮反传前先清空记录，避免上一炮或上一阶段的记录残留
    memset(p->init_record[0], 0, sizeof(float) * p->nt * p->nx);

    // 读取地震记录，放在sz对应的位置
    int i, j;
    for (i = 0; i < p->nt; i++)
    {
        for (j = p->nbc; j < p->nx - p->nbc; j++)
        {
            p->init_record[j][i] = pr[j - p->nbc][p->nt - i - 1];
        }
    }

}



void visco_wave_2d_seismic_generate_TVM_PSandCS(visco_wave_2d_t *p, int type_compute_Laplace, char *PS_or_CS)
{

   int i, j, k;
    int it, ix, iz;
    float K2;

    float c0, c1, c2;
    c0 = 411.7 + 4.36;
    c1 = -51.64;
    c2 = 4.366 / 2;

    int isx = (int)(p->sx);
    int isz = (int)(p->sz);
    int irz = (int)(p->rz);

    float **C1, **C2, **C3, **C4;
    C1 = alloc2float(p->nz, p->nx);
    C2 = alloc2float(p->nz, p->nx);
    C3 = alloc2float(p->nz, p->nx);
    C4 = alloc2float(p->nz, p->nx);

    for (i = 0; i < p->nx; i++)
    {
        for (j = 0; j < p->nz; j++)
        {
            C1[i][j] = (1 - 4 * c2 / (PI * p->Q[i][j]) + 2 * log(p->Omega0) / PI / p->Q[i][j])/ (p->vp[i][j] * p->vp[i][j]);
            C2[i][j] = 2 * c1 / (PI * p->Q[i][j]* p->vp[i][j]);
            C3[i][j] = 2 * (c0 - c2) / (PI * p->Q[i][j]* p->vp[i][j] * p->vp[i][j]);
            C4[i][j] = 1 / (p->vp[i][j] * p->Q[i][j]);
        }
    }
   

    float **P_current, **P_past, **P_next, **P_1t;

    se_complex **P_fft, **P_fft_temp, **P_1t_fft;

    P_current = alloc2float(p->nz, p->nx);
    P_past = alloc2float(p->nz, p->nx);
    P_next = alloc2float(p->nz, p->nx);
    P_1t = alloc2float(p->nz, p->nx);

    // 置0初始化
    memset(P_current[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_past[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_next[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_1t[0], 0, sizeof(float) * p->nx * p->nz);


    P_fft = alloc2secomplex(p->nz, p->nx);
    P_fft_temp = alloc2secomplex(p->nz, p->nx);
    P_1t_fft = alloc2secomplex(p->nz, p->nx);


    memset(P_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_fft_temp[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_1t_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);

    float **P_la = alloc2float(p->nz, p->nx);
    se_complex **P_la_fft = alloc2secomplex(p->nz, p->nx);

    float **amp = alloc2float(p->nz, p->nx);
    se_complex **amp_fft = alloc2secomplex(p->nz, p->nx);

    float **KP = alloc2float(p->nz, p->nx);
    se_complex **KP_fft = alloc2secomplex(p->nz, p->nx);


    // 置0初始化
    memset(P_la[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_la_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(amp[0], 0, sizeof(float) * p->nx * p->nz);
    memset(amp_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(KP[0], 0, sizeof(float) * p->nx * p->nz);
    memset(KP_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);


    // 计算波数k
    float *kx = (float *)malloc(p->nx * sizeof(float));
    float *ky = (float *)malloc(p->nz * sizeof(float));
    calculate_wavenumbers(kx, ky, p->nx, p->nz, p->dxz);


    for (it = 0; it < p->nt; it++)
    {

        P_current[isx][isz] = P_current[isx][isz] + p->init_record[isx][it]; // 添加雷克子波震源

        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {
                P_1t[ix][iz] = (P_current[ix][iz] - P_past[ix][iz]) / p->dt;
            }
        }

        // 将P_current放到P_fft的实部
        for (i = 0; i < p->nx; i++)
        {
            for (j = 0; j < p->nz; j++)
            {
                P_fft[i][j].r = P_current[i][j];
                P_fft[i][j].i = 0.0;
                P_1t_fft[i][j].r = P_1t[i][j];
                P_1t_fft[i][j].i = 0.0;
            }
        }

        // 进行二维FFT
        se_fftw_fft_2d_quick(P_fft, P_fft_temp, p->nx, p->nz, 0);
        memcpy(P_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
        se_fftw_fft_2d_quick(P_1t_fft, P_fft_temp, p->nx, p->nz, 0);
        memcpy(P_1t_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

        //伪谱或差分计算Laplace
        if (type_compute_Laplace == 0)
        {
            //有限差分八阶计算Laplace
            laplace_2d(P_current, P_la, p->nx, p->nz, p->dxz, 8);
        }
        else if (type_compute_Laplace == 1)
        {
            // 伪谱法计算Laplace
            for (ix = 0; ix < p->nx; ix++)
            {
                for (iz = 0; iz < p->nz; iz++)
                {
                    K2 = kx[ix] * kx[ix] + ky[iz] * ky[iz];
                    P_la_fft[ix][iz].r = -K2 * P_fft[ix][iz].r;
                    P_la_fft[ix][iz].i = -K2 * P_fft[ix][iz].i;

                }
            }
            se_fftw_fft_2d_quick(P_la_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(P_la_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

            for (i = 0; i < p->nx; i++)
            {
                for (j = 0; j < p->nz; j++)
                {
                    P_la[i][j] = P_la_fft[i][j].r;
                }
            }
        }
        else
        {
            err("error in type_compute_Laplace!");
        }


   
        if (strcmp(PS_or_CS, "PS") == 0)
        {
            // 伪谱法计算amp和KP
            for (ix = 0; ix < p->nx; ix++)
            {
                for (iz = 0; iz < p->nz; iz++)
                {
                    K2 = kx[ix] * kx[ix] + ky[iz] * ky[iz];

                    amp_fft[ix][iz].r = sqrt(K2) * P_1t_fft[ix][iz].r;
                    amp_fft[ix][iz].i = sqrt(K2) * P_1t_fft[ix][iz].i;

                    KP_fft[ix][iz].r = sqrt(K2) * P_fft[ix][iz].r;
                    KP_fft[ix][iz].i = sqrt(K2) * P_fft[ix][iz].i;

                }
            }

            se_fftw_fft_2d_quick(KP_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(KP_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
            se_fftw_fft_2d_quick(amp_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(amp_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

            for (i = 0; i < p->nx; i++)
            {
                for (j = 0; j < p->nz; j++)
                {
                    amp[i][j] = amp_fft[i][j].r;
                    KP[i][j] = KP_fft[i][j].r;
                }
            }
        }
        else if (strcmp(PS_or_CS, "CS") == 0)
        {
            // CS延拓计算amp和KP
            caffarelli_extension_2d(P_1t,amp,p->nx,p->nz,p->dxz);
            caffarelli_extension_2d(P_current,KP,p->nx,p->nz,p->dxz);

        }

        

        if (it == 10)
        {
            // 打印p->type
            printf("p->type = %d\n", p->type);
        }

        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {

                if (p->type == 1)
                {
                    // 解耦的Aki模型波场传播，利用切比雪夫近似
                    P_next[ix][iz] = 2 * P_current[ix][iz] - P_past[ix][iz] + (pow(p->dt, 2) / C1[ix][iz]) * (C4[ix][iz]* (-amp[ix][iz]) + P_la[ix][iz] - C2[ix][iz] * KP[ix][iz] - C3[ix][iz] * P_current[ix][iz]);
                }
                else if (p->type == 0)
                {
                    // 声波方程
                    P_next[ix][iz] = 2 * P_current[ix][iz] - P_past[ix][iz] + pow(p->vp[ix][iz], 2) * pow(p->dt, 2) * P_la[ix][iz];
                }

                else
                {
                    err("error in type!");
                }
                
            }
        }

        // 吸收边界调用函数
        apply_absorbing_boundary_conditions(P_next, P_current, P_past, p);


        memcpy(P_past[0], P_current[0], sizeof(float) * p->nx * p->nz);
        memcpy(P_current[0], P_next[0], sizeof(float) * p->nx * p->nz);

        // 记录地震数据
        for (i = p->nbc; i < p->nx - p->nbc; i++)
        {
            p->record[i - p->nbc][it] = P_current[i][irz];
        }


        if (it % 300 == 0)
        {
            SE_MESSAGE("Thread ID : %d :: Processing %d/%d\n", p->myid, it, p->nt);
        }


//        // 过程波场输出
//        int tfix;
//        tfix=500;
//        output_snapshot(it, tfix, p->nx, p->nz, P_current, "snap500.bin");
//
//        tfix=830;
//        output_snapshot(it, tfix, p->nx, p->nz, P_current, "snap830.bin");
//
//        tfix=850;
//        output_snapshot(it, tfix, p->nx, p->nz, P_current, "snap850.bin");

    }



    free(kx);
    free(ky);

    free2float(amp);
    free2secomplex(amp_fft);
    free2float(P_current);
    free2float(P_past);
    free2float(P_next);

    free2secomplex(P_fft);
    free2secomplex(P_fft_temp);

    free2secomplex(P_1t_fft);
    free2float(P_1t);

    free2float(C1);
    free2float(C2);
    free2float(C3);
    free2float(C4);

    free2float(KP);
    free2secomplex(KP_fft);
    free2float(P_la);
    free2secomplex(P_la_fft);

    // free2float(damp);


}

void visco_wave_2d_propagation_TVM_PSandCS(visco_wave_2d_t *p, int flag, int type_compute_Laplace, char *PS_or_CS, int amp_compensation_sign, float ***P_all_snapshots, float **Image)
{

    int i, j, k;
    int it, ix, iz;
    float K2;

    float c0, c1, c2;
    c0 = 411.7 + 4.36;
    c1 = -51.64;
    c2 = 4.366 / 2;

    int isx = (int)(p->sx);
    int isz = (int)(p->sz);
    int irz = (int)(p->rz);

    float **C1, **C2, **C3, **C4;
    C1 = alloc2float(p->nz, p->nx);
    C2 = alloc2float(p->nz, p->nx);
    C3 = alloc2float(p->nz, p->nx);
    C4 = alloc2float(p->nz, p->nx);

    for (i = 0; i < p->nx; i++)
    {
        for (j = 0; j < p->nz; j++)
        {
            C1[i][j] = (1 - 4 * c2 / (PI * p->Q[i][j]) + 2 * log(p->Omega0) / PI / p->Q[i][j])/ (p->vp[i][j] * p->vp[i][j]);
            C2[i][j] = 2 * c1 / (PI * p->Q[i][j]* p->vp[i][j]);
            C3[i][j] = 2 * (c0 - c2) / (PI * p->Q[i][j]* p->vp[i][j] * p->vp[i][j]);
            C4[i][j] = 1 / (p->vp[i][j] * p->Q[i][j]);
        }
    }
   

    float **P_current, **P_past, **P_next, **P_1t;

    se_complex **P_fft, **P_fft_temp, **P_1t_fft;

    P_current = alloc2float(p->nz, p->nx);
    P_past = alloc2float(p->nz, p->nx);
    P_next = alloc2float(p->nz, p->nx);
    P_1t = alloc2float(p->nz, p->nx);

    // 置0初始化
    memset(P_current[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_past[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_next[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_1t[0], 0, sizeof(float) * p->nx * p->nz);


    P_fft = alloc2secomplex(p->nz, p->nx);
    P_fft_temp = alloc2secomplex(p->nz, p->nx);
    P_1t_fft = alloc2secomplex(p->nz, p->nx);


    memset(P_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_fft_temp[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_1t_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);

    float **P_la = alloc2float(p->nz, p->nx);
    se_complex **P_la_fft = alloc2secomplex(p->nz, p->nx);

    float **amp = alloc2float(p->nz, p->nx);
    se_complex **amp_fft = alloc2secomplex(p->nz, p->nx);

    float **KP = alloc2float(p->nz, p->nx);
    se_complex **KP_fft = alloc2secomplex(p->nz, p->nx);


    // 置0初始化
    memset(P_la[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_la_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(amp[0], 0, sizeof(float) * p->nx * p->nz);
    memset(amp_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(KP[0], 0, sizeof(float) * p->nx * p->nz);
    memset(KP_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);


    // 计算波数k
    float *kx = (float *)malloc(p->nx * sizeof(float));
    float *ky = (float *)malloc(p->nz * sizeof(float));
    calculate_wavenumbers(kx, ky, p->nx, p->nz, p->dxz);


    // Butterworth低通滤波
    float **H = alloc2float(p->nz, p->nx);
    memset(H[0], 0, sizeof(float) * p->nx * p->nz);
    butterworth_lowpass_filter(H, kx, ky, p->nx, p->nz, p->vp, p->fc, p->order);


    for (it = 0; it < p->nt; it++)
    {
        if (flag == 0)
        {
            P_current[isx][isz] = p->init_record[isx][it]; // 添加雷克子波震源
        }

        else
        {
            for (ix = p->nbc; ix < p->nx - p->nbc; ix++)
            {

                P_current[ix][isz] =  p->init_record[ix][it]; // 添加检波器震源
            }
        }

        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {
                P_1t[ix][iz] = (P_current[ix][iz] - P_past[ix][iz]) / p->dt;
            }
        }

        // 将P_current放到P_fft的实部
        for (i = 0; i < p->nx; i++)
        {
            for (j = 0; j < p->nz; j++)
            {
                P_fft[i][j].r = P_current[i][j];
                P_fft[i][j].i = 0.0;
                P_1t_fft[i][j].r = P_1t[i][j];
                P_1t_fft[i][j].i = 0.0;
            }
        }

        // 进行二维FFT
        se_fftw_fft_2d_quick(P_fft, P_fft_temp, p->nx, p->nz, 0);
        memcpy(P_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
        se_fftw_fft_2d_quick(P_1t_fft, P_fft_temp, p->nx, p->nz, 0);
        memcpy(P_1t_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

        //伪谱或差分计算Laplace
        if (type_compute_Laplace == 0)
        {
            //有限差分八阶计算Laplace
            laplace_2d(P_current, P_la, p->nx, p->nz, p->dxz, 8);
        }
        else if (type_compute_Laplace == 1)
        {
            // 伪谱法计算Laplace
            for (ix = 0; ix < p->nx; ix++)
            {
                for (iz = 0; iz < p->nz; iz++)
                {
                    K2 = kx[ix] * kx[ix] + ky[iz] * ky[iz];
                    P_la_fft[ix][iz].r = -K2 * P_fft[ix][iz].r;
                    P_la_fft[ix][iz].i = -K2 * P_fft[ix][iz].i;

                }
            }
            se_fftw_fft_2d_quick(P_la_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(P_la_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

            for (i = 0; i < p->nx; i++)
            {
                for (j = 0; j < p->nz; j++)
                {
                    P_la[i][j] = P_la_fft[i][j].r;
                }
            }
        }
        else
        {
            err("error in type_compute_Laplace!");
        }



        if (strcmp(PS_or_CS, "PS") == 0)
        {
            // printf("Using PS extension!\n");
            // 伪谱法计算amp和KP
            for (ix = 0; ix < p->nx; ix++)
            {
                for (iz = 0; iz < p->nz; iz++)
                {
                    K2 = kx[ix] * kx[ix] + ky[iz] * ky[iz];

                    amp_fft[ix][iz].r = sqrt(K2) * P_1t_fft[ix][iz].r;
                    amp_fft[ix][iz].i = sqrt(K2) * P_1t_fft[ix][iz].i;

                    KP_fft[ix][iz].r = sqrt(K2) * P_fft[ix][iz].r;
                    KP_fft[ix][iz].i = sqrt(K2) * P_fft[ix][iz].i;

                    // 对PS延拓后的场滤波
                    amp_fft[ix][iz].r = amp_fft[ix][iz].r * H[ix][iz];
                    amp_fft[ix][iz].i = amp_fft[ix][iz].i * H[ix][iz];

                    KP_fft[ix][iz].r = KP_fft[ix][iz].r * H[ix][iz];
                    KP_fft[ix][iz].i = KP_fft[ix][iz].i * H[ix][iz];

                }
            }

            se_fftw_fft_2d_quick(KP_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(KP_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
            se_fftw_fft_2d_quick(amp_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(amp_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

            for (i = 0; i < p->nx; i++)
            {
                for (j = 0; j < p->nz; j++)
                {
                    amp[i][j] = amp_fft[i][j].r;
                    KP[i][j] = KP_fft[i][j].r;
                }
            }
        }
        else if (strcmp(PS_or_CS, "CS") == 0)
        {
            // printf("Using CS extension!\n");
            // CS延拓计算amp和KP
            caffarelli_extension_2d(P_1t,amp,p->nx,p->nz,p->dxz);
            caffarelli_extension_2d(P_current,KP,p->nx,p->nz,p->dxz);

            // 对CS延拓后的场滤波
            for (i = 0; i < p->nx; i++)
            {
                for (j = 0; j < p->nz; j++)
                {
                    amp_fft[i][j].r = amp[i][j];
                    amp_fft[i][j].i = 0.0;
                    KP_fft[i][j].r = KP[i][j];
                    KP_fft[i][j].i = 0.0;
                }
            }

            // 进行二维FFT
            se_fftw_fft_2d_quick(amp_fft, P_fft_temp, p->nx, p->nz, 0);
            memcpy(amp_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
            se_fftw_fft_2d_quick(KP_fft, P_fft_temp, p->nx, p->nz, 0);
            memcpy(KP_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

            // 低通滤波
            for (ix = 0; ix < p->nx; ix++)
            {
                for (iz = 0; iz < p->nz; iz++)
                {
                    amp_fft[ix][iz].r = amp_fft[ix][iz].r * H[ix][iz];
                    amp_fft[ix][iz].i = amp_fft[ix][iz].i * H[ix][iz];

                    KP_fft[ix][iz].r = KP_fft[ix][iz].r * H[ix][iz];
                    KP_fft[ix][iz].i = KP_fft[ix][iz].i * H[ix][iz];
                }
            }


            se_fftw_fft_2d_quick(KP_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(KP_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
            se_fftw_fft_2d_quick(amp_fft, P_fft_temp, p->nx, p->nz, 1);
            memcpy(amp_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);  

            for (i = 0; i < p->nx; i++)
            {
                for (j = 0; j < p->nz; j++)
                {
                    amp[i][j] = amp_fft[i][j].r;
                    KP[i][j] = KP_fft[i][j].r;
                }
            }  

        }


        if (it == 10)
        {
            // 打印p->type
            printf("p->type = %d\n", p->type);
        }

        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {

                if (p->type == 1)
                {
                    // 解耦的Aki模型波场传播，利用切比雪夫近似
                    P_next[ix][iz] = 2 * P_current[ix][iz] - P_past[ix][iz] + (pow(p->dt, 2) / C1[ix][iz]) * (C4[ix][iz] *amp_compensation_sign* -amp[ix][iz] + P_la[ix][iz] - C2[ix][iz] * KP[ix][iz] - C3[ix][iz] * P_current[ix][iz]);
                }
                else if (p->type == 0)
                {
                    // 声波方程
                    P_next[ix][iz] = 2 * P_current[ix][iz] - P_past[ix][iz] + pow(p->vp[ix][iz], 2) * pow(p->dt, 2) * P_la[ix][iz];
                }

                else
                {
                    err("error in type!");
                }
                
            }
        }

        // 吸收边界调用函数
        apply_absorbing_boundary_conditions(P_next, P_current, P_past, p);


        memcpy(P_past[0], P_current[0], sizeof(float) * p->nx * p->nz);
        memcpy(P_current[0], P_next[0], sizeof(float) * p->nx * p->nz);

        // 记录地震数据
        for (i = p->nbc; i < p->nx - p->nbc; i++)
        {
            p->record[i - p->nbc][it] = P_current[i][irz];
        }


        if (flag == 0)
        {
            for (ix = 0; ix < p->nx; ix++)
            {
                for (iz = 0; iz < p->nz; iz++)
                {
                    P_all_snapshots[ix][iz][it] = P_current[ix][iz];
                }
            }

        }
        else
        {
            // 成像
            for (ix = 0; ix < p->nx - 2*p->nbc; ix++)
            {
                for (iz = 0; iz < p->nz - 2*p->nbc; iz++)
                {
                    // // 互相关成像条件
                    Image[ix][iz] = Image[ix][iz] + P_current[ix + p->nbc][iz + p->nbc] * P_all_snapshots[ix + p->nbc][iz + p->nbc][p->nt - it - 1];

                    // 归一化互相关成像条件
                    // Image[ix][iz] = Image[ix][iz] + (P_current[ix + p->nbc][iz + p->nbc] * P_all_snapshots[ix + p->nbc][iz + p->nbc][p->nt - it - 1])/(P_all_snapshots[ix + p->nbc][iz + p->nbc][p->nt - it - 1]*P_all_snapshots[ix + p->nbc][iz + p->nbc][p->nt - it - 1]+1e-12);
                }
            }
        }


        if (it % 300 == 0)
        {
            SE_MESSAGE("Thread ID : %d :: Processing %d/%d\n", p->myid, it, p->nt);
        }


//        // 过程波场输出
//        int tfix;
//        tfix=500;
//        output_snapshot(it, tfix, p->nx, p->nz, P_current, "snap500.bin");
//
//        tfix=1000;
//        output_snapshot(it, tfix, p->nx, p->nz, P_current, "snap1000.bin");
//
//        tfix=1500;
//        output_snapshot(it, tfix, p->nx, p->nz, P_current, "snap1500.bin");

    }



    free(kx);
    free(ky);

    free2float(amp);
    free2secomplex(amp_fft);
    free2float(P_current);
    free2float(P_past);
    free2float(P_next);

    free2secomplex(P_fft);
    free2secomplex(P_fft_temp);

    free2secomplex(P_1t_fft);
    free2float(P_1t);

    free2float(C1);
    free2float(C2);
    free2float(C3);
    free2float(C4);

    free2float(KP);
    free2secomplex(KP_fft);
    free2float(P_la);
    free2secomplex(P_la_fft);

    free2float(H);

}


void visco_wave_2d_propagation_TVM_CS(visco_wave_2d_t *p, int flag, int flag_denoise, int amp_compensation_sign, float ***P_all_snapshots, float **Image)
{

    int i, j, k;
    int it, ix, iz;
    float K2;

    float w, s, t1, t2, t3;

    float c0, c1, c2;
    c0 = 411.7 + 4.36;
    c1 = -51.64;
    c2 = 4.366 / 2;

    int isx = (int)(p->sx);
    int isz = (int)(p->sz);
    int irz = (int)(p->rz);

    float **C1, **C2, **C3, **C4;
    C1 = alloc2float(p->nz, p->nx);
    C2 = alloc2float(p->nz, p->nx);
    C3 = alloc2float(p->nz, p->nx);
    C4 = alloc2float(p->nz, p->nx);

    for (i = 0; i < p->nx; i++)
    {
        for (j = 0; j < p->nz; j++)
        {
            C1[i][j] = (1 - 4 * c2 / (PI * p->Q[i][j]) + 2 * log(p->Omega0) / PI / p->Q[i][j])/ (p->vp[i][j] * p->vp[i][j]);
            C2[i][j] = 2 * c1 / (PI * p->Q[i][j]* p->vp[i][j]);
            C3[i][j] = 2 * (c0 - c2) / (PI * p->Q[i][j]* p->vp[i][j] * p->vp[i][j]);
            C4[i][j] = 1 / (p->vp[i][j] * p->Q[i][j]);
        }
    }



    float **P_current, **P_past, **P_next, **P_psm, **P_1t, **r_current, **r_past, **P_current_denoise, **P_current_Butterfilt;

    se_complex **P_fft_last, **P_fft, **P_fft_o, **P_fft_temp, **P_1t_fft, **p_fft_shift, **P_current_Butterfilt_fft;

    P_current = alloc2float(p->nz, p->nx);
    P_past = alloc2float(p->nz, p->nx);
    P_next = alloc2float(p->nz, p->nx);
    P_psm = alloc2float(p->nz, p->nx);
    P_1t = alloc2float(p->nz, p->nx);
    r_current = alloc2float(p->nz, p->nx);
    r_past = alloc2float(p->nz, p->nx);
    P_current_denoise = alloc2float(p->nz, p->nx);
    P_current_Butterfilt = alloc2float(p->nz, p->nx);

    // 置0初始化
    memset(P_current[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_past[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_next[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_psm[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_1t[0], 0, sizeof(float) * p->nx * p->nz);
    memset(r_current[0], 0, sizeof(float) * p->nx * p->nz);
    memset(r_past[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_current_denoise[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_current_Butterfilt[0], 0, sizeof(float) * p->nx * p->nz);

    P_fft_last = alloc2secomplex(p->nz, p->nx);
    P_fft = alloc2secomplex(p->nz, p->nx);
    P_fft_o = alloc2secomplex(p->nz, p->nx);
    P_fft_temp = alloc2secomplex(p->nz, p->nx);
    P_1t_fft = alloc2secomplex(p->nz, p->nx);
    P_current_Butterfilt_fft = alloc2secomplex(p->nz, p->nx);
    p_fft_shift = alloc2secomplex(p->nz, p->nx);

    memset(P_fft_last[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_fft_o[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_fft_temp[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_1t_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(P_current_Butterfilt_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(p_fft_shift[0], 0, sizeof(se_complex) * p->nx * p->nz);


    float **P_la = alloc2float(p->nz, p->nx);
    se_complex **P_la_fft = alloc2secomplex(p->nz, p->nx);

    float **amp = alloc2float(p->nz, p->nx);
    se_complex **amp_fft = alloc2secomplex(p->nz, p->nx);

    float **KP = alloc2float(p->nz, p->nx);
    se_complex **KP_fft = alloc2secomplex(p->nz, p->nx);

    // 置0初始化
    memset(P_la[0], 0, sizeof(float) * p->nx * p->nz);
    memset(P_la_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(amp[0], 0, sizeof(float) * p->nx * p->nz);
    memset(amp_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);
    memset(KP[0], 0, sizeof(float) * p->nx * p->nz);
    memset(KP_fft[0], 0, sizeof(se_complex) * p->nx * p->nz);

    //计算波数k
    float kmax = PI / p->dxz;
    float delta_kx = kmax / (p->nx / 2);
    float delta_ky = kmax / (p->nz / 2);
    float *kx = (float *)malloc(p->nx * sizeof(float));
    float *ky = (float *)malloc(p->nz * sizeof(float));

    for (i = 0; i < p->nx / 2; i++)
    {
        kx[i] = (i + 1) * delta_kx;
        kx[p->nx / 2 + i] = -kmax + (i + 1) * delta_kx;
    }
    for (i = 0; i < p->nz / 2; i++)
    {
        ky[i] = (i + 1) * delta_ky;
        ky[p->nz / 2 + i] = -kmax + (i + 1) * delta_ky;
    }
    
    // for (i = 0; i < p->nx; i++)
    // {
    //     printf("kx[%d] = %f\n", i, kx[i]);
    // }
    // exit(0);  
    

    // Butterworth低通滤波
    // 计算最大速度vmax
    float vmax=0;
    for (i = 0; i < p->nx; i++)
    {
        for (j = 0; j < p->nz; j++)
        {
            if (p->vp[i][j] > vmax)
            {
                vmax = p->vp[i][j];
            }
        }
    }
    float fc = 120; // 截止频率 
    float kc = fc * 2 * PI / vmax; // 截止波数

    int order = 2;      // 滤波器阶数
    float **k_matrix, **H;
    k_matrix = alloc2float(p->nz, p->nx);
    H = alloc2float(p->nz, p->nx);
    memset(k_matrix[0], 0, sizeof(float) * p->nx * p->nz);
    memset(H[0], 0, sizeof(float) * p->nx * p->nz);

    for (ix = 0; ix < p->nx; ix++)
    {
        for (iz = 0; iz < p->nz; iz++)
        {
            K2 = kx[ix] * kx[ix] + ky[iz] * ky[iz];

            k_matrix[ix][iz] = sqrt(K2);                                       // 计算波数矩阵
            H[ix][iz] = 1.0 / (1.0 + pow((k_matrix[ix][iz] / kc), 2 * order)); // 构造 Butterworth低通滤波器（理想圆形滤波器）
        }

        // exit(0);
    }

    // write_binary_file_2d("k_matrix.bin", k_matrix, p->nx, p->nz); // 输出k_matrix到二进制文件

    // 将H进行fftshift
    // fftshift_2d(H, p->nx, p->nz);


    // write_binary_file_2d("H_new.bin", H, p->nx, p->nz); // 输出滤波器到二进制文件
    // exit(0);

    for (it = 0; it < p->nt; it++)
    {
        if (flag == 0)
        {
            P_current[isx][isz] = p->init_record[isx][it]; // 添加雷克子波震源
        }

        else
        {
            for (ix = p->nbc; ix < p->nx - p->nbc; ix++)
            {

                P_current[ix][isz] =  p->init_record[ix][it]; // 添加检波器震源
            }
        }

        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {
                P_1t[ix][iz] = (P_current[ix][iz] - P_past[ix][iz]) / p->dt;
            }
        }

        // 将P_current放到P_fft的实部
        for (i = 0; i < p->nx; i++)
        {
            for (j = 0; j < p->nz; j++)
            {
                P_fft[i][j].r = P_current[i][j];
                P_fft[i][j].i = 0.0;
                P_1t_fft[i][j].r = P_1t[i][j];
                P_1t_fft[i][j].i = 0.0;
            }
        }

        // 进行二维FFT
        se_fftw_fft_2d_quick(P_fft, P_fft_o, p->nx, p->nz, 0);
        memcpy(P_fft[0], P_fft_o[0], sizeof(se_complex) * p->nx * p->nz);
        // se_fftw_fft_2d_quick(P_1t_fft, P_fft_temp, p->nx, p->nz, 0);
        // memcpy(P_1t_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {
                K2 = kx[ix] * kx[ix] + ky[iz] * ky[iz];
                P_la_fft[ix][iz].r = -K2 * P_fft[ix][iz].r;
                P_la_fft[ix][iz].i = -K2 * P_fft[ix][iz].i;

                // amp_fft[ix][iz].r = sqrt(K2) * P_1t_fft[ix][iz].r;
                // amp_fft[ix][iz].i = sqrt(K2) * P_1t_fft[ix][iz].i;

                // KP_fft[ix][iz].r = sqrt(K2) * P_fft[ix][iz].r;
                // KP_fft[ix][iz].i = sqrt(K2) * P_fft[ix][iz].i;
            }
        }

        
        // // Butterworth低通滤波
        // // 对P_fft低通滤波
        // for (ix = 0; ix < p->nx; ix++)
        // {
        //     for (iz = 0; iz < p->nz; iz++)
        //     {
        //         amp_fft[ix][iz].r = amp_fft[ix][iz].r * H[ix][iz];
        //         amp_fft[ix][iz].i = amp_fft[ix][iz].i * H[ix][iz];

        //         KP_fft[ix][iz].r = KP_fft[ix][iz].r * H[ix][iz];
        //         KP_fft[ix][iz].i = KP_fft[ix][iz].i * H[ix][iz];
        //     }
        // }


        se_fftw_fft_2d_quick(P_la_fft, P_fft_temp, p->nx, p->nz, 1);
        memcpy(P_la_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
        // se_fftw_fft_2d_quick(KP_fft, P_fft_temp, p->nx, p->nz, 1);
        // memcpy(KP_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
        // se_fftw_fft_2d_quick(amp_fft, P_fft_temp, p->nx, p->nz, 1);
        // memcpy(amp_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

        for (i = 0; i < p->nx; i++)
        {
            for (j = 0; j < p->nz; j++)
            {
                P_la[i][j] = P_la_fft[i][j].r;
                // amp[i][j] = amp_fft[i][j].r;
                // KP[i][j] = KP_fft[i][j].r;
            }
        }

        // laplace_2d(P_current,P_la,p->nx,p->nz,p->dxz,6);

        caffarelli_extension_2d(P_1t,amp,p->nx,p->nz,p->dxz);
        caffarelli_extension_2d(P_current,KP,p->nx,p->nz,p->dxz);

        // 对CS延拓后的场滤波
        for (i = 0; i < p->nx; i++)
        {
            for (j = 0; j < p->nz; j++)
            {
                amp_fft[i][j].r = amp[i][j];
                amp_fft[i][j].i = 0.0;
                KP_fft[i][j].r = KP[i][j];
                KP_fft[i][j].i = 0.0;
            }
        }

        // 进行二维FFT
        se_fftw_fft_2d_quick(amp_fft, P_fft_temp, p->nx, p->nz, 0);
        memcpy(amp_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
        se_fftw_fft_2d_quick(KP_fft, P_fft_temp, p->nx, p->nz, 0);
        memcpy(KP_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);

        // 低通滤波
        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {
                amp_fft[ix][iz].r = amp_fft[ix][iz].r * H[ix][iz];
                amp_fft[ix][iz].i = amp_fft[ix][iz].i * H[ix][iz];

                KP_fft[ix][iz].r = KP_fft[ix][iz].r * H[ix][iz];
                KP_fft[ix][iz].i = KP_fft[ix][iz].i * H[ix][iz];
            }
        }


        se_fftw_fft_2d_quick(KP_fft, P_fft_temp, p->nx, p->nz, 1);
        memcpy(KP_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);
        se_fftw_fft_2d_quick(amp_fft, P_fft_temp, p->nx, p->nz, 1);
        memcpy(amp_fft[0], P_fft_temp[0], sizeof(se_complex) * p->nx * p->nz);  

        for (i = 0; i < p->nx; i++)
        {
            for (j = 0; j < p->nz; j++)
            {
                // P_la[i][j] = P_la_fft[i][j].r;
                amp[i][j] = amp_fft[i][j].r;
                KP[i][j] = KP_fft[i][j].r;
            }
        }  

        
        if (it == 10)
        {
            // 打印p->type
            printf("p->type = %d\n", p->type);
        }

        for (ix = 0; ix < p->nx; ix++)
        {
            for (iz = 0; iz < p->nz; iz++)
            {


                if (p->type == 1)
                {
                    // 解耦的Aki模型波场传播，利用切比雪夫近似
                    P_next[ix][iz] = 2 * P_current[ix][iz] - P_past[ix][iz] + (pow(p->dt, 2) / C1[ix][iz]) * (C4[ix][iz] *amp_compensation_sign* -amp[ix][iz] + P_la[ix][iz] - C2[ix][iz] * KP[ix][iz] - C3[ix][iz] * P_current[ix][iz]);
                }
                else if (p->type == 0)
                {
                    // 声波方程
                    P_next[ix][iz] = 2 * P_current[ix][iz] - P_past[ix][iz] + pow(p->vp[ix][iz], 2) * pow(p->dt, 2) * P_la[ix][iz];
                }

                else
                {
                    err("error in type!");
                }

                // 吸收边界
                // left
                if (ix >= 0 && ix < p->L - 4 && iz >= 0 && iz < p->nz)
                {
                    w = 1 - 1.0 * ix / p->L;
                    s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                    t1 = (2 - s) * (1 - s) / 2;
                    t2 = s * (2 - s);
                    t3 = s * (s - 1) / 2;

                    P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                     t2 * P_current[ix + 1][iz] + t3 * P_current[ix + 2][iz]) +
                                          (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                      2 * t1 * t2 * P_past[ix + 1][iz] +
                                                      (2 * t1 * t3 + t2 * t2) * P_past[ix + 2][iz] +
                                                      2 * t2 * t3 * P_past[ix + 3][iz] +
                                                      t3 * t3 * P_past[ix + 4][iz])) +
                                     (1 - w) * P_next[ix][iz];
                }

                // right
                if (ix > p->nx - p->L + 3 && ix < p->nx && iz >= 0 && iz < p->nz)
                {
                    w = 1 - 1.0 * (p->nx - ix - 1) / p->L;
                    s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                    t1 = (2 - s) * (1 - s) / 2;
                    t2 = s * (2 - s);
                    t3 = s * (s - 1) / 2;

                    P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                     t2 * P_current[ix - 1][iz] + t3 * P_current[ix - 2][iz]) +
                                          (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                      2 * t1 * t2 * P_past[ix - 1][iz] +
                                                      (2 * t1 * t3 + t2 * t2) * P_past[ix - 2][iz] +
                                                      2 * t2 * t3 * P_past[ix - 3][iz] +
                                                      t3 * t3 * P_past[ix - 4][iz])) +
                                     (1 - w) * P_next[ix][iz];
                }

                // top
                if (iz >= 0 && iz < p->L - 4 && ix >= 0 && ix < p->nx)
                {
                    w = 1 - 1.0 * iz / p->L;
                    s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                    t1 = (2 - s) * (1 - s) / 2;
                    t2 = s * (2 - s);
                    t3 = s * (s - 1) / 2;

                    P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                     t2 * P_current[ix][iz + 1] + t3 * P_current[ix][iz + 2]) +
                                          (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                      2 * t1 * t2 * P_past[ix][iz + 1] +
                                                      (2 * t1 * t3 + t2 * t2) * P_past[ix][iz + 2] +
                                                      2 * t2 * t3 * P_past[ix][iz + 3] +
                                                      t3 * t3 * P_past[ix][iz + 4])) +
                                     (1 - w) * P_next[ix][iz];
                }

                // bottom
                if (iz > p->nz - p->L + 3 && iz < p->nz && ix >= 0 && ix < p->nx)
                {
                    w = 1 - 1.0 * (p->nz - iz - 1) / p->L;
                    s = p->alpha * p->vp[ix][iz] * p->dt / p->dxz;
                    t1 = (2 - s) * (1 - s) / 2;
                    t2 = s * (2 - s);
                    t3 = s * (s - 1) / 2;

                    P_next[ix][iz] = w * ((1 * 2) * (t1 * P_current[ix][iz] +
                                                     t2 * P_current[ix][iz - 1] + t3 * P_current[ix][iz - 2]) +
                                          (-1 * 1) * (t1 * t1 * P_past[ix][iz] +
                                                      2 * t1 * t2 * P_past[ix][iz - 1] +
                                                      (2 * t1 * t3 + t2 * t2) * P_past[ix][iz - 2] +
                                                      2 * t2 * t3 * P_past[ix][iz - 3] +
                                                      t3 * t3 * P_past[ix][iz - 4])) +
                                     (1 - w) * P_next[ix][iz];
                }
            }
        }

        memcpy(P_past[0], P_current[0], sizeof(float) * p->nx * p->nz);
        memcpy(P_current[0], P_next[0], sizeof(float) * p->nx * p->nz);

        memcpy(r_past[0], r_current[0], sizeof(float) * p->nx * p->nz);

        // 记录地震数据
        for (i = p->nbc; i < p->nx - p->nbc; i++)
        {
            p->record[i - p->nbc][it] = P_current[i][irz];
        }

        // if (it==550)
        // {
        //     write_binary_file_2d("../result/homo/snap550_DSLS_S_Nodenoise.bin", P_current, p->nx, p->nz);
        // }

        // if (flag_denoise == 1)
        // {
        //     // laplace滤波
        //     laplace_denoise(P_current, P_current_denoise, p->nx, p->nz);
        //     for (ix = 0; ix < p->nx; ix++)
        //     {
        //         for (iz = 0; iz < p->nz; iz++)
        //         {
        //             P_all_snapshots[ix][iz][it] = P_current_denoise[ix][iz];
        //         }
        //     }
        // }
        // else
        // {
        //     for (ix = 0; ix < p->nx; ix++)
        //     {
        //         for (iz = 0; iz < p->nz; iz++)
        //         {
        //             P_all_snapshots[ix][iz][it] = P_current[ix][iz];
        //         }
        //     }
        // }


        if (flag == 0)
        {
            for (ix = 0; ix < p->nx; ix++)
            {
                for (iz = 0; iz < p->nz; iz++)
                {
                    P_all_snapshots[ix][iz][it] = P_current[ix][iz];
                }
            }

        }
        else
        {
            // 成像
            for (ix = 0; ix < p->nx - 2*p->nbc; ix++)
            {
                for (iz = 0; iz < p->nz - 2*p->nbc; iz++)
                {
                    Image[ix][iz] = Image[ix][iz] + P_current[ix + p->nbc][iz + p->nbc] * P_all_snapshots[ix + p->nbc][iz + p->nbc][p->nt - it - 1];
                }
            }
        }


        if (it % 100 == 0)
        {
            SE_MESSAGE("Thread ID : %d :: Processing %d/%d\n", p->myid, it, p->nt);
        }


//        // 过程波场输出
//        if (it == 250)
//        {
//            printf("Snap Output!! nx=%d,nz=%d\n", p->nx, p->nz);
//            FILE *fp = fopen("snap250_CS.bin", "wb");
//            if (fp == NULL)
//            {
//                err("error in opening file snap");
//            }
//            
//
//            for (i = 0; i < p->nx; i++)
//            {
//                for (j = 0; j < p->nz; j++)
//                {
//                    fwrite(&P_current[i][j], sizeof(float), 1, fp);
//                }
//            }
//            fclose(fp);
//        }
//
//        if (it == 500)
//        {
//            printf("Snap Output!! nx=%d,nz=%d\n", p->nx, p->nz);
//            FILE *fp = fopen("snap500_CS.bin", "wb");
//            if (fp == NULL)
//            {
//                err("error in opening file snap");
//            }
//            
//
//            for (i = 0; i < p->nx; i++)
//            {
//                for (j = 0; j < p->nz; j++)
//                {
//                    fwrite(&P_current[i][j], sizeof(float), 1, fp);
//                }
//            }
//            fclose(fp);
//        }
//
//        if (it == 1000)
//        {
//            printf("Snap Output!! nx=%d,nz=%d\n", p->nx, p->nz);
//            FILE *fp = fopen("snap1000_CS.bin", "wb");
//            if (fp == NULL)
//            {
//                err("error in opening file snap");
//            }
//            
//
//            for (i = 0; i < p->nx; i++)
//            {
//                for (j = 0; j < p->nz; j++)
//                {
//                    fwrite(&P_current[i][j], sizeof(float), 1, fp);
//                }
//            }
//            fclose(fp);
//        }
//
//        if (it == 1500)
//        {
//            printf("Snap Output!! nx=%d,nz=%d\n", p->nx, p->nz);
//            FILE *fp = fopen("snap1500_CS.bin", "wb");
//            if (fp == NULL)
//            {
//                err("error in opening file snap");
//            }
//            
//
//            for (i = 0; i < p->nx; i++)
//            {
//                for (j = 0; j < p->nz; j++)
//                {
//                    fwrite(&P_current[i][j], sizeof(float), 1, fp);
//                }
//            }
//            fclose(fp);
//        }
//
//
//        if (it == 2450)
//        {
//            printf("Snap Output!! nx=%d,nz=%d\n", p->nx, p->nz);
//            FILE *fp = fopen("snap2450_CS.bin", "wb");
//            if (fp == NULL)
//            {
//                err("error in opening file snap");
//            }
//
//            for (i = 0; i < p->nx; i++)
//            {
//                for (j = 0; j < p->nz; j++)
//                {
//                    fwrite(&P_current[i][j], sizeof(float), 1, fp);
//                }
//            }
//            fclose(fp);
//        }

        

    }

    free(kx);
    free(ky);


    free2float(amp);
    free2secomplex(amp_fft);
    free2float(P_current);
    free2float(P_past);
    free2float(P_next);
    free2float(P_psm);
    free2float(P_current_denoise);

    free2float(r_current);
    free2float(r_past);

    free2secomplex(P_fft_last);
    free2secomplex(P_fft);
    free2secomplex(P_fft_o);
    free2secomplex(P_fft_temp);

    free2secomplex(P_1t_fft);
    free2float(P_1t);

    free2float(C1);
    free2float(C2);
    free2float(C3);
    free2float(C4);

    free2float(KP);
    free2secomplex(KP_fft);
    free2float(P_la);
    free2secomplex(P_la_fft);


    free2float(P_current_Butterfilt);
    free2secomplex(P_current_Butterfilt_fft);
    free2float(H);
    free2float(k_matrix);

}



