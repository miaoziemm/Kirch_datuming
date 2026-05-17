# 将时间域 Kirchhoff 平面外推代码重构为“蝶形算法友好”的频率域形式：实现指导文件

本文档用于指导另一个 AI 或程序员，将当前的时间域波场外推代码重构为适合后续直接替换为蝶形算法的频率域代码。目标不是马上实现蝶形算法，而是先把现有代码改造成标准的频率域空间积分算子形式：

\[
U_{\mathrm{out}}(r,\omega)=\sum_s K(r,s,\omega)U_{\mathrm{in}}(s,\omega),
\]

然后将其中的空间求和接口替换为蝶形算法。重构后的代码必须保持现有参数、输入、输出文件格式不变，并且数值结果应能与原时间域结果对上。

---

## 1. 必须保持不变的内容

修改代码时，不得改变外部使用方式。以下内容必须保持一致：

1. 命令行参数名称与含义不变：
   - `input_file`
   - `output_file`
   - `sgreen_file`
   - `rgreen_file`
   - `model_file`
   - `verb`
   - `cmp`
   - `aperture`
   - `taper`
   - `length`
   - `interm`

2. 输入数据维度不变：

\[
\text{input}(t,h,s), \qquad n_1=nt,\quad n_2=nh,\quad n_3=ns.
\]

3. 输出数据维度不变：

\[
\text{output}(t,h,s), \qquad n_1=nt,\quad n_2=nh,\quad n_3=ns.
\]

4. 输出 RSF/SEP 头信息必须保持与原代码一致：

```c
out->headers->ndim = 3;
out->headers->n[0] = nt;
out->headers->n[1] = nh;
out->headers->n[2] = ns;
out->headers->d[0] = dt;
out->headers->d[1] = dh;
out->headers->d[2] = ds;
out->headers->o[0] = 0.0f;
out->headers->o[1] = h0;
out->headers->o[2] = s0;
```

5. 原代码中的两次时间反转逻辑必须保持：

```c
reverse_trace(nt, nh, ns, tr_in);
...
reverse_trace(nt, nh, ns, tr_in);
```

也就是说，频率域实现应在第一次 `reverse_trace` 之后进行频率域外推，并在最终写出前执行第二次 `reverse_trace`，以保证输出和原时间域代码的时间方向一致。

6. `interm` 中间结果的行为也要保持一致。原代码在检波点端外推完成后直接写出 `tr_out`，没有对 `tr_out` 再做时间反转。因此频率域版本如果指定了 `interm`，也应在检波点端外推完成并 IFFT 回时间域后写出同样方向的 `tr_out`。

---

## 2. 当前时间域代码的实际计算结构

当前程序分两步进行平面到平面的外推。

### 2.1 检波点端外推

固定震源编号 `is`，对检波点方向进行外推：

```c
for (is = 0; is < ns; is++) {
    for (ih = 0; ih < nh; ih++) {
        for (ic = left; ic <= right; ic++) {
            tau = rtable[cc][c];
            dist = rdatum*rdatum + (ic-ih)*dh*(ic-ih)*dh;
            tr_out[is][ih][it] += coef/M_PI
                * dh*rdatum*tau/dist
                * kirdat_pick(delta, tr_in[is][ic], shift);
        }
    }
}
```

对应的单频空间积分形式应为：

\[
U_r(i_s,i_h,\omega)
=
\sum_{i_c\in\mathcal A_r(i_h)}
K_r(i_h,i_c,i_s,\omega)
U_{\mathrm{in}}(i_s,i_c,\omega).
\]

其中 \(i_h\) 是目标检波点索引，\(i_c\) 是原平面上的输入检波点索引。

### 2.2 震源端外推

检波点端外推完成后，程序再对震源端进行外推。

#### `cmp=0`：检波点坐标为绝对坐标

```c
for (ih = 0; ih < nh; ih++) {
    for (is = 0; is < ns; is++) {
        for (ic = left; ic <= right; ic++) {
            tau = stable[cc][c];
            dist = sdatum*sdatum + (ic-is)*ds*(ic-is)*ds;
            tr_in[is][ih][it] += coef/M_PI
                * ds*sdatum*tau/dist
                * kirdat_pick(delta, tr_out[ic][ih], shift);
        }
    }
}
```

对应频率域形式为：

\[
U_{\mathrm{out}}(i_s,i_h,\omega)
=
\sum_{i_c\in\mathcal A_s(i_s)}
K_s(i_s,i_c,\omega)
U_r(i_c,i_h,\omega).
\]

#### `cmp=1`：检波点坐标为相对震源的 offset

原代码不是简单的固定 `ih` 做源端求和，而是先按绝对检波点坐标 `r` 重组共检波点道集：

```c
r = ir*dr + ...;
ih = (r - is*ds)/dh + 0.5f;
hh = (r - ic*ds)/dh + 0.5f;
tr_in[is][ih][it] += ... * kirdat_pick(delta, tr_out[ic][hh], shift);
```

对应频率域形式为：

\[
U_{\mathrm{out}}(i_s,i_h(i_s,r),\omega)
=
\sum_{i_c\in\mathcal A_s(i_s,r)}
K_s(i_s,i_c,r,\omega)
U_r(i_c,h(i_c,r),\omega).
\]

因此 `cmp=1` 分支后续改蝶形算法时不能简单按二维规则数组处理，需要先保留原来的 `ir`、`is`、`ic`、`ih`、`hh` 映射逻辑。

---

## 3. 原时间域滤波器的严格离散形式

原代码中的 `filt_set(tau)` 构造了一个依赖走时 \(\tau\) 的二阶差分滤波核。

令时间采样间隔为：

\[
\Delta t = dt.
\]

滤波长度对应采样点数：

\[
N = \left\lfloor \frac{\texttt{length}}{\Delta t} \right\rfloor +2.
\]

先定义：

\[
G_k(\tau)=
\sqrt{
\left(
\frac{\tau+k\Delta t}{\tau}
\right)^2-1
},
\qquad k=0,1,\cdots,N-1.
\]

因为 \(G_0=0\)，代码中先令 `filt[0]=0`。

最终滤波核为：

\[
F_0(\tau)=G_1(\tau)-G_0(\tau),
\]

\[
F_k(\tau)=G_{k+1}(\tau)-2G_k(\tau)+G_{k-1}(\tau),
\qquad k=1,2,\cdots,N-2.
\]

这个 \(F_k\) 对应代码中最终的 `filt[its][k]`。

---

## 4. 与时间域代码完全等价的频率域单对道核

这是重构时最重要的部分。不要把时间域算子简单写成 \(\exp(i\omega\tau)\)。原代码使用的是：

1. 走时延迟；
2. 有限长度二阶差分滤波核；
3. 由 `delta` 决定的线性插值；
4. 截断的因果卷积。

### 4.1 时间域代码中的 `shift` 和 `delta`

原代码对每一个 \(\tau\) 使用：

```c
shift = 0;
delta = 0.0f;
for (it = 0; it < nt; it++) {
    if (((float)it)*dt < tau)
        continue;
    else if (shift == 0)
        delta = (((float)it*dt)-tau)/dt;

    out[it] += weight * kirdat_pick(delta, trace, shift);
    shift++;
}
```

因此第一个参与计算的时间采样点为：

\[
n_0 = \left\lceil \frac{\tau}{\Delta t} \right\rceil.
\]

线性插值参数为：

\[
\delta = \frac{n_0\Delta t - \tau}{\Delta t}
= n_0 - \frac{\tau}{\Delta t}.
\]

注意这里不是 \(\tau/\Delta t-\lfloor\tau/\Delta t\rfloor\)，而是 `ceil` 之后的剩余量。这个符号方向必须与原代码保持一致。

### 4.2 时间域等价卷积核

对一对输入道 \(x[n]\) 和输出道 \(y[n]\)，原代码等价于：

\[
y[n]
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}
F_k(\tau)
\left[
(1-\delta)x[n-n_0-k]
+
\delta x[n-n_0-k+1]
\right],
\]

其中越界索引按零处理。空间权重尚未包含在上式中。

因此，单对道的严格离散频率响应为：

\[
\boxed{
H(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}
F_k(\tau)
\left[
(1-\delta)e^{-i\omega(n_0+k)\Delta t}
+
\delta e^{-i\omega(n_0+k-1)\Delta t}
\right]
}
\]

这是在 FFTW 默认正变换约定下应使用的形式。

FFTW 的正变换为：

\[
X_m=\sum_{n=0}^{N_{\mathrm{fft}}-1}x_n e^{-i2\pi mn/N_{\mathrm{fft}}},
\]

反变换为：

\[
x_n=\frac{1}{N_{\mathrm{fft}}}\sum_{m=0}^{N_{\mathrm{fft}}-1}X_m e^{i2\pi mn/N_{\mathrm{fft}}}.
\]

在该约定下，延迟 \(x[n-p]\) 对应频率域乘子 \(e^{-i\omega p\Delta t}\)。

### 4.3 写成 \(\exp(i\omega t(s,r))\) 的形式

如果论文或蝶形算法推导希望写成：

\[
K(r,s,\omega)=A(r,s,\omega)\exp[i\omega t(s,r)],
\]

可以这样写，但在实际 FFTW 代码中要注意符号约定。对于 FFTW 正变换，代码中直接使用的核更自然地写成：

\[
K(r,s,\omega)=\widetilde A(r,s,\omega)\exp[-i\omega t(s,r)].
\]

若必须写成 \(\exp[i\omega t]\)，则可以将相位符号并入频率定义或使用共轭核。实现时建议以 FFTW 约定为准，直接使用上面的 \(H(\omega,\tau)\) 公式，避免符号错误。

---

## 5. 完整的频率域空间外推核

### 5.1 通用单端外推公式

设：

- \(s\)：原平面上的输入点；
- \(r\)：外推后平面上的目标点；
- \(z_d\)：两个平面之间的外推距离；
- \(\Delta s\)：空间采样间隔；
- \(t(s,r)=\tau\)：由 Green 函数表给出的走时；
- \(w(r,s)\)：孔径 taper 系数。

则频率域外推为：

\[
\boxed{
U_{\mathrm{out}}(r,\omega)
=
\sum_{s\in\mathcal A(r)}
K(r,s,\omega)U_{\mathrm{in}}(s,\omega)
}
\]

其中：

\[
\boxed{
K(r,s,\omega)=
w(r,s)
\frac{\Delta s}{\pi}
\frac{z_d\,t(s,r)}{z_d^2+|s-r|^2}
H\left(\omega,t(s,r)\right)
}
\]

这里的 \(H(\omega,\tau)\) 必须使用第 4 节中的完整离散频率响应。

### 5.2 检波点端核

对检波点端外推：

\[
\boxed{
K_r(i_h,i_c,i_s,\omega)=
q_r(i_h,i_c)
\frac{\Delta h}{\pi}
\frac{z_r\,\tau_r(i_c,i_h;i_s)}
{z_r^2+(i_c-i_h)^2\Delta h^2}
H\left(\omega,\tau_r(i_c,i_h;i_s)\right)
}
\]

其中：

\[
z_r=\texttt{rdatum},
\qquad
\Delta h=\texttt{dh}.
\]

当 `cmp=1` 时：

\[
\tau_r(i_c,i_h;i_s)=
\texttt{rtable[cc][c]},
\]

\[
c=\mathrm{round}\left(\frac{s_0+i_s\Delta s+h_0+i_h\Delta h-rg_0}{drg}\right),
\]

\[
cc=\mathrm{round}\left(\frac{s_0+i_s\Delta s+h_0+i_c\Delta h-rg_0}{drg}\right).
\]

当 `cmp=0` 时：

\[
c=\mathrm{round}\left(\frac{h_0+i_h\Delta h-rg_0}{drg}\right),
\]

\[
cc=\mathrm{round}\left(\frac{h_0+i_c\Delta h-rg_0}{drg}\right).
\]

这里的 `round` 必须模仿 C 代码中的：

```c
(int)(x + 0.5f)
```

不要随意换成银行家舍入。

### 5.3 震源端核，`cmp=0`

对于绝对检波点坐标，震源端外推为：

\[
\boxed{
K_s(i_s,i_c,\omega)=
q_s(i_s,i_c)
\frac{\Delta s}{\pi}
\frac{z_s\,\tau_s(i_c,i_s)}
{z_s^2+(i_c-i_s)^2\Delta s^2}
H\left(\omega,\tau_s(i_c,i_s)\right)
}
\]

其中：

\[
z_s=\texttt{sdatum},
\qquad
\Delta s=\texttt{ds},
\]

\[
\tau_s(i_c,i_s)=\texttt{stable[cc][c]},
\]

\[
c=\mathrm{round}\left(\frac{s_0+i_s\Delta s-sg_0}{dsg}\right),
\]

\[
cc=\mathrm{round}\left(\frac{s_0+i_c\Delta s-sg_0}{dsg}\right).
\]

### 5.4 震源端核，`cmp=1`

对于相对 offset 坐标，应保持原来的共检波点循环。对每一个绝对检波点坐标 \(r\)，有：

\[
ih(i_s,r)=\mathrm{round}\left(\frac{r-i_s\Delta s}{\Delta h}\right),
\]

\[
hh(i_c,r)=\mathrm{round}\left(\frac{r-i_c\Delta s}{\Delta h}\right).
\]

频率域形式为：

\[
\boxed{
U_{\mathrm{out}}(i_s,ih(i_s,r),\omega)
=
\sum_{i_c\in\mathcal A_s(i_s,r)}
K_s(i_s,i_c,r,\omega)
U_r(i_c,hh(i_c,r),\omega)
}
\]

其中：

\[
\boxed{
K_s(i_s,i_c,r,\omega)=
q_s(i_s,i_c,r)
\frac{\Delta s}{\pi}
\frac{z_s\,\tau_s(i_c,i_s)}
{z_s^2+(i_c-i_s)^2\Delta s^2}
H\left(\omega,\tau_s(i_c,i_s)\right)
}
\]

这部分暂时不建议直接做蝶形算法，应先保留原来的 `ir` gather 结构，只把时间滤波换成频率域核乘法。等 `cmp=0` 的规则轴版本验证通过后，再处理 `cmp=1` 的非规则映射。

---

## 6. 推荐的新代码结构

当前代码是“每个输入点–输出点对”在时间域逐点求和。重构后应变成下面的结构。

### 6.1 主流程

```text
1. 读取参数，保持原代码完全一致。
2. 读取 input_file、model_file、sgreen_file、rgreen_file，保持原代码完全一致。
3. 分配 tr_in[ns][nh][nt]、tr_out[ns][nh][nt]。
4. 读取输入数据到 tr_in。
5. reverse_trace(nt, nh, ns, tr_in)。
6. 选择 nfft，要求能避免循环卷积污染。
7. 对所有输入道做 FFT：
      tr_in[is][ih][t] -> U_in[is][ih][iw]
8. 频率域检波点端外推：
      U_tmp = ApplyReceiverOperator(U_in)
9. 若指定 interm：
      IFFT U_tmp -> tr_out
      写出 tr_out，保持原代码中 interm 的时间方向。
10. 频率域震源端外推：
      U_out = ApplySourceOperator(U_tmp)
11. IFFT U_out -> tr_in
12. reverse_trace(nt, nh, ns, tr_in)
13. 写 output_file。
14. 释放内存。
```

### 6.2 需要新增的核心函数

建议将代码拆成以下函数，方便后续直接替换为空间蝶形算法：

```c
typedef struct {
    int nt, nh, ns;
    int nfft, nw;
    float dt, dh, ds;
    float h0, s0;
    float rdatum, sdatum;
    float length;
    int nsam;
    int aper, tap;
    int cmp;
    int nsg, nrg;
    float sg0, dsg, rg0, drg;
    float **stable;
    float **rtable;
} ExtrapCtx;
```

```c
int choose_nfft(const ExtrapCtx *ctx);
```

```c
void forward_fft_all_traces(
    const ExtrapCtx *ctx,
    float ***tr,
    fftwf_complex ***U
);
```

```c
void inverse_fft_all_traces(
    const ExtrapCtx *ctx,
    fftwf_complex ***U,
    float ***tr
);
```

```c
void compute_filter_F(
    float tau,
    float dt,
    int nsam,
    float *F
);
```

```c
fftwf_complex compute_H_tau(
    float omega,
    float tau,
    float dt,
    int nsam,
    const float *F
);
```

```c
fftwf_complex receiver_kernel(
    const ExtrapCtx *ctx,
    int is,
    int ih,
    int ic,
    int iw
);
```

```c
fftwf_complex source_kernel_cmp0(
    const ExtrapCtx *ctx,
    int is,
    int ic,
    int iw
);
```

```c
fftwf_complex source_kernel_cmp1(
    const ExtrapCtx *ctx,
    int is,
    int ic,
    int ir,
    int iw
);
```

```c
void apply_receiver_direct(
    const ExtrapCtx *ctx,
    fftwf_complex ***U_in,
    fftwf_complex ***U_tmp
);
```

```c
void apply_source_direct_cmp0(
    const ExtrapCtx *ctx,
    fftwf_complex ***U_tmp,
    fftwf_complex ***U_out
);
```

```c
void apply_source_direct_cmp1(
    const ExtrapCtx *ctx,
    fftwf_complex ***U_tmp,
    fftwf_complex ***U_out
);
```

后续实现蝶形算法时，只替换：

```c
apply_receiver_direct(...)
apply_source_direct_cmp0(...)
apply_source_direct_cmp1(...)
```

其中 `compute_H_tau()` 和 `receiver_kernel()` / `source_kernel()` 不应被破坏，因为它们保证与时间域公式一致。

---

## 7. FFT 长度与线性卷积

原时间域代码是因果截断求和，不是循环卷积。因此频率域实现不能直接使用长度为 `nt` 的 FFT 做循环卷积，否则高延迟部分会从尾部绕回到前面。

建议：

\[
N_{\mathrm{fft}} \ge nt + n_{\max},
\]

其中：

\[
n_{\max}=\left\lceil \frac{\tau_{\max}}{\Delta t} \right\rceil + N - 2.
\]

\(\tau_{\max}\) 应从 `stable` 和 `rtable` 中实际可能使用的走时里取最大值。为了简单稳妥，可以先扫描两个表得到全局最大值。

然后将 `nfft` 取为不小于该值的高效 FFT 长度，例如 2 的幂：

```c
int nfft = 1;
while (nfft < nt + max_delay) nfft <<= 1;
```

使用 real-to-complex FFT 时：

\[
nw=\frac{nfft}{2}+1.
\]

IFFT 后只保留前 `nt` 个采样点，后面的补零部分丢弃。

---

## 8. `compute_H_tau()` 的实现要点

### 8.1 公式

对 FFTW 正变换约定，应实现：

\[
H(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}
F_k(\tau)
\left[
(1-\delta)e^{-i\omega(n_0+k)\Delta t}
+
\delta e^{-i\omega(n_0+k-1)\Delta t}
\right].
\]

其中：

\[
n_0=\left\lceil \frac{\tau}{\Delta t}\right\rceil,
\qquad
\delta=n_0-\frac{\tau}{\Delta t}.
\]

### 8.2 伪代码

```c
static inline complexf compute_H_tau(float omega, float tau,
                                     float dt, int nsam)
{
    int n0 = (int)ceilf(tau / dt);
    float delta = n0 - tau / dt;

    // 为了避免 tau 恰好落在采样点时出现 -0 或非常小误差
    if (fabsf(delta) < 1e-6f) delta = 0.0f;
    if (fabsf(delta - 1.0f) < 1e-6f) delta = 1.0f;

    complexf H = 0.0f + 0.0f * I;

    for (int k = 0; k <= nsam-2; k++) {
        float Fk = compute_Fk(tau, dt, k);  // 或从临时 F[k] 数组读取

        float p1 = (float)(n0 + k) * dt;
        float p2 = (float)(n0 + k - 1) * dt;

        complexf e1 = cexpf(-I * omega * p1);
        complexf e2 = cexpf(-I * omega * p2);

        H += Fk * ((1.0f - delta) * e1 + delta * e2);
    }

    H /= dt;
    return H;
}
```

注意：

- 这里用的是 `-I * omega * p`，因为 FFTW 正变换是 \(e^{-i\omega t}\)。
- 如果后续代码使用自定义相反号的傅里叶变换，则需要整体共轭。
- 不能把 `delta` 改成 `tau/dt - floor(tau/dt)`。
- 不能遗漏 `1/dt`。

---

## 9. 频率域直接求和版本的伪代码

### 9.1 检波点端外推

```c
void apply_receiver_direct(const ExtrapCtx *ctx,
                           complexf ***U_in,
                           complexf ***U_tmp)
{
    zero_complex_array(U_tmp, ctx->ns, ctx->nh, ctx->nw);

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int is = 0; is < ctx->ns; is++) {
        for (int ih = 0; ih < ctx->nh; ih++) {

            int c = (int)(((ctx->cmp ? (ctx->s0 + is*ctx->ds) : 0.0f)
                         + ctx->h0 + ih*ctx->dh - ctx->rg0) / ctx->drg + 0.5f);
            check_receiver_table(c);

            int left  = max(0, ih - ctx->aper);
            int right = min(ctx->nh - 1, ih + ctx->aper);

            for (int ic = left; ic <= right; ic++) {

                int cc = (int)(((ctx->cmp ? (ctx->s0 + is*ctx->ds) : 0.0f)
                              + ctx->h0 + ic*ctx->dh - ctx->rg0) / ctx->drg + 0.5f);
                check_receiver_table(cc);

                float coef = 1.0f;
                coef *= (ic-left  >= ctx->tap) ? 1.0f : (float)(ic-left) / ctx->tap;
                coef *= (right-ic >= ctx->tap) ? 1.0f : (float)(right-ic) / ctx->tap;

                float tau = ctx->rtable[cc][c];
                float dist = ctx->rdatum*ctx->rdatum
                           + (ic-ih)*ctx->dh*(ic-ih)*ctx->dh;

                float geom = coef / M_PI * ctx->dh * ctx->rdatum * tau / dist;

                for (int iw = 0; iw < ctx->nw; iw++) {
                    float omega = 2.0f * M_PI * iw / (ctx->nfft * ctx->dt);
                    complexf H = compute_H_tau(omega, tau, ctx->dt, ctx->nsam);
                    complexf K = geom * H;
                    U_tmp[is][ih][iw] += K * U_in[is][ic][iw];
                }
            }
        }
    }
}
```

### 9.2 震源端外推，`cmp=0`

```c
void apply_source_direct_cmp0(const ExtrapCtx *ctx,
                              complexf ***U_tmp,
                              complexf ***U_out)
{
    zero_complex_array(U_out, ctx->ns, ctx->nh, ctx->nw);

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int ih = 0; ih < ctx->nh; ih++) {
        for (int is = 0; is < ctx->ns; is++) {

            int c = (int)((ctx->s0 + is*ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
            check_source_table(c);

            int left  = max(0, is - ctx->aper);
            int right = min(ctx->ns - 1, is + ctx->aper);

            for (int ic = left; ic <= right; ic++) {

                int cc = (int)((ctx->s0 + ic*ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
                check_source_table(cc);

                float coef = 1.0f;
                coef *= (ic-left  >= ctx->tap) ? 1.0f : (float)(ic-left) / ctx->tap;
                coef *= (right-ic >= ctx->tap) ? 1.0f : (float)(right-ic) / ctx->tap;

                float tau = ctx->stable[cc][c];
                float dist = ctx->sdatum*ctx->sdatum
                           + (ic-is)*ctx->ds*(ic-is)*ctx->ds;

                float geom = coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;

                for (int iw = 0; iw < ctx->nw; iw++) {
                    float omega = 2.0f * M_PI * iw / (ctx->nfft * ctx->dt);
                    complexf H = compute_H_tau(omega, tau, ctx->dt, ctx->nsam);
                    complexf K = geom * H;
                    U_out[is][ih][iw] += K * U_tmp[ic][ih][iw];
                }
            }
        }
    }
}
```

### 9.3 震源端外推，`cmp=1`

这部分必须保留原代码的 `ir` 共检波点循环和索引映射。

```c
void apply_source_direct_cmp1(const ExtrapCtx *ctx,
                              complexf ***U_tmp,
                              complexf ***U_out)
{
    zero_complex_array(U_out, ctx->ns, ctx->nh, ctx->nw);

    float s = fabsf((ctx->ns - 1) * ctx->ds);
    float h = fabsf((ctx->nh - 1) * ctx->dh);

    float dr;
    int jump;
    if (fabsf(ctx->ds) >= fabsf(ctx->dh)) {
        dr = fabsf(ctx->dh);
        jump = 1;
    } else {
        dr = fabsf(ctx->ds);
        jump = (int)(ctx->dh / ctx->ds + 0.5f);
    }

    int nr = (int)((s + h) / dr + 1.5f);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ir = 0; ir < nr; ir++) {

        float r = ir*dr
                + ((ctx->ds <= 0.0f) ? -1.0f : 0.0f) * s
                + ((ctx->dh <= 0.0f) ? -1.0f : 0.0f) * h;

        int sleft  = (int)((ir*dr
                    + ((ctx->ds <= 0.0f) ? -1.0f : 0.0f) * s
                    + ((ctx->ds <= 0.0f) ?  0.0f : -1.0f) * h)
                    / ctx->ds + 0.5f);

        int sright = (int)((ir*dr
                    + ((ctx->ds <= 0.0f) ? -1.0f : 0.0f) * s
                    + ((ctx->ds <= 0.0f) ? -1.0f : 0.0f) * h)
                    / ctx->ds + 0.5f);

        if (sleft < 0) sleft = 0;
        if (sright > ctx->ns - 1) sright = ctx->ns - 1;

        int left_check = (int)((r - sleft*ctx->ds) / ctx->dh + 0.5f);
        if (left_check < 0 || left_check > ctx->nh - 1) sleft++;

        int right_check = (int)((r - sright*ctx->ds) / ctx->dh + 0.5f);
        if (right_check < 0 || right_check > ctx->nh - 1) sright--;

        for (int is = sleft; is <= sright; is += jump) {

            int c = (int)((ctx->s0 + is*ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
            check_source_table(c);

            int ih = (int)((r - is*ctx->ds) / ctx->dh + 0.5f);

            int left  = max(sleft,  is - jump*ctx->aper);
            int right = min(sright, is + jump*ctx->aper);

            for (int ic = left; ic <= right; ic += jump) {

                int cc = (int)((ctx->s0 + ic*ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
                check_source_table(cc);

                int hh = (int)((r - ic*ctx->ds) / ctx->dh + 0.5f);

                float coef = 1.0f;
                coef *= (ic-left  >= ctx->tap) ? 1.0f : (float)(ic-left) / jump / ctx->tap;
                coef *= (right-ic >= ctx->tap) ? 1.0f : (float)(right-ic) / jump / ctx->tap;

                float tau = ctx->stable[cc][c];
                float dist = ctx->sdatum*ctx->sdatum
                           + (ic-is)*ctx->ds*(ic-is)*ctx->ds;

                float geom = coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;

                for (int iw = 0; iw < ctx->nw; iw++) {
                    float omega = 2.0f * M_PI * iw / (ctx->nfft * ctx->dt);
                    complexf H = compute_H_tau(omega, tau, ctx->dt, ctx->nsam);
                    complexf K = geom * H;
                    U_out[is][ih][iw] += K * U_tmp[ic][hh][iw];
                }
            }
        }
    }
}
```

注意：`cmp=1` 分支中多个 `ir` 循环理论上可能写同一个 `U_out[is][ih][iw]`。原时间域代码也这样写，但因为几何映射通常是一对一或近似一对一，未必出现明显冲突。并行化时应仔细确认是否存在写冲突；如有冲突，应加 atomic、局部缓存或调整并行层级。

---

## 10. 为蝶形算法预留的接口

完成频率域直接求和后，代码应该把空间求和抽象成接口。后续蝶形算法只替换这个接口。

### 10.1 建议的统一算子接口

```c
typedef complexf (*KernelEval)(
    const ExtrapCtx *ctx,
    int target,
    int source,
    int fixed_index,
    int iw
);
```

对于检波点端：

- `target = ih`
- `source = ic`
- `fixed_index = is`

对于 `cmp=0` 震源端：

- `target = is`
- `source = ic`
- `fixed_index = ih`

然后写一个直接求和版本：

```c
void apply_spatial_operator_direct(
    const ExtrapCtx *ctx,
    complexf **U_in_2d,
    complexf **U_out_2d,
    int n_target,
    int n_source,
    int iw,
    KernelEval eval_kernel
);
```

后续替换为蝶形算法时，只需要新增：

```c
void apply_spatial_operator_butterfly(
    const ExtrapCtx *ctx,
    complexf **U_in_2d,
    complexf **U_out_2d,
    int n_target,
    int n_source,
    int iw,
    KernelEval eval_kernel
);
```

### 10.2 蝶形算法需要看到的核形式

蝶形算法通常处理如下形式：

\[
U(r,\omega)=\sum_s A(r,s,\omega)e^{\pm i\omega t(s,r)}U(s,\omega).
\]

本问题中：

\[
A(r,s,\omega)=
w(r,s)
\frac{\Delta s}{\pi}
\frac{z_d t(s,r)}{z_d^2+|s-r|^2}
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(t)
\left[
(1-\delta)e^{\mp i\omega(n_0+k)\Delta t \pm i\omega t}
+
\delta e^{\mp i\omega(n_0+k-1)\Delta t \pm i\omega t}
\right].
\]

实现上不必强行拆相位和振幅。为了保证正确性，`eval_kernel()` 可以先直接返回完整复核：

\[
K(r,s,\omega)=
w(r,s)
\frac{\Delta s}{\pi}
\frac{z_d t(s,r)}{z_d^2+|s-r|^2}
H(\omega,t(s,r)).
\]

等直接求和版本验证通过后，再由蝶形算法内部决定如何拆分主相位和缓慢振幅。

---

## 11. 与原时间域结果对齐的验证流程

必须分阶段验证，不能直接上蝶形算法。

### 11.1 单对道验证

构造一个只含单个输入道和单个输出道的小测试：

1. 使用同一个 \(\tau\)、`dt`、`length`；
2. 时间域用原来的 `filt_set()` + `kirdat_pick()` 得到输出；
3. 频率域用 `compute_H_tau()` 乘输入 FFT 后 IFFT；
4. 对比两者。

误差应主要来自浮点精度和 FFT roundoff。

相对误差建议计算：

\[
\epsilon=\frac{\|y_{\mathrm{freq}}-y_{\mathrm{time}}\|_2}{\|y_{\mathrm{time}}\|_2+10^{-20}}.
\]

### 11.2 检波点端外推验证

只替换第一步检波点端外推，输出 `interm`，与原时间域代码写出的 `interm` 比较。

这一步验证：

- `rtable` 索引；
- `cmp` 下的接收点坐标；
- receiver aperture；
- receiver taper；
- 几何权重；
- 频率域滤波核。

### 11.3 震源端外推验证，`cmp=0`

在 `cmp=0` 情况下比较完整输出。

这一步较简单，因为共检波点道集是规则的：

\[
U_{\mathrm{out}}(i_s,i_h,\omega)=\sum_{i_c}K_s(i_s,i_c,\omega)U_r(i_c,i_h,\omega).
\]

### 11.4 震源端外推验证，`cmp=1`

最后再验证 `cmp=1`。

重点检查：

- `ir` 循环范围；
- `dr` 和 `jump`；
- `sleft`、`sright`；
- `ih=(r-is*ds)/dh+0.5`；
- `hh=(r-ic*ds)/dh+0.5`；
- 可能的并行写冲突。

### 11.5 最终可接受误差

如果频率域直接求和严格使用完整的 \(H(\omega,\tau)\)，并且 FFT 长度足够避免循环卷积，结果应与时间域版本高度一致。单精度情况下，完整数据的相对误差通常应接近：

\[
10^{-5}\sim 10^{-4}
\]

如果误差明显更大，优先检查：

1. FFT 正负号；
2. IFFT 是否除以 `nfft`；
3. `delta` 是否使用了 `ceil` 定义；
4. 是否遗漏了 `1/dt`；
5. `nfft` 是否太短导致循环卷积；
6. taper 是否完全复制原代码；
7. `cmp=1` 的 `ih/hh` 映射是否一致；
8. 是否保留了两次 `reverse_trace`。

---

## 12. 常见错误清单

### 错误 1：把核写成纯相位项

错误写法：

\[
K(r,s,\omega)=
\frac{\Delta s}{\pi}
\frac{z_d t(s,r)}{z_d^2+|s-r|^2}
e^{-i\omega t(s,r)}.
\]

这会遗漏原代码中的二阶差分滤波核和线性插值，因此不能与时间域代码对上。

正确写法：

\[
K(r,s,\omega)=
\frac{\Delta s}{\pi}
\frac{z_d t(s,r)}{z_d^2+|s-r|^2}
H(\omega,t(s,r)).
\]

### 错误 2：`delta` 定义反了

错误写法：

\[
\delta=\frac{\tau}{\Delta t}-\left\lfloor\frac{\tau}{\Delta t}\right\rfloor.
\]

正确写法：

\[
\delta=\left\lceil\frac{\tau}{\Delta t}\right\rceil-\frac{\tau}{\Delta t}.
\]

### 错误 3：FFT 长度用 `nt`

如果使用 `nfft=nt`，频域乘法会变成循环卷积，延迟尾波可能绕回到前部，导致结果不一致。

### 错误 4：忘记 IFFT 归一化

FFTW 的反变换不会自动除以 `nfft`。必须手动：

```c
tr[it] = real(ifft_result[it]) / nfft;
```

### 错误 5：直接在 `accumulate_pair_fft()` 里加蝶形算法

蝶形算法加速的是固定频率下的空间求和：

\[
U(r,\omega)=\sum_s K(r,s,\omega)U(s,\omega).
\]

不能在每一对输入道和输出道内部做蝶形算法。必须先统一 FFT 所有道，然后按频率切片做空间算子。

---

## 13. 最终目标代码结构

最终代码应呈现如下结构：

```text
main()
 ├── read_parameters_same_as_original()
 ├── open_and_read_files_same_as_original()
 ├── reverse_trace(tr_in)
 ├── setup_context()
 ├── choose_nfft()
 ├── forward_fft_all_traces(tr_in, U_in)
 │
 ├── apply_receiver_operator(U_in, U_tmp)
 │     ├── direct version now
 │     └── butterfly version later
 │
 ├── if (interm)
 │     ├── inverse_fft_all_traces(U_tmp, tr_out)
 │     └── write interm exactly as original
 │
 ├── if (cmp == 0)
 │     └── apply_source_operator_cmp0(U_tmp, U_out)
 │            ├── direct version now
 │            └── butterfly version later
 │
 ├── if (cmp == 1)
 │     └── apply_source_operator_cmp1(U_tmp, U_out)
 │            ├── direct version now
 │            └── butterfly/nonuniform version later
 │
 ├── inverse_fft_all_traces(U_out, tr_in)
 ├── reverse_trace(tr_in)
 ├── write output_file with same headers
 └── free memory
```

---

## 14. 建议的实际修改顺序

### 第一步：保留原代码，新增频率域单对道验证函数

不要马上改主流程。先写 `compute_H_tau()`，并用单对道测试确认与 `kirdat_pick()` 一致。

### 第二步：实现全局 FFT 和 IFFT

将所有输入道一次性 FFT 到频率域，再 IFFT 回来，确认误差只来自 FFT roundoff。

### 第三步：实现频率域检波点端直接求和

只替换 common-shot gather 部分，写出 `interm`，与原代码 `interm` 对比。

### 第四步：实现 `cmp=0` 的震源端直接求和

在绝对接收坐标情况下完成全流程验证。

### 第五步：实现 `cmp=1` 的震源端直接求和

保留原 `ir` gather 映射，验证完整输出。

### 第六步：抽象空间算子接口

将 `for target` / `for source` 的空间求和封装成可替换接口。

### 第七步：替换为蝶形算法

只替换空间求和函数，不改：

- 参数读取；
- 文件读写；
- FFT/IFFT；
- Green 函数索引；
- 几何权重；
- \(H(\omega,\tau)\)；
- `cmp` 坐标逻辑。

---

## 15. 一句话总结

当前代码不应直接在时间循环或单对道卷积内部加入蝶形算法。正确路线是：

\[
\boxed{
\text{时间域逐点外推}
\rightarrow
\text{全局 FFT}
\rightarrow
\text{频率域空间直接求和}
\rightarrow
\text{空间求和替换为蝶形算法}
}
\]

在这个过程中，必须使用完整的离散频率响应：

\[
\boxed{
H(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}
F_k(\tau)
\left[
(1-\delta)e^{-i\omega(n_0+k)\Delta t}
+
\delta e^{-i\omega(n_0+k-1)\Delta t}
\right]
}
\]

否则频率域结果无法与原时间域代码严格对齐。
