# 将当前频率域外推代码修改为“蝶形算法友好结构”的指导说明

本文档用于指导另一个 AI 或程序员继续修改当前的频率域波场外推代码，使其：

1. 保持当前参数、输入文件、输出文件和数据格式不变；
2. 保持与原时间域 Kirchhoff 外推代码的数值结果可对齐；
3. 将代码整理成适合后续直接替换为蝶形算法的频率域空间积分结构；
4. 先实现“频率域直接求和”等价版本，再替换空间求和为蝶形算法。

---

## 1. 当前代码的总体判断

当前代码的大方向是正确的。它已经从原来的“每一对输入道和输出道单独做 FFT 卷积”的低效形式，重构为：

```text
读取时间域输入记录
        ↓
所有道统一 FFT 到频率域
        ↓
对每个频率切片执行检波点端空间积分
        ↓
对每个频率切片执行震源端空间积分
        ↓
所有道统一 IFFT 回时间域
        ↓
写出结果
```

这正是后续使用蝶形算法所需要的形式。

当前代码中最关键的结构是：

```c
fft_traces_to_freq(tr_in, &ctx, p_f, rpad, spec, Uin);

for (int iw = 0; iw < ctx.nw; iw++) {
    apply_receiver_operator_freq_iw(&ctx, Uin, Utmp, iw);
}

for (int iw = 0; iw < ctx.nw; iw++) {
    apply_source_operator_freq_iw(&ctx, Utmp, Uout, iw);
}

ifft_traces_from_freq(Uout, &ctx, p_b, spec, rpad, tr_in);
```

这已经符合“固定频率下的空间积分矩阵作用”：

\[
U_{\rm out}(r,\omega)
=
\sum_s K(r,s,\omega)U_{\rm in}(s,\omega).
\]

因此，这套代码可以作为后续改成蝶形算法的基础。

但是，它目前还不能认为已经与原时间域代码完全等价。主要问题有两个：

1. 缺少原时间域 `kirdat_pick()` 中由早期样点截断导致的边界修正项；
2. 当前 `nfft` 对连续两次频率域外推不够稳妥，可能存在循环卷积污染风险。

因此，需要先完成本文档中的修正，再进行蝶形算法替换。

---

## 2. 必须保持不变的内容

修改代码时，以下内容必须保持不变。

### 2.1 命令行参数名称不能改变

必须继续支持：

```text
input_file=
output_file=
sgreen_file=
rgreen_file=
model_file=
cmp=
aperture=
taper=
length=
verb=
interm=
```

不能要求用户额外修改现有运行脚本。

可以新增调试参数，但新增参数必须有默认值，并且不影响原来的运行方式。例如：

```text
use_direct_freq=1
use_butterfly=0
debug_compare=0
```

### 2.2 输入数据格式不能改变

输入仍然是三维地震记录：

\[
d(t,h,s)
\]

对应头信息：

```text
n1 = nt
n2 = nh
n3 = ns
d1 = dt
d2 = dh
d3 = ds
o1 = time origin
o2 = h0
o3 = s0
```

代码中仍然按照：

```c
ctx.nt = in->headers->n[0];
ctx.nh = in->headers->n[1];
ctx.ns = in->headers->n[2];
ctx.dt = in->headers->d[0];
ctx.h0 = in->headers->o[1];
ctx.dh = in->headers->d[1];
ctx.s0 = in->headers->o[2];
ctx.ds = in->headers->d[2];
```

读取。

### 2.3 输出数据格式不能改变

输出仍然是三维数据：

\[
d_{\rm out}(t,h,s)
\]

输出头必须保持：

```c
out->headers->ndim = 3;
out->headers->n[0] = ctx.nt;
out->headers->n[1] = ctx.nh;
out->headers->n[2] = ctx.ns;
out->headers->d[0] = ctx.dt;
out->headers->d[1] = ctx.dh;
out->headers->d[2] = ctx.ds;
out->headers->o[0] = 0.0f;
out->headers->o[1] = ctx.h0;
out->headers->o[2] = ctx.s0;
```

### 2.4 `reverse_trace()` 的行为不能随意改变

当前代码读取输入后执行：

```c
reverse_trace(ctx.nt, ctx.nh, ctx.ns, tr_in);
```

写出最终结果前再次执行：

```c
reverse_trace(ctx.nt, ctx.nh, ctx.ns, tr_in);
```

如果原时间域代码也是这样处理的，则必须保留。不要在频率域修改过程中擅自删除这两次反转，否则结果会与原时间域代码不一致。

---

## 3. 当前频率域公式

### 3.1 单端平面外推公式

设：

- \(s\)：原始平面上的点；
- \(r\)：外推后平面上的点；
- \(\tau=t(s,r)\)：从 \(s\) 到 \(r\) 的走时；
- \(z_d\)：两个平面之间的外推距离；
- \(\Delta s\)：空间采样间隔；
- \(w(r,s)\)：孔径 taper 权重；
- \(R^2(s,r)=z_d^2+|s-r|^2\)。

单端外推的频率域形式为：

\[
U(r,\omega)
=
\sum_{s\in\mathcal A(r)}
K(r,s,\omega)U(s,\omega).
\]

其中：

\[
K(r,s,\omega)
=
W(r,s)H(\omega,\tau),
\]

\[
W(r,s)
=
w(r,s)
\frac{\Delta s}{\pi}
\frac{z_d \tau}{z_d^2+|s-r|^2}.
\]

这里的 \(H(\omega,\tau)\) 是与原时间域滤波器对应的频率响应。

---

## 4. 与原时间域 `kirdat_pick()` 对应的滤波核

原时间域代码中，滤波器首先构造：

\[
G_k(\tau)
=
\sqrt{
\left(
\frac{\tau+k\Delta t}{\tau}
\right)^2
-1
}.
\]

然后做一次差分：

\[
D_k=G_{k+1}-G_k.
\]

再做二次差分：

\[
F_k =
\begin{cases}
G_1-G_0, & k=0, \\
G_{k+1}-2G_k+G_{k-1}, & k=1,2,\ldots,N-2.
\end{cases}
\]

当前代码中的 `build_filter_F()` 基本正确：

```c
static void build_filter_F(float tau, float dt, int nsam, float *F)
{
    float *G = (float *)malloc((size_t)nsam * sizeof(float));
    G[0] = 0.0f;
    for (int k = 1; k < nsam; k++) {
        float t = (tau + k * dt) / tau;
        G[k] = sqrtf(t * t - 1.0f);
    }
    for (int k = 0; k < nsam - 1; k++) G[k] = G[k + 1] - G[k];
    for (int k = nsam - 2; k > 0; k--) G[k] = G[k] - G[k - 1];
    for (int k = 0; k < nsam - 1; k++) F[k] = G[k];
    free(G);
}
```

不要改变这个函数的基本逻辑，除非要做性能优化或缓存优化。

---

## 5. FFTW 符号约定

FFTW 的 `r2c` 正变换约定为：

\[
X(\omega)
=
\sum_{n=0}^{N-1}
x[n]e^{-i2\pi kn/N}.
\]

因此，如果时间域中存在延迟：

\[
x[n-p],
\]

则频率域中对应：

\[
e^{-i2\pi kp/N}X[k].
\]

所以代码中使用：

```c
cosf(-w * p)
sinf(-w * p)
```

是合理的。

注意：论文或理论说明中经常写成：

\[
\exp[i\omega t(s,r)].
\]

但在 FFTW 原生实现中，对应的代码形式通常是：

\[
\exp[-i\omega t(s,r)].
\]

不要因为公式中的正号而把代码中的 `-w*p` 改成 `+w*p`。除非同时改变整个 FFT/IFFT 的符号约定，否则会导致相位方向错误。

---

## 6. 当前代码的主要问题 1：缺少早期边界修正

### 6.1 原时间域 `kirdat_pick()` 的真实行为

原时间域代码中，`kirdat_pick()` 的核心循环是：

```c
for (isam = 0; (isam < nsam-1) && (shift-isam >= 0); isam++) {
    value += ((1.-delta)*trace[shift-isam] + delta*trace[shift-isam+1])
        *filt[its][isam]/dt;
}
```

其中：

\[
q=\frac{\tau}{\Delta t},
\]

\[
m=\lceil q\rceil,
\]

\[
\delta=m-q.
\]

对于某一个输出时间样点 \(n\)，可以理解为：

\[
\text{shift}=n-m.
\]

因此时间域结果为：

\[
y[n]
=
\sum_{\substack{k=0\\ n-m-k\ge 0}}^{N-2}
\frac{F_k}{\Delta t}
\left[
(1-\delta)x[n-m-k]
+
\delta x[n-m-k+1]
\right].
\]

由于条件 `shift-isam >= 0` 的存在，求和只允许：

\[
n-m-k \ge 0.
\]

这意味着早期样点处存在截断。

### 6.2 当前频率域代码实现的行为

当前 `build_H_iw()` 实际构造了：

\[
H(\omega,\tau)
=
\sum_{k=0}^{N-2}
\frac{F_k}{\Delta t}
\left[
(1-\delta)e^{-i\omega(m+k)\Delta t}
+
\delta e^{-i\omega(m+k-1)\Delta t}
\right].
\]

这对应完整卷积：

\[
y_{\rm freq}[n]
=
\sum_{k=0}^{N-2}
\frac{F_k}{\Delta t}
\left[
(1-\delta)x[n-m-k]
+
\delta x[n-m-k+1]
\right].
\]

如果输入 \(x[n]\) 在负时间处用零延拓，则大部分项是合理的；但是当 \(n=m+k-1\) 时，第二项会取到：

\[
x[0].
\]

而原时间域 `kirdat_pick()` 在 `shift-isam < 0` 时会完全跳过该项，因此不会保留这个 \(\delta x[0]\) 贡献。

因此，当前频率域代码比原时间域代码多了一个早期边界项。

### 6.3 需要补充的边界修正项

为了与原时间域代码严格一致，需要在频率域输出中减去这个多出来的边界贡献。

对每一对 \((r,s)\)，需要加入：

\[
E(\omega,\tau)x[0],
\]

其中：

\[
E(\omega,\tau)
=
-
\frac{\delta}{\Delta t}
\sum_{k=0}^{N-2}
F_k(\tau)
e^{-i\omega(m+k-1)\Delta t}.
\]

因此严格等价的单端频率域公式应为：

\[
U_{\rm out}(r,\omega)
=
\sum_s
W(r,s)
\left[
H(\omega,\tau)U_{\rm in}(s,\omega)
+
E(\omega,\tau)x_{\rm in}(s,0)
\right].
\]

其中：

- \(U_{\rm in}(s,\omega)\) 是输入道的频率谱；
- \(x_{\rm in}(s,0)\) 是输入道第一个时间样点；
- \(H(\omega,\tau)\) 是当前 `build_H_iw()` 已经构造的频率响应；
- \(E(\omega,\tau)\) 是新增的边界修正频率响应。

### 6.4 建议新增函数

建议新增：

```c
static void build_HE_iw(float tau,
                        const FreqContext *ctx,
                        int iw,
                        float *F,
                        float *hr, float *hi,
                        float *er, float *ei)
{
    if (tau <= 0.0f || !isfinite(tau)) {
        *hr = *hi = 0.0f;
        *er = *ei = 0.0f;
        return;
    }

    build_filter_F(tau, ctx->dt, ctx->nsam, F);

    float q = tau / ctx->dt;
    int m = (int)ceilf(q);
    float delta = (float)m - q;

    float omega = 2.0f * (float)M_PI * iw / ctx->nfft;

    float Hr = 0.0f, Hi = 0.0f;
    float Er = 0.0f, Ei = 0.0f;

    for (int k = 0; k <= ctx->nsam - 2; k++) {
        float fk = F[k] / ctx->dt;

        float p1 = (float)(m + k);
        float p2 = (float)(m + k - 1);

        /* H term */
        Hr += fk * ((1.0f - delta) * cosf(-omega * p1)
                  + delta        * cosf(-omega * p2));
        Hi += fk * ((1.0f - delta) * sinf(-omega * p1)
                  + delta        * sinf(-omega * p2));

        /* Edge correction term:
         * E = - delta * fk * exp(-i omega p2)
         */
        Er += -delta * fk * cosf(-omega * p2);
        Ei += -delta * fk * sinf(-omega * p2);
    }

    *hr = Hr;
    *hi = Hi;
    *er = Er;
    *ei = Ei;
}
```

然后原来的 `build_H_iw()` 可以替换为 `build_HE_iw()`。

### 6.5 在空间累加中使用边界修正

原来的累加是：

```c
float ar = Uin[ia + iw][0], ai = Uin[ia + iw][1];
float br = hr, bi = hi;
accumulate_pair_at_freq(ar, ai, br, bi, w,
                        &Utmp[ib + iw][0],
                        &Utmp[ib + iw][1]);
```

应改为：

```c
float xr0 = input_time_first_sample_for_this_trace;

/* H * U */
accumulate_pair_at_freq(ar, ai, hr, hi, w,
                        &Utmp[ib + iw][0],
                        &Utmp[ib + iw][1]);

/* E * x[0] */
Utmp[ib + iw][0] += w * er * xr0;
Utmp[ib + iw][1] += w * ei * xr0;
```

注意：对检波点端外推，`xr0` 应取原始输入道：

```c
xr0 = tr_in[is][ic][0];
```

但是当前代码进入频率域后，`tr_in` 后面可能被复用。为了安全，建议在 FFT 之前保存所有输入道首样点：

```c
float *Xin0 = malloc(ns * nh * sizeof(float));
Xin0[is * nh + ih] = tr_in[is][ih][0];
```

对 source pass，边界修正中的 \(x[0]\) 应该是 receiver pass 后中间道的第一个时间样点。由于 source pass 在频率域中直接使用 `Utmp`，如果要严格加入边界修正，有两种做法。

#### 做法 A：中间 IFFT 取首样点，再重新 FFT

这是最严格、最容易对齐时间域的做法：

```text
receiver pass in frequency domain
        ↓
IFFT Utmp to time-domain tr_mid
        ↓
保存 Xtmp0[is,ih] = tr_mid[is][ih][0]
        ↓
FFT tr_mid back to Utmp
        ↓
source pass with edge correction
```

这种方法最接近原时间域分步外推。

#### 做法 B：从频率域计算中间道首样点

根据 IFFT：

\[
x[0]
=
\frac{1}{N_{\rm fft}}
\sum_{k=0}^{N_{\rm fft}-1}X[k].
\]

对于 `r2c` 存储，需要考虑 Hermitian 对称：

\[
x[0]
=
\frac{1}{N_{\rm fft}}
\left[
X_0
+
X_{N/2}
+
2\sum_{k=1}^{N/2-1}\operatorname{Re}(X_k)
\right].
\]

可以写函数从 `Utmp` 直接计算中间道首样点：

```c
static float first_sample_from_r2c(const fftwf_complex *U,
                                   int base,
                                   const FreqContext *ctx)
{
    float sum = U[base + 0][0];

    if (ctx->nfft % 2 == 0) {
        sum += U[base + ctx->nw - 1][0];
        for (int iw = 1; iw < ctx->nw - 1; iw++) {
            sum += 2.0f * U[base + iw][0];
        }
    } else {
        for (int iw = 1; iw < ctx->nw; iw++) {
            sum += 2.0f * U[base + iw][0];
        }
    }

    return sum / ctx->nfft;
}
```

推荐优先使用做法 A 验证数值，再考虑做法 B 优化。

---

## 7. 当前代码的主要问题 2：`nfft` 对双端外推不够稳妥

### 7.1 当前 `choose_nfft()` 的逻辑

当前代码大致使用：

```c
int pmax = (int)ceilf(max_tau / dt) + (nsam - 2);
return next_pow2(nt + pmax + 2);
```

这对一次外推基本够用。

### 7.2 双端频率域外推的风险

当前流程是：

\[
U_{\rm out}
=
K_s(\omega)K_r(\omega)U_{\rm in}.
\]

这相当于连续做两次频率域滤波。如果 `nfft` 只按照一次滤波长度设置，第二次滤波产生的长尾可能发生 circular convolution 卷绕，从而污染前 \(nt\) 个输出样点。

### 7.3 建议修改

将 `choose_nfft()` 改成：

```c
static int choose_nfft(int nt, int nsam, float dt,
                       int nrg, int nsg,
                       float **rtable, float **stable)
{
    float max_tau = 0.0f;

    for (int i = 0; i < nrg * nrg; i++) {
        if (rtable[0][i] > max_tau) max_tau = rtable[0][i];
    }

    for (int i = 0; i < nsg * nsg; i++) {
        if (stable[0][i] > max_tau) max_tau = stable[0][i];
    }

    int pmax = (int)ceilf(max_tau / dt) + (nsam - 2);

    /*
     * Use 2*pmax because receiver-side and source-side extrapolations
     * are both applied in the same frequency-domain pipeline.
     */
    return next_pow2(nt + 2 * pmax + 2);
}
```

如果后续选择在 receiver pass 后 IFFT 并截断到 `nt`，再 FFT 进入 source pass，则可以退回到一次滤波长度。但是为了安全和简单，建议先使用 `2*pmax`。

---

## 8. 当前代码距离“可直接替换为蝶形算法”的差距

当前代码已经有蝶形算法所需的数学结构，但还没有把核函数和空间应用算子彻底解耦。

现在的 `apply_receiver_operator_freq_iw()` 内部同时做了：

1. 遍历输出点；
2. 遍历输入点；
3. 查走时表；
4. 构造核函数；
5. 执行复数乘法；
6. 累加输出。

为了后续更容易替换为蝶形算法，建议拆成下面结构：

```text
kernel evaluation
        +
spatial operator application
```

也就是：

```c
receiver_kernel(...)
source_kernel(...)
apply_direct_operator(...)
apply_butterfly_operator(...)
```

---

## 9. 建议的新代码结构

### 9.1 上下文结构体

当前 `FreqContext` 可以保留，并建议增加：

```c
typedef struct {
    int nt, nh, ns;
    int nsg, nrg;
    int aper, tap, cmp, verb;
    int nfft, nw, nsam;
    float dt, h0, dh, s0, ds, sg0, dsg, rg0, drg;
    float sdatum, rdatum;
    float **stable, **rtable;

    int use_butterfly;
    int use_edge_correction;
} FreqContext;
```

默认：

```c
ctx.use_butterfly = se_have_par("use_butterfly") ? se_get_par_int("use_butterfly") : 0;
ctx.use_edge_correction = se_have_par("use_edge_correction") ? se_get_par_int("use_edge_correction") : 1;
```

这样可以先用直接求和验证结果，再打开蝶形算法。

---

### 9.2 核函数返回结构

定义：

```c
typedef struct {
    float kr;
    float ki;
    float er;
    float ei;
    float x0;
    int valid;
} KernelValue;
```

其中：

- `kr, ki` 是 \(W H\) 的实部和虚部；
- `er, ei` 是 \(W E\) 的实部和虚部；
- `x0` 是边界修正需要的输入道首样点；
- `valid=0` 表示该输入输出对不在孔径内或表索引无效。

---

### 9.3 检波点端核函数

建议写成：

```c
static KernelValue receiver_kernel(const FreqContext *ctx,
                                   int iw,
                                   int is,
                                   int ih,
                                   int ic,
                                   const float *Xin0,
                                   HCache *cache)
{
    KernelValue kv;
    kv.kr = kv.ki = kv.er = kv.ei = 0.0f;
    kv.x0 = 0.0f;
    kv.valid = 0;

    int c = (int)((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f)
                   + ctx->h0 + ih * ctx->dh - ctx->rg0)
                   / ctx->drg) + 0.5f);

    int cc = (int)((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f)
                    + ctx->h0 + ic * ctx->dh - ctx->rg0)
                    / ctx->drg) + 0.5f);

    if (c < 0 || c >= ctx->nrg || cc < 0 || cc >= ctx->nrg) {
        return kv;
    }

    int left = (ih - ctx->aper < 0) ? 0 : ih - ctx->aper;
    int right = (ih + ctx->aper > ctx->nh - 1) ? ctx->nh - 1 : ih + ctx->aper;

    float coef = taper_weight(left, ic, right, ctx->tap, 1);
    float tau = ctx->rtable[cc][c];

    float dist = ctx->rdatum * ctx->rdatum
               + (ic - ih) * (ic - ih) * ctx->dh * ctx->dh;

    float W = coef / M_PI * ctx->dh * ctx->rdatum * tau / dist;

    if (fabsf(W) < 1e-20f) {
        return kv;
    }

    float hr, hi, er, ei;
    get_HE_from_cache_or_build(cache, tau, ctx, iw, &hr, &hi, &er, &ei);

    kv.kr = W * hr;
    kv.ki = W * hi;
    kv.er = W * er;
    kv.ei = W * ei;
    kv.x0 = Xin0[is * ctx->nh + ic];
    kv.valid = 1;

    return kv;
}
```

---

### 9.4 震源端核函数

`cmp=0` 分支比较简单：

```c
static KernelValue source_kernel_cmp0(const FreqContext *ctx,
                                      int iw,
                                      int ih,
                                      int is,
                                      int ic,
                                      const float *Xtmp0,
                                      HCache *cache)
{
    KernelValue kv;
    kv.kr = kv.ki = kv.er = kv.ei = 0.0f;
    kv.x0 = 0.0f;
    kv.valid = 0;

    int c = (int)((ctx->s0 + is * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);
    int cc = (int)((ctx->s0 + ic * ctx->ds - ctx->sg0) / ctx->dsg + 0.5f);

    if (c < 0 || c >= ctx->nsg || cc < 0 || cc >= ctx->nsg) {
        return kv;
    }

    int left = (is - ctx->aper < 0) ? 0 : is - ctx->aper;
    int right = (is + ctx->aper > ctx->ns - 1) ? ctx->ns - 1 : is + ctx->aper;

    float coef = taper_weight(left, ic, right, ctx->tap, 1);
    float tau = ctx->stable[cc][c];

    float dist = ctx->sdatum * ctx->sdatum
               + (ic - is) * (ic - is) * ctx->ds * ctx->ds;

    float W = coef / M_PI * ctx->ds * ctx->sdatum * tau / dist;

    if (fabsf(W) < 1e-20f) {
        return kv;
    }

    float hr, hi, er, ei;
    get_HE_from_cache_or_build(cache, tau, ctx, iw, &hr, &hi, &er, &ei);

    kv.kr = W * hr;
    kv.ki = W * hi;
    kv.er = W * er;
    kv.ei = W * ei;
    kv.x0 = Xtmp0[ic * ctx->nh + ih];
    kv.valid = 1;

    return kv;
}
```

`cmp=1` 分支涉及绝对检波点坐标重排：

```c
hh = (int)((r - ic * ctx->ds) / ctx->dh + 0.5f);
```

这个分支应优先保持原代码逻辑。不要为了蝶形算法先改变它。建议先完成 `cmp=0` 的标准频率域直接求和和蝶形接口，再处理 `cmp=1`。

---

## 10. 直接求和算子接口

### 10.1 检波点端直接求和

建议把 `apply_receiver_operator_freq_iw()` 改成类似：

```c
static void apply_receiver_direct_iw(const FreqContext *ctx,
                                     const fftwf_complex *Uin,
                                     fftwf_complex *Utmp,
                                     int iw,
                                     const float *Xin0)
{
    HCache cache;
    hcache_init_receiver(&cache, ctx, iw);

    for (int is = 0; is < ctx->ns; is++) {
        for (int ih = 0; ih < ctx->nh; ih++) {

            int left = (ih - ctx->aper < 0) ? 0 : ih - ctx->aper;
            int right = (ih + ctx->aper > ctx->nh - 1) ? ctx->nh - 1 : ih + ctx->aper;

            size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);

            for (int ic = left; ic <= right; ic++) {

                KernelValue kv = receiver_kernel(ctx, iw, is, ih, ic, Xin0, &cache);
                if (!kv.valid) continue;

                size_t ia = IDX3(is, ic, 0, ctx->nh, ctx->nw);

                float ar = Uin[ia + iw][0];
                float ai = Uin[ia + iw][1];

                Utmp[ib + iw][0] += kv.kr * ar - kv.ki * ai;
                Utmp[ib + iw][1] += kv.kr * ai + kv.ki * ar;

                if (ctx->use_edge_correction) {
                    Utmp[ib + iw][0] += kv.er * kv.x0;
                    Utmp[ib + iw][1] += kv.ei * kv.x0;
                }
            }
        }
    }

    hcache_free(&cache);
}
```

这样以后替换为蝶形算法时，只需要替换 `for ih / for ic` 的空间积分部分。

---

### 10.2 震源端直接求和

`cmp=0` 分支应写成：

```c
static void apply_source_direct_cmp0_iw(const FreqContext *ctx,
                                        const fftwf_complex *Utmp,
                                        fftwf_complex *Uout,
                                        int iw,
                                        const float *Xtmp0)
{
    HCache cache;
    hcache_init_source(&cache, ctx, iw);

    for (int ih = 0; ih < ctx->nh; ih++) {
        for (int is = 0; is < ctx->ns; is++) {

            int left = (is - ctx->aper < 0) ? 0 : is - ctx->aper;
            int right = (is + ctx->aper > ctx->ns - 1) ? ctx->ns - 1 : is + ctx->aper;

            size_t ib = IDX3(is, ih, 0, ctx->nh, ctx->nw);

            for (int ic = left; ic <= right; ic++) {

                KernelValue kv = source_kernel_cmp0(ctx, iw, ih, is, ic, Xtmp0, &cache);
                if (!kv.valid) continue;

                size_t ia = IDX3(ic, ih, 0, ctx->nh, ctx->nw);

                float ar = Utmp[ia + iw][0];
                float ai = Utmp[ia + iw][1];

                Uout[ib + iw][0] += kv.kr * ar - kv.ki * ai;
                Uout[ib + iw][1] += kv.kr * ai + kv.ki * ar;

                if (ctx->use_edge_correction) {
                    Uout[ib + iw][0] += kv.er * kv.x0;
                    Uout[ib + iw][1] += kv.ei * kv.x0;
                }
            }
        }
    }

    hcache_free(&cache);
}
```

`cmp=1` 先保持原循环结构，只把核函数和边界修正拆出来。

---

## 11. 蝶形算法的替换位置

完成上述重构后，直接求和部分对应：

\[
U(r,\omega)
=
\sum_s K(r,s,\omega)U(s,\omega).
\]

后续蝶形算法应替换的就是这个空间求和，而不是 FFT、IFFT、走时表读取或文件输入输出。

### 11.1 检波点端

直接求和：

```c
for (int ih = 0; ih < nh; ih++) {
    for (int ic = left; ic <= right; ic++) {
        Utmp[ih] += K_r(ih,ic,iw) * Uin[ic];
    }
}
```

应替换成：

```c
butterfly_apply_receiver_iw(
    ctx,
    iw,
    Uin_for_one_shot,
    Utmp_for_one_shot,
    receiver_kernel_callback
);
```

数学形式：

\[
U_r(i_h,\omega)
=
\sum_{i_c}
K_r(i_h,i_c,\omega)U_{\rm in}(i_c,\omega).
\]

### 11.2 震源端

直接求和：

```c
for (int is = 0; is < ns; is++) {
    for (int ic = left; ic <= right; ic++) {
        Uout[is] += K_s(is,ic,iw) * Utmp[ic];
    }
}
```

应替换成：

```c
butterfly_apply_source_iw(
    ctx,
    iw,
    Utmp_for_one_receiver,
    Uout_for_one_receiver,
    source_kernel_callback
);
```

数学形式：

\[
U_{\rm out}(i_s,\omega)
=
\sum_{i_c}
K_s(i_s,i_c,\omega)U_r(i_c,\omega).
\]

---

## 12. 推荐的修改顺序

不要一步到位直接上蝶形算法。建议严格按下面顺序修改。

### 步骤 1：保持当前频率域直接求和结构

先不要引入蝶形算法。保持：

```text
FFT all traces
receiver direct summation
source direct summation
IFFT all traces
```

只修正边界项和 `nfft`。

### 步骤 2：加入边界修正

实现 `build_HE_iw()`，并在 receiver/source 两个 pass 中加入：

\[
E(\omega,\tau)x[0].
\]

重点验证：

1. 单道单 pair 的频率域结果是否与原 `kirdat_pick()` 输出一致；
2. receiver pass 的中间结果是否与时间域 receiver pass 一致；
3. 完整双端结果是否与时间域结果一致。

### 步骤 3：扩大 `nfft`

将：

```c
return next_pow2(nt + pmax + 2);
```

改为：

```c
return next_pow2(nt + 2 * pmax + 2);
```

然后重新验证结果。

### 步骤 4：拆出核函数接口

把核函数构造和空间求和拆开。完成后，结果应与步骤 3 完全一致。

### 步骤 5：只替换 `cmp=0` 的检波点端为蝶形算法

先不要动 source pass 和 `cmp=1`。只替换：

\[
U_r(i_h,\omega)
=
\sum_{i_c}
K_r(i_h,i_c,\omega)U_{\rm in}(i_c,\omega).
\]

验证 receiver pass 中间结果。

### 步骤 6：替换 `cmp=0` 的震源端为蝶形算法

验证完整结果。

### 步骤 7：处理 `cmp=1`

`cmp=1` 涉及绝对检波点坐标重排，建议最后处理。先保证直接求和版本完全正确，再设计非规则索引或分组版本的蝶形算法。

---

## 13. 数值验证标准

### 13.1 单 pair 验证

选择一个固定输入道 \(x(t)\)，固定一个 \(\tau\)，比较：

1. 原时间域 `kirdat_pick()` 对所有输出时间样点得到的 \(y_{\rm time}(t)\)；
2. `H U + E x[0]` 经 IFFT 得到的 \(y_{\rm freq}(t)\)。

要求：

\[
\frac{\|y_{\rm freq}-y_{\rm time}\|_2}{\|y_{\rm time}\|_2}
\]

应接近单精度误差水平，通常应在：

```text
1e-5 ~ 1e-4
```

范围内，具体取决于数据量、FFT 长度和浮点累加顺序。

### 13.2 receiver pass 验证

只执行检波点端外推，写出 `interm`，与原时间域 receiver pass 的中间结果比较。

建议输出：

```text
max_abs_error
relative_l2_error
sample-by-sample correlation
```

### 13.3 full pass 验证

执行完整 receiver + source 外推，与原时间域最终结果比较。

由于频率域和时间域累加顺序不同，不要求逐样点 bitwise 一致，但应满足：

```text
relative_l2_error < 1e-4 或接近该量级
correlation > 0.9999
主要反射事件位置和振幅一致
```

如果误差集中在最前几个时间样点，优先检查边界修正项。

如果误差呈周期性或在末尾卷绕到前面，优先检查 `nfft`。

如果整体相位反了，优先检查 `exp(-i omega p)` 的符号。

---

## 14. 性能注意事项

### 14.1 不要在最内层频繁 malloc/free

当前代码在每个 `apply_*_freq_iw()` 中分配缓存是可以接受的，但后续优化时应避免在最内层 pair 循环中 malloc/free。

### 14.2 \(H(\omega,\tau)\) 应缓存

对于同一个 `iw`，同一个表索引 `(cc,c)` 对应的 \(\tau\) 会被重复使用。应该缓存：

```text
hr_cache
hi_cache
er_cache
ei_cache
cached
```

不要对每个 pair 反复调用 `build_filter_F()`。

### 14.3 蝶形算法前先保证直接求和正确

蝶形算法会引入低秩近似、插值误差和树结构误差。如果直接求和版本尚未与时间域对齐，后续无法判断误差来自公式、实现还是蝶形近似。

---

## 15. 最终目标代码结构

最终建议代码组织如下：

```text
main()
 ├── read_parameters_and_headers()
 ├── read_input_traces()
 ├── reverse_trace()
 ├── read_green_tables()
 ├── initialize_context()
 ├── choose_nfft()
 ├── save_first_samples(Xin0)
 ├── fft_traces_to_freq()
 │
 ├── for iw in frequencies:
 │       if use_butterfly:
 │           apply_receiver_butterfly_iw()
 │       else:
 │           apply_receiver_direct_iw()
 │
 ├── compute_or_reconstruct_Xtmp0()
 │
 ├── for iw in frequencies:
 │       if use_butterfly:
 │           apply_source_butterfly_iw()
 │       else:
 │           apply_source_direct_iw()
 │
 ├── ifft_traces_from_freq()
 ├── reverse_trace()
 ├── write_output()
 └── free_all()
```

其中：

```text
apply_receiver_direct_iw()
apply_receiver_butterfly_iw()
```

必须共享同一个 `receiver_kernel()`。

```text
apply_source_direct_iw()
apply_source_butterfly_iw()
```

必须共享同一个 `source_kernel()`。

这样才能保证直接求和版本和蝶形算法版本只差空间求和方法，而不是核函数定义不同。

---

## 16. 总结

当前代码已经符合昨天提出的主要方向：

\[
\text{统一 FFT}
\rightarrow
\text{固定频率下的空间积分}
\rightarrow
\text{统一 IFFT}.
\]

它已经是“蝶形算法前置重构版”。

但是，为了使结果和原时间域代码真正对上，必须先完成：

1. 加入 `kirdat_pick()` 早期样点截断对应的边界修正项；
2. 将 `nfft` 改为双端外推安全长度；
3. 把核函数构造和空间求和解耦；
4. 先验证频率域直接求和版本；
5. 再把 receiver/source 的空间求和替换成蝶形算法。

不要直接把当前 `apply_receiver_operator_freq_iw()` 或 `apply_source_operator_freq_iw()` 整体替换为蝶形算法。正确做法是先抽象出：

\[
K(r,s,\omega)
\]

然后只替换：

\[
\sum_s K(r,s,\omega)U(s,\omega)
\]

这一层空间求和。
