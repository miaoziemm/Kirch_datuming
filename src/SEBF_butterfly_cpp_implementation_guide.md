# SEBF 蝶形算法 C++ 实现指南（第二版）

本文档用于指导另一个 AI 或程序员，将已有 MATLAB 蝶形算法流程改写为 C++ 版本，并集成到项目的 `SEBF` 文件夹下，用于加速当前频率域 Kirchhoff 型平面外推代码中的空间求和。

本版指南特别强调一个关键问题：

> MATLAB 版本中的蝶形算法流程可以借鉴，但 MATLAB 示例中的外推核函数、权重函数和相位符号不能直接照搬。C++/SEBF 版本必须使用当前频率域外推代码中的核函数，即 `W(r,s) * H_fast(omega,tau)`。否则 butterfly 结果无法与当前频率域直接求和结果、以及时间域外推结果对齐。

---

## 1. 总体目标

当前频率域外推代码已经被整理为固定频率下的空间积分形式：

\[
U_{\rm out}(r,\omega)
=
\sum_s K(r,s,\omega)U_{\rm in}(s,\omega).
\]

蝶形算法只负责替换这一层空间求和：

\[
\sum_s K(r,s,\omega)U(s,\omega).
\]

蝶形算法不能改变以下内容：

1. 输入/输出文件格式；
2. FFT/IFFT 流程；
3. Green 表读取方式；
4. 走时表索引方式；
5. 几何权重；
6. taper；
7. 当前频率域外推核函数；
8. `cmp=0/1` 的坐标逻辑；
9. 当前 `nfft` 策略；
10. 当前快速版中是否使用边界修正。

换句话说，direct 版本和 butterfly 版本必须只在“空间求和方法”上不同：

```text
direct:     显式双重循环求和
butterfly:  蝶形算法近似求和
```

其余内容必须一致。

---

## 2. MATLAB 代码能借鉴什么，不能照搬什么

### 2.1 可以借鉴的部分

MATLAB 版本的蝶形算法通常包含以下步骤，这些流程可以迁移到 C++：

1. 构造 source tree；
2. 构造 target tree；
3. 在叶子节点进行初始化；
4. 使用 Chebyshev 节点；
5. 使用 barycentric Lagrange 插值；
6. 第一阶段递推；
7. 中间层 switch；
8. 第二阶段递推；
9. 叶子节点重构；
10. 与 direct summation 对比误差。

这些是蝶形算法的结构，可以照着改写。

### 2.2 不能照搬的部分

MATLAB 代码中的测试核函数一般类似：

\[
K_{\rm matlab}(r,s,\omega)
=
c_\theta(r,s)\exp(i\omega t(s,r)).
\]

该核只是 MATLAB 示例中用于演示蝶形算法的振荡核，不是当前 C/C++ 频率域 Kirchhoff 外推使用的真实核。

不能直接照搬：

\[
c_\theta(r,s),
\]

不能直接照搬：

\[
\exp(i\omega t(s,r)),
\]

也不能直接照搬 MATLAB 中可能出现的：

\[
f(s) \leftarrow i\omega f(s)/v.
\]

这些操作会改变外推核函数，导致结果与当前 C/C++ 频率域外推结果不一致。

---

## 3. 当前频率域外推必须使用的核函数

当前快速频率域外推代码使用的单端核函数应写为：

\[
K(r,s,\omega)
=
W(r,s)H_{\rm fast}(\omega,t(s,r)).
\]

其中：

\[
W(r,s)
=
w_{\rm taper}(r,s)
\frac{\Delta s}{\pi}
\frac{z_d\,t(s,r)}
{z_d^2+|s-r|^2}.
\]

这里：

- \(r\)：外推后平面上的输出点；
- \(s\)：原始平面上的输入点；
- \(t(s,r)\)：从输入点 \(s\) 到输出点 \(r\) 的走时；
- \(z_d\)：平面间距，receiver pass 中为 `rdatum`，source pass 中为 `sdatum`；
- \(\Delta s\)：输入点所在平面的采样间隔；
- \(w_{\rm taper}\)：当前代码中的孔径 taper；
- \(H_{\rm fast}\)：当前快速频率域外推代码中的时间滤波频率响应。

---

## 4. 当前快速版中的 \(H_{\rm fast}\)

当前快速版对应的频率响应为：

\[
H_{\rm fast}(\omega,\tau)
=
\left[
(1-\delta)e^{-i\omega_d m}
+
\delta e^{-i\omega_d(m-1)}
\right]
\sum_{k=0}^{N-2}
\frac{F_k(\tau)}{\Delta t}
e^{-i\omega_d k}.
\]

其中：

\[
\omega_d = \frac{2\pi i_\omega}{N_{\rm fft}},
\]

\[
q=\frac{\tau}{\Delta t},
\]

\[
m=\lceil q\rceil,
\]

\[
\delta=m-q.
\]

滤波器 \(F_k(\tau)\) 来自原时间域代码中的二阶差分核：

\[
G_k(\tau)
=
\sqrt{
\left(
\frac{\tau+k\Delta t}{\tau}
\right)^2-1
},
\]

\[
F_0(\tau)=G_1(\tau)-G_0(\tau),
\]

\[
F_k(\tau)=G_{k+1}(\tau)-2G_k(\tau)+G_{k-1}(\tau),
\qquad k=1,2,\ldots,N-2.
\]

注意：这里使用的是 FFTW 的正变换约定，因此时间延迟对应：

\[
e^{-i\omega t}.
\]

不要把 MATLAB 里的：

\[
e^{+i\omega t}
\]

直接照搬到 C++ 中，否则相位会反。

---

## 5. direct summation 与 butterfly 的验证关系

正确的验证流程不是直接比较：

```text
MATLAB butterfly kernel result  vs  time-domain extrapolation result
```

而是：

```text
Step 1:
C++ direct summation using K_cpp
        vs
C++ butterfly summation using the same K_cpp

Step 2:
C++ frequency-domain direct/butterfly result after IFFT
        vs
current fast frequency-domain extrapolation output

Step 3:
current fast frequency-domain output
        vs
original time-domain extrapolation output
```

其中最重要的是 Step 1。

蝶形算法必须先证明：

\[
U_{\rm out}^{\rm butterfly}(r,\omega)
\approx
U_{\rm out}^{\rm direct}(r,\omega)
\]

对于同一个核函数：

\[
K_{\rm cpp}(r,s,\omega)
=
W(r,s)H_{\rm fast}(\omega,t(s,r)).
\]

如果这一步对不上，说明蝶形算法实现或核函数接口有问题。

---

## 6. 误差指标

建议在测试中输出：

### 6.1 复数频率域误差

\[
\epsilon_{\rm rel}
=
\frac{
\|U_{\rm butterfly}-U_{\rm direct}\|_2
}{
\|U_{\rm direct}\|_2+\epsilon
}.
\]

建议阈值：

```text
1e-3 ~ 1e-5
```

具体取决于蝶形算法的插值阶数、分块层数和 Chebyshev 点数。

### 6.2 最大绝对误差

\[
\epsilon_{\max}
=
\max_r
|U_{\rm butterfly}(r)-U_{\rm direct}(r)|.
\]

### 6.3 时间域误差

所有频率完成后，经 IFFT 得到时间域结果，计算：

\[
\epsilon_{\rm time}
=
\frac{
\|d_{\rm butterfly}(t,r)-d_{\rm direct}(t,r)\|_2
}{
\|d_{\rm direct}(t,r)\|_2+\epsilon
}.
\]

---

## 7. 推荐 SEBF 文件结构

所有文件放到项目的 `SEBF` 文件夹下。

建议结构如下：

```text
SEBF/
├── include/
│   ├── sebf_complex.h
│   ├── sebf_context.h
│   ├── sebf_kernel.h
│   ├── sebf_tree.h
│   ├── sebf_chebyshev.h
│   ├── sebf_butterfly.h
│   └── sebf_direct.h
├── src/
│   ├── sebf_kernel.cpp
│   ├── sebf_tree.cpp
│   ├── sebf_chebyshev.cpp
│   ├── sebf_butterfly.cpp
│   └── sebf_direct.cpp
└── tests/
    └── test_sebf_direct_vs_butterfly.cpp
```

下面给出每个文件的建议内容。

---

# 8. `SEBF/include/sebf_complex.h`

```cpp
#ifndef SEBF_COMPLEX_H
#define SEBF_COMPLEX_H

#include <cmath>

namespace sebf {

struct Complexf {
    float r;
    float i;

    Complexf() : r(0.0f), i(0.0f) {}
    Complexf(float real, float imag) : r(real), i(imag) {}

    Complexf& operator+=(const Complexf& b) {
        r += b.r;
        i += b.i;
        return *this;
    }

    Complexf& operator-=(const Complexf& b) {
        r -= b.r;
        i -= b.i;
        return *this;
    }
};

inline Complexf operator+(const Complexf& a, const Complexf& b) {
    return Complexf(a.r + b.r, a.i + b.i);
}

inline Complexf operator-(const Complexf& a, const Complexf& b) {
    return Complexf(a.r - b.r, a.i - b.i);
}

inline Complexf operator*(const Complexf& a, const Complexf& b) {
    return Complexf(a.r * b.r - a.i * b.i,
                    a.r * b.i + a.i * b.r);
}

inline Complexf operator*(float a, const Complexf& b) {
    return Complexf(a * b.r, a * b.i);
}

inline Complexf conj(const Complexf& a) {
    return Complexf(a.r, -a.i);
}

inline float abs2(const Complexf& a) {
    return a.r * a.r + a.i * a.i;
}

inline Complexf exp_i(float phase) {
    return Complexf(std::cos(phase), std::sin(phase));
}

} // namespace sebf

#endif
```

---

# 9. `SEBF/include/sebf_context.h`

```cpp
#ifndef SEBF_CONTEXT_H
#define SEBF_CONTEXT_H

namespace sebf {

struct FreqContextView {
    int nt;
    int nh;
    int ns;
    int nsg;
    int nrg;
    int aper;
    int tap;
    int cmp;
    int nfft;
    int nw;
    int nsam;

    float dt;
    float h0;
    float dh;
    float s0;
    float ds;
    float sg0;
    float dsg;
    float rg0;
    float drg;
    float sdatum;
    float rdatum;

    float** stable;
    float** rtable;
};

struct ButterflyOptions {
    int order;          // Chebyshev interpolation order
    int leaf_size;      // leaf box size
    int max_level;      // if <=0, determine automatically
    int use_aperture;   // 1: kernel returns zero outside aperture
    int verbose;
};

} // namespace sebf

#endif
```

---

# 10. `SEBF/include/sebf_kernel.h`

```cpp
#ifndef SEBF_KERNEL_H
#define SEBF_KERNEL_H

#include "sebf_complex.h"
#include "sebf_context.h"

namespace sebf {

enum class KernelType {
    Receiver,
    SourceCmp0,
    SourceCmp1
};

struct KernelRequest {
    KernelType type;
    int iw;

    /*
     * gather_id meaning:
     * Receiver:   gather_id = is
     * SourceCmp0: gather_id = ih
     * SourceCmp1: gather_id = ir or user-defined common receiver index
     */
    int gather_id;

    /*
     * iout and iin are local operator indices.
     * Receiver:
     *     iout = ih, iin = ic
     * SourceCmp0:
     *     iout = is, iin = ic
     */
    int iout;
    int iin;
};

class KernelEvaluator {
public:
    explicit KernelEvaluator(const FreqContextView* ctx);

    Complexf eval_receiver(int iw, int is, int ih, int ic) const;

    Complexf eval_source_cmp0(int iw, int ih, int is, int ic) const;

    /*
     * Source cmp=1 is irregular. It should be implemented later after
     * cmp=0 has been verified.
     */
    Complexf eval_source_cmp1(int iw,
                              int ir,
                              int is,
                              int ih,
                              int ic,
                              int hh,
                              int left,
                              int right,
                              int jump) const;

private:
    const FreqContextView* ctx_;

    float taper_weight(int left, int center, int right, int tap, int stride) const;

    void build_filter_F(float tau, float* F) const;

    Complexf build_H_fast(float tau, int iw) const;
};

} // namespace sebf

#endif
```

---

# 11. `SEBF/src/sebf_kernel.cpp`

```cpp
#include "sebf_kernel.h"

#include <cmath>
#include <cstdlib>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace sebf {

KernelEvaluator::KernelEvaluator(const FreqContextView* ctx)
    : ctx_(ctx)
{
}

float KernelEvaluator::taper_weight(int left,
                                    int center,
                                    int right,
                                    int tap,
                                    int stride) const
{
    if (tap <= 0) return 1.0f;

    /*
     * Keep integer-division semantics used by the current C/C++ code.
     */
    float wl = (center - left >= tap)
             ? 1.0f
             : static_cast<float>(((center - left) / stride) / tap);

    float wr = (right - center >= tap)
             ? 1.0f
             : static_cast<float>(((right - center) / stride) / tap);

    return wl * wr;
}

void KernelEvaluator::build_filter_F(float tau, float* F) const
{
    const int nsam = ctx_->nsam;
    const float dt = ctx_->dt;

    std::vector<float> G(nsam, 0.0f);

    G[0] = 0.0f;

    for (int k = 1; k < nsam; k++) {
        float x = (tau + k * dt) / tau;
        G[k] = std::sqrt(x * x - 1.0f);
    }

    for (int k = 0; k < nsam - 1; k++) {
        G[k] = G[k + 1] - G[k];
    }

    for (int k = nsam - 2; k > 0; k--) {
        G[k] = G[k] - G[k - 1];
    }

    for (int k = 0; k < nsam - 1; k++) {
        F[k] = G[k];
    }
}

Complexf KernelEvaluator::build_H_fast(float tau, int iw) const
{
    if (tau <= 0.0f || !std::isfinite(tau)) {
        return Complexf(0.0f, 0.0f);
    }

    const int nsam = ctx_->nsam;
    const float dt = ctx_->dt;
    const int nfft = ctx_->nfft;

    std::vector<float> F(nsam - 1, 0.0f);
    build_filter_F(tau, F.data());

    const float q = tau / dt;
    const int m = static_cast<int>(std::ceil(q));
    const float delta = static_cast<float>(m) - q;

    const float wd = 2.0f * static_cast<float>(M_PI) * iw / nfft;

    /*
     * H_fast =
     * [(1-delta) exp(-i wd m) + delta exp(-i wd (m-1))]
     * * sum_k F_k/dt * exp(-i wd k)
     */

    Complexf sumF(0.0f, 0.0f);

    /*
     * Use recurrence for exp(-i wd k).
     */
    Complexf step(std::cos(-wd), std::sin(-wd));
    Complexf phase_k(1.0f, 0.0f);

    for (int k = 0; k <= nsam - 2; k++) {
        float fk = F[k] / dt;
        sumF += fk * phase_k;
        phase_k = phase_k * step;
    }

    Complexf phase_m(std::cos(-wd * m), std::sin(-wd * m));
    Complexf phase_m1(std::cos(-wd * (m - 1)), std::sin(-wd * (m - 1)));

    Complexf interp = (1.0f - delta) * phase_m + delta * phase_m1;

    return interp * sumF;
}

Complexf KernelEvaluator::eval_receiver(int iw, int is, int ih, int ic) const
{
    const FreqContextView* ctx = ctx_;

    int left = (ih - ctx->aper < 0) ? 0 : ih - ctx->aper;
    int right = (ih + ctx->aper > ctx->nh - 1) ? ctx->nh - 1 : ih + ctx->aper;

    if (ic < left || ic > right) {
        return Complexf(0.0f, 0.0f);
    }

    int c = static_cast<int>((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f)
                              + ctx->h0 + ih * ctx->dh - ctx->rg0)
                              / ctx->drg) + 0.5f);

    int cc = static_cast<int>((((ctx->cmp ? (ctx->s0 + is * ctx->ds) : 0.0f)
                               + ctx->h0 + ic * ctx->dh - ctx->rg0)
                               / ctx->drg) + 0.5f);

    if (c < 0 || c >= ctx->nrg || cc < 0 || cc >= ctx->nrg) {
        return Complexf(0.0f, 0.0f);
    }

    float coef = taper_weight(left, ic, right, ctx->tap, 1);
    float tau = ctx->rtable[cc][c];

    float dist = ctx->rdatum * ctx->rdatum
               + (ic - ih) * (ic - ih) * ctx->dh * ctx->dh;

    float W = coef / static_cast<float>(M_PI)
            * ctx->dh * ctx->rdatum * tau / dist;

    Complexf H = build_H_fast(tau, iw);

    return W * H;
}

Complexf KernelEvaluator::eval_source_cmp0(int iw, int ih, int is, int ic) const
{
    const FreqContextView* ctx = ctx_;

    int left = (is - ctx->aper < 0) ? 0 : is - ctx->aper;
    int right = (is + ctx->aper > ctx->ns - 1) ? ctx->ns - 1 : is + ctx->aper;

    if (ic < left || ic > right) {
        return Complexf(0.0f, 0.0f);
    }

    int c = static_cast<int>((ctx->s0 + is * ctx->ds - ctx->sg0)
                             / ctx->dsg + 0.5f);

    int cc = static_cast<int>((ctx->s0 + ic * ctx->ds - ctx->sg0)
                              / ctx->dsg + 0.5f);

    if (c < 0 || c >= ctx->nsg || cc < 0 || cc >= ctx->nsg) {
        return Complexf(0.0f, 0.0f);
    }

    float coef = taper_weight(left, ic, right, ctx->tap, 1);
    float tau = ctx->stable[cc][c];

    float dist = ctx->sdatum * ctx->sdatum
               + (ic - is) * (ic - is) * ctx->ds * ctx->ds;

    float W = coef / static_cast<float>(M_PI)
            * ctx->ds * ctx->sdatum * tau / dist;

    Complexf H = build_H_fast(tau, iw);

    return W * H;
}

Complexf KernelEvaluator::eval_source_cmp1(int iw,
                                           int ir,
                                           int is,
                                           int ih,
                                           int ic,
                                           int hh,
                                           int left,
                                           int right,
                                           int jump) const
{
    const FreqContextView* ctx = ctx_;

    int c = static_cast<int>((ctx->s0 + is * ctx->ds - ctx->sg0)
                             / ctx->dsg + 0.5f);

    int cc = static_cast<int>((ctx->s0 + ic * ctx->ds - ctx->sg0)
                              / ctx->dsg + 0.5f);

    if (c < 0 || c >= ctx->nsg || cc < 0 || cc >= ctx->nsg) {
        return Complexf(0.0f, 0.0f);
    }

    float coef = taper_weight(left, ic, right, ctx->tap, jump);
    float tau = ctx->stable[cc][c];

    float dist = ctx->sdatum * ctx->sdatum
               + (ic - is) * (ic - is) * ctx->ds * ctx->ds;

    float W = coef / static_cast<float>(M_PI)
            * ctx->ds * ctx->sdatum * tau / dist;

    Complexf H = build_H_fast(tau, iw);

    return W * H;
}

} // namespace sebf
```

---

# 12. `SEBF/include/sebf_direct.h`

```cpp
#ifndef SEBF_DIRECT_H
#define SEBF_DIRECT_H

#include "sebf_complex.h"
#include "sebf_context.h"
#include "sebf_kernel.h"

namespace sebf {

void direct_apply_receiver_iw(const FreqContextView* ctx,
                              int iw,
                              const Complexf* Uin,
                              Complexf* Utmp);

void direct_apply_source_cmp0_iw(const FreqContextView* ctx,
                                 int iw,
                                 const Complexf* Utmp,
                                 Complexf* Uout);

} // namespace sebf

#endif
```

---

# 13. `SEBF/src/sebf_direct.cpp`

```cpp
#include "sebf_direct.h"

namespace sebf {

static inline size_t idx3(int is, int ih, int iw, int nh, int nw)
{
    return (static_cast<size_t>(is) * nh + ih) * nw + iw;
}

void direct_apply_receiver_iw(const FreqContextView* ctx,
                              int iw,
                              const Complexf* Uin,
                              Complexf* Utmp)
{
    KernelEvaluator kernel(ctx);

    for (int is = 0; is < ctx->ns; is++) {
        for (int ih = 0; ih < ctx->nh; ih++) {
            Complexf sum(0.0f, 0.0f);

            int left = (ih - ctx->aper < 0) ? 0 : ih - ctx->aper;
            int right = (ih + ctx->aper > ctx->nh - 1) ? ctx->nh - 1 : ih + ctx->aper;

            for (int ic = left; ic <= right; ic++) {
                Complexf K = kernel.eval_receiver(iw, is, ih, ic);
                Complexf x = Uin[idx3(is, ic, iw, ctx->nh, ctx->nw)];
                sum += K * x;
            }

            Utmp[idx3(is, ih, iw, ctx->nh, ctx->nw)] += sum;
        }
    }
}

void direct_apply_source_cmp0_iw(const FreqContextView* ctx,
                                 int iw,
                                 const Complexf* Utmp,
                                 Complexf* Uout)
{
    KernelEvaluator kernel(ctx);

    for (int ih = 0; ih < ctx->nh; ih++) {
        for (int is = 0; is < ctx->ns; is++) {
            Complexf sum(0.0f, 0.0f);

            int left = (is - ctx->aper < 0) ? 0 : is - ctx->aper;
            int right = (is + ctx->aper > ctx->ns - 1) ? ctx->ns - 1 : is + ctx->aper;

            for (int ic = left; ic <= right; ic++) {
                Complexf K = kernel.eval_source_cmp0(iw, ih, is, ic);
                Complexf x = Utmp[idx3(ic, ih, iw, ctx->nh, ctx->nw)];
                sum += K * x;
            }

            Uout[idx3(is, ih, iw, ctx->nh, ctx->nw)] += sum;
        }
    }
}

} // namespace sebf
```

---

# 14. `SEBF/include/sebf_chebyshev.h`

```cpp
#ifndef SEBF_CHEBYSHEV_H
#define SEBF_CHEBYSHEV_H

#include <vector>

namespace sebf {

void chebyshev_nodes(int order, float a, float b, std::vector<float>& nodes);

void barycentric_weights(int order, std::vector<float>& weights);

float lagrange_basis_value(const std::vector<float>& nodes,
                           const std::vector<float>& bary_w,
                           int j,
                           float x);

} // namespace sebf

#endif
```

---

# 15. `SEBF/src/sebf_chebyshev.cpp`

```cpp
#include "sebf_chebyshev.h"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace sebf {

void chebyshev_nodes(int order, float a, float b, std::vector<float>& nodes)
{
    nodes.resize(order);

    float c = 0.5f * (a + b);
    float h = 0.5f * (b - a);

    for (int j = 0; j < order; j++) {
        float x = std::cos(static_cast<float>(M_PI) * (2.0f * j + 1.0f) / (2.0f * order));
        nodes[j] = c + h * x;
    }
}

void barycentric_weights(int order, std::vector<float>& weights)
{
    weights.resize(order);

    for (int j = 0; j < order; j++) {
        /*
         * Chebyshev first-kind nodes have alternating barycentric weights.
         * For this implementation, a normalized alternating form is sufficient.
         */
        weights[j] = (j % 2 == 0) ? 1.0f : -1.0f;
    }
}

float lagrange_basis_value(const std::vector<float>& nodes,
                           const std::vector<float>& bary_w,
                           int j,
                           float x)
{
    const int n = static_cast<int>(nodes.size());
    const float eps = 1e-7f;

    for (int k = 0; k < n; k++) {
        if (std::fabs(x - nodes[k]) < eps) {
            return (k == j) ? 1.0f : 0.0f;
        }
    }

    float denom = 0.0f;
    for (int k = 0; k < n; k++) {
        denom += bary_w[k] / (x - nodes[k]);
    }

    float numer = bary_w[j] / (x - nodes[j]);

    return numer / denom;
}

} // namespace sebf
```

---

# 16. `SEBF/include/sebf_tree.h`

```cpp
#ifndef SEBF_TREE_H
#define SEBF_TREE_H

#include <vector>

namespace sebf {

struct Box1D {
    int begin;
    int end;        // exclusive
    float xmin;
    float xmax;
    int level;
    int parent;
    int left;
    int right;

    bool is_leaf() const {
        return left < 0 && right < 0;
    }

    int size() const {
        return end - begin;
    }
};

class Tree1D {
public:
    void build_uniform(int n, float x0, float dx, int leaf_size);

    const std::vector<Box1D>& boxes() const { return boxes_; }
    const std::vector<int>& leaves() const { return leaves_; }
    int root() const { return 0; }

private:
    int build_recursive(int begin,
                        int end,
                        float x0,
                        float dx,
                        int level,
                        int parent,
                        int leaf_size);

private:
    std::vector<Box1D> boxes_;
    std::vector<int> leaves_;
};

} // namespace sebf

#endif
```

---

# 17. `SEBF/src/sebf_tree.cpp`

```cpp
#include "sebf_tree.h"

namespace sebf {

void Tree1D::build_uniform(int n, float x0, float dx, int leaf_size)
{
    boxes_.clear();
    leaves_.clear();
    build_recursive(0, n, x0, dx, 0, -1, leaf_size);
}

int Tree1D::build_recursive(int begin,
                            int end,
                            float x0,
                            float dx,
                            int level,
                            int parent,
                            int leaf_size)
{
    int id = static_cast<int>(boxes_.size());

    Box1D box;
    box.begin = begin;
    box.end = end;
    box.xmin = x0 + begin * dx;
    box.xmax = x0 + (end - 1) * dx;
    box.level = level;
    box.parent = parent;
    box.left = -1;
    box.right = -1;

    boxes_.push_back(box);

    if (end - begin <= leaf_size) {
        leaves_.push_back(id);
        return id;
    }

    int mid = begin + (end - begin) / 2;

    int left_id = build_recursive(begin, mid, x0, dx, level + 1, id, leaf_size);
    int right_id = build_recursive(mid, end, x0, dx, level + 1, id, leaf_size);

    boxes_[id].left = left_id;
    boxes_[id].right = right_id;

    return id;
}

} // namespace sebf
```

---

# 18. `SEBF/include/sebf_butterfly.h`

```cpp
#ifndef SEBF_BUTTERFLY_H
#define SEBF_BUTTERFLY_H

#include "sebf_complex.h"
#include "sebf_context.h"
#include "sebf_kernel.h"

namespace sebf {

/*
 * This is the main interface expected by the extrapolation code.
 * It computes y = K * x for one fixed frequency and one fixed gather.
 */
void butterfly_apply_1d(
    int n_in,
    int n_out,
    float x_in0,
    float dx_in,
    float x_out0,
    float dx_out,
    int iw,
    int gather_id,
    const Complexf* x,
    Complexf* y,
    const FreqContextView* ctx,
    KernelType kernel_type,
    const ButterflyOptions& opt);

/*
 * Receiver-side wrapper:
 * For fixed iw and fixed is:
 * y[ih] = sum_ic K_receiver(is, ih, ic, iw) x[ic].
 */
void butterfly_apply_receiver_gather_iw(
    const FreqContextView* ctx,
    int iw,
    int is,
    const Complexf* x,
    Complexf* y,
    const ButterflyOptions& opt);

/*
 * Source-side cmp=0 wrapper:
 * For fixed iw and fixed ih:
 * y[is] = sum_ic K_source_cmp0(ih, is, ic, iw) x[ic].
 */
void butterfly_apply_source_cmp0_gather_iw(
    const FreqContextView* ctx,
    int iw,
    int ih,
    const Complexf* x,
    Complexf* y,
    const ButterflyOptions& opt);

} // namespace sebf

#endif
```

---

# 19. `SEBF/src/sebf_butterfly.cpp`

下面给出一个“接口完整、可直接替换 direct 的骨架版本”。

注意：为了让其他 AI 能先完成可运行版本，这里先提供一个 fallback direct implementation。真正的蝶形递推应在此基础上替换 `butterfly_apply_1d()` 内部实现，但函数接口不变。

```cpp
#include "sebf_butterfly.h"

#include <cstring>
#include <vector>

namespace sebf {

static Complexf eval_kernel_by_type(const KernelEvaluator& kernel,
                                    KernelType type,
                                    int iw,
                                    int gather_id,
                                    int iout,
                                    int iin)
{
    if (type == KernelType::Receiver) {
        int is = gather_id;
        int ih = iout;
        int ic = iin;
        return kernel.eval_receiver(iw, is, ih, ic);
    }

    if (type == KernelType::SourceCmp0) {
        int ih = gather_id;
        int is = iout;
        int ic = iin;
        return kernel.eval_source_cmp0(iw, ih, is, ic);
    }

    /*
     * SourceCmp1 is intentionally not implemented in the generic regular
     * 1D interface. It requires an irregular gather interface.
     */
    return Complexf(0.0f, 0.0f);
}

void butterfly_apply_1d(
    int n_in,
    int n_out,
    float x_in0,
    float dx_in,
    float x_out0,
    float dx_out,
    int iw,
    int gather_id,
    const Complexf* x,
    Complexf* y,
    const FreqContextView* ctx,
    KernelType kernel_type,
    const ButterflyOptions& opt)
{
    /*
     * Stage 0: fallback direct implementation.
     *
     * This must be kept as a debug mode even after the real butterfly
     * algorithm is implemented, because it provides an exact reference
     * with the same interface and same kernel evaluator.
     */
    KernelEvaluator kernel(ctx);

    for (int iout = 0; iout < n_out; iout++) {
        y[iout] = Complexf(0.0f, 0.0f);

        for (int iin = 0; iin < n_in; iin++) {
            Complexf K = eval_kernel_by_type(kernel,
                                             kernel_type,
                                             iw,
                                             gather_id,
                                             iout,
                                             iin);

            y[iout] += K * x[iin];
        }
    }

    /*
     * Real butterfly implementation should replace the double loop above.
     *
     * Required algorithmic stages:
     *
     * 1. Build source tree and target tree.
     * 2. For each source leaf and target root, initialize equivalent
     *    coefficients on Chebyshev nodes.
     * 3. Upward recursion on the source tree while descending the target tree.
     * 4. Switch at the middle level.
     * 5. Downward recursion on the target tree while source boxes coarsen.
     * 6. At target leaves, evaluate the interpolant at physical output points.
     *
     * Important:
     *     All kernel evaluations must call eval_kernel_by_type() or
     *     KernelEvaluator directly. Do not use MATLAB's ctheta*exp(iwt) kernel.
     */
}

void butterfly_apply_receiver_gather_iw(
    const FreqContextView* ctx,
    int iw,
    int is,
    const Complexf* x,
    Complexf* y,
    const ButterflyOptions& opt)
{
    butterfly_apply_1d(ctx->nh,
                       ctx->nh,
                       ctx->h0,
                       ctx->dh,
                       ctx->h0,
                       ctx->dh,
                       iw,
                       is,
                       x,
                       y,
                       ctx,
                       KernelType::Receiver,
                       opt);
}

void butterfly_apply_source_cmp0_gather_iw(
    const FreqContextView* ctx,
    int iw,
    int ih,
    const Complexf* x,
    Complexf* y,
    const ButterflyOptions& opt)
{
    butterfly_apply_1d(ctx->ns,
                       ctx->ns,
                       ctx->s0,
                       ctx->ds,
                       ctx->s0,
                       ctx->ds,
                       iw,
                       ih,
                       x,
                       y,
                       ctx,
                       KernelType::SourceCmp0,
                       opt);
}

} // namespace sebf
```

---

## 20. 如何把真实蝶形递推填入 `butterfly_apply_1d()`

这一节说明真实 butterfly 算法应该怎么写。这里使用数学流程，不强制给出唯一实现。

### 20.1 算子

固定一个频率和一个 gather 后，目标是：

\[
y(r_i)=\sum_{s_j}K(r_i,s_j,\omega)x(s_j).
\]

### 20.2 树结构

构造输入点树 \(T_s\) 和输出点树 \(T_r\)。

输入点：

\[
s_j=x_{\rm in,0}+j\Delta x_{\rm in}.
\]

输出点：

\[
r_i=x_{\rm out,0}+i\Delta x_{\rm out}.
\]

### 20.3 Chebyshev 节点

每个 box 内使用 \(p\) 个 Chebyshev 节点：

\[
x_q = \frac{a+b}{2} + \frac{b-a}{2}
\cos\left(\frac{(2q+1)\pi}{2p}\right).
\]

### 20.4 初始化

对 source leaf \(A\) 和 target root \(B\)，计算源方向等效系数。

初始化必须使用当前核函数：

\[
K_{\rm cpp}(r,s,\omega)
=
W(r,s)H_{\rm fast}(\omega,t(s,r)).
\]

不能使用：

\[
c_\theta(r,s)e^{i\omega t(s,r)}.
\]

### 20.5 递推

蝶形算法递推的基本思想是：

- source box 逐层变大；
- target box 逐层变小；
- 每一对互补尺度的 box 使用低秩插值表示；
- 中间层进行一次表示切换；
- 最终在 target leaf 上重构输出。

### 20.6 终止重构

对于每个 target leaf 中的真实输出点 \(r_i\)，使用 target box Chebyshev 节点上的系数插值得到：

\[
y(r_i).
\]

### 20.7 实现要求

真实 butterfly 实现必须满足：

1. direct 和 butterfly 使用同一个 `KernelEvaluator`；
2. 不能改变 `H_fast` 的相位符号；
3. 不能使用 MATLAB 示例权重；
4. 不允许在 butterfly 内部另写一套走时或几何权重；
5. 允许使用 Chebyshev 插值近似空间求和；
6. 支持 `order` 参数控制精度；
7. 支持 `leaf_size` 参数控制树叶大小。

---

# 21. `SEBF/tests/test_sebf_direct_vs_butterfly.cpp`

这个测试用于验证 butterfly 接口。初始阶段由于 `butterfly_apply_1d()` fallback 到 direct，误差应接近 0。等真实 butterfly 替换进去后，误差应保持较低水平。

```cpp
#include "sebf_complex.h"
#include "sebf_context.h"
#include "sebf_direct.h"
#include "sebf_butterfly.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace sebf;

static float rel_l2_error(const std::vector<Complexf>& a,
                          const std::vector<Complexf>& b)
{
    double num = 0.0;
    double den = 0.0;

    for (size_t i = 0; i < a.size(); i++) {
        double dr = static_cast<double>(a[i].r - b[i].r);
        double di = static_cast<double>(a[i].i - b[i].i);
        num += dr * dr + di * di;

        double br = static_cast<double>(b[i].r);
        double bi = static_cast<double>(b[i].i);
        den += br * br + bi * bi;
    }

    return static_cast<float>(std::sqrt(num / (den + 1e-30)));
}

static float** alloc_table(int n)
{
    float** a = static_cast<float**>(std::malloc(n * sizeof(float*)));
    a[0] = static_cast<float*>(std::malloc(n * n * sizeof(float)));

    for (int i = 1; i < n; i++) {
        a[i] = a[0] + i * n;
    }

    return a;
}

static void free_table(float** a)
{
    if (!a) return;
    std::free(a[0]);
    std::free(a);
}

int main()
{
    FreqContextView ctx;

    ctx.nt = 512;
    ctx.nh = 128;
    ctx.ns = 64;
    ctx.nsg = 128;
    ctx.nrg = 128;
    ctx.aper = 32;
    ctx.tap = 10;
    ctx.cmp = 0;
    ctx.nfft = 1024;
    ctx.nw = ctx.nfft / 2 + 1;
    ctx.nsam = 32;

    ctx.dt = 0.001f;
    ctx.h0 = 0.0f;
    ctx.dh = 10.0f;
    ctx.s0 = 0.0f;
    ctx.ds = 10.0f;
    ctx.sg0 = 0.0f;
    ctx.dsg = 10.0f;
    ctx.rg0 = 0.0f;
    ctx.drg = 10.0f;
    ctx.sdatum = 100.0f;
    ctx.rdatum = 100.0f;

    ctx.stable = alloc_table(ctx.nsg);
    ctx.rtable = alloc_table(ctx.nrg);

    /*
     * Build simple smooth travel-time tables.
     */
    const float vel = 2000.0f;

    for (int i = 0; i < ctx.nrg; i++) {
        for (int j = 0; j < ctx.nrg; j++) {
            float dx = (i - j) * ctx.drg;
            float r = std::sqrt(ctx.rdatum * ctx.rdatum + dx * dx);
            ctx.rtable[i][j] = r / vel;
        }
    }

    for (int i = 0; i < ctx.nsg; i++) {
        for (int j = 0; j < ctx.nsg; j++) {
            float dx = (i - j) * ctx.dsg;
            float r = std::sqrt(ctx.sdatum * ctx.sdatum + dx * dx);
            ctx.stable[i][j] = r / vel;
        }
    }

    const int iw = 80;
    const int is = 10;

    std::vector<Complexf> x(ctx.nh);
    std::vector<Complexf> y_direct(ctx.nh);
    std::vector<Complexf> y_bf(ctx.nh);

    for (int i = 0; i < ctx.nh; i++) {
        float phase = 0.05f * i;
        x[i] = Complexf(std::cos(phase), std::sin(phase));
    }

    /*
     * Direct receiver-side gather.
     */
    KernelEvaluator kernel(&ctx);

    for (int ih = 0; ih < ctx.nh; ih++) {
        y_direct[ih] = Complexf(0.0f, 0.0f);

        int left = (ih - ctx.aper < 0) ? 0 : ih - ctx.aper;
        int right = (ih + ctx.aper > ctx.nh - 1) ? ctx.nh - 1 : ih + ctx.aper;

        for (int ic = left; ic <= right; ic++) {
            Complexf K = kernel.eval_receiver(iw, is, ih, ic);
            y_direct[ih] += K * x[ic];
        }
    }

    ButterflyOptions opt;
    opt.order = 8;
    opt.leaf_size = 16;
    opt.max_level = 0;
    opt.use_aperture = 1;
    opt.verbose = 1;

    butterfly_apply_receiver_gather_iw(&ctx,
                                       iw,
                                       is,
                                       x.data(),
                                       y_bf.data(),
                                       opt);

    float err = rel_l2_error(y_bf, y_direct);

    std::printf("receiver gather direct vs butterfly rel_l2_error = %.8e\n", err);

    if (err > 1e-3f) {
        std::printf("ERROR: butterfly result is not accurate enough.\n");
        free_table(ctx.stable);
        free_table(ctx.rtable);
        return 1;
    }

    free_table(ctx.stable);
    free_table(ctx.rtable);

    std::printf("PASS\n");
    return 0;
}
```

---

## 22. 如何接入当前频率域外推代码

当前主程序中原来有：

```c
for (int iw = 0; iw < ctx.nw; iw++) {
    apply_receiver_direct_operator_freq_iw(&ctx, Uin, Utmp, iw);
}

for (int iw = 0; iw < ctx.nw; iw++) {
    apply_source_direct_operator_freq_iw(&ctx, Utmp, Uout, iw);
}
```

后续可以改成：

```c
for (int iw = 0; iw < ctx.nw; iw++) {
    if (ctx.use_butterfly) {
        apply_receiver_butterfly_operator_freq_iw(&ctx, Uin, Utmp, iw);
    } else {
        apply_receiver_direct_operator_freq_iw(&ctx, Uin, Utmp, iw);
    }
}

for (int iw = 0; iw < ctx.nw; iw++) {
    if (ctx.use_butterfly && ctx.cmp == 0) {
        apply_source_butterfly_operator_freq_iw(&ctx, Utmp, Uout, iw);
    } else {
        apply_source_direct_operator_freq_iw(&ctx, Utmp, Uout, iw);
    }
}
```

建议先只支持：

```text
receiver butterfly: yes
source cmp=0 butterfly: yes
source cmp=1 butterfly: no, fallback direct
```

---

## 23. `apply_receiver_butterfly_operator_freq_iw()` 示例

```cpp
void apply_receiver_butterfly_operator_freq_iw(
    const FreqContextView* ctx,
    const sebf::Complexf* Uin,
    sebf::Complexf* Utmp,
    int iw)
{
    sebf::ButterflyOptions opt;
    opt.order = 8;
    opt.leaf_size = 32;
    opt.max_level = 0;
    opt.use_aperture = 1;
    opt.verbose = 0;

    std::vector<sebf::Complexf> x(ctx->nh);
    std::vector<sebf::Complexf> y(ctx->nh);

    for (int is = 0; is < ctx->ns; is++) {
        for (int ic = 0; ic < ctx->nh; ic++) {
            size_t idx = (static_cast<size_t>(is) * ctx->nh + ic) * ctx->nw + iw;
            x[ic] = Uin[idx];
        }

        sebf::butterfly_apply_receiver_gather_iw(ctx,
                                                 iw,
                                                 is,
                                                 x.data(),
                                                 y.data(),
                                                 opt);

        for (int ih = 0; ih < ctx->nh; ih++) {
            size_t idx = (static_cast<size_t>(is) * ctx->nh + ih) * ctx->nw + iw;
            Utmp[idx] += y[ih];
        }
    }
}
```

---

## 24. `apply_source_butterfly_operator_freq_iw()` 示例，仅 cmp=0

```cpp
void apply_source_butterfly_operator_freq_iw(
    const FreqContextView* ctx,
    const sebf::Complexf* Utmp,
    sebf::Complexf* Uout,
    int iw)
{
    if (ctx->cmp != 0) {
        /*
         * Source cmp=1 is irregular and should fall back to direct summation.
         */
        return;
    }

    sebf::ButterflyOptions opt;
    opt.order = 8;
    opt.leaf_size = 32;
    opt.max_level = 0;
    opt.use_aperture = 1;
    opt.verbose = 0;

    std::vector<sebf::Complexf> x(ctx->ns);
    std::vector<sebf::Complexf> y(ctx->ns);

    for (int ih = 0; ih < ctx->nh; ih++) {
        for (int ic = 0; ic < ctx->ns; ic++) {
            size_t idx = (static_cast<size_t>(ic) * ctx->nh + ih) * ctx->nw + iw;
            x[ic] = Utmp[idx];
        }

        sebf::butterfly_apply_source_cmp0_gather_iw(ctx,
                                                    iw,
                                                    ih,
                                                    x.data(),
                                                    y.data(),
                                                    opt);

        for (int is = 0; is < ctx->ns; is++) {
            size_t idx = (static_cast<size_t>(is) * ctx->nh + ih) * ctx->nw + iw;
            Uout[idx] += y[is];
        }
    }
}
```

---

## 25. `cmp=1` 的处理建议

`cmp=1` source pass 中存在：

\[
hh = \frac{r-i_c\Delta s}{\Delta h}+0.5.
\]

也就是说，输入不是规则的：

\[
U_{\rm tmp}(i_c,i_h,\omega),
\]

而是：

\[
U_{\rm tmp}(i_c,h_h(i_c),\omega).
\]

因此 `cmp=1` 的 source-side 不是标准规则一维矩阵乘法。建议分阶段处理：

```text
Phase 1:
receiver pass 使用 butterfly

Phase 2:
source cmp=0 使用 butterfly

Phase 3:
source cmp=1 保持 direct

Phase 4:
单独设计 irregular butterfly 接口处理 cmp=1
```

不要一开始就强行把 `cmp=1` 改成规则蝶形算法。

---

## 26. 精度参数建议

初始建议：

```text
order = 8
leaf_size = 32
```

如果误差过大：

```text
order = 10 或 12
leaf_size = 16
```

如果速度不够：

```text
order = 6
leaf_size = 64
```

验证时先使用小模型，例如：

```text
nh = 128
ns = 64
nw 只选若干频率
aperture = 32
```

通过后再上大模型。

---

## 27. 关键提醒

必须在指南和代码注释中反复强调：

```text
MATLAB 版本中的 ctheta * exp(1j * omega * tt) 只是测试核。
当前 C++ 实现不能使用这个核。
当前 C++ 实现必须使用 W(r,s) * H_fast(omega,tau)。
```

如果后续另一个 AI 照搬 MATLAB 核函数，结果很可能：

1. 与频率域直接求和不一致；
2. 与时间域外推不一致；
3. 相位方向错误；
4. 振幅权重错误；
5. 无法用于实际 Kirchhoff 外推。

---

## 28. 最终验收标准

该 SEBF 实现完成后，必须满足：

### 28.1 receiver-side 单频测试

\[
\frac{
\|U_{\rm bf}-U_{\rm direct}\|_2
}{
\|U_{\rm direct}\|_2
}
< 10^{-3}
\]

优先达到：

```text
1e-4 或更低
```

### 28.2 source-side cmp=0 单频测试

同样要求：

```text
relative L2 error < 1e-3
```

### 28.3 多频测试

选取若干频率，例如：

```text
iw = 10, 20, 40, 80, 120
```

逐频比较。

### 28.4 完整外推测试

将 butterfly 接入主程序后，比较：

```text
use_butterfly=0
use_butterfly=1
```

输出结果的相对误差。

推荐阈值：

```text
relative L2 error < 1e-3
correlation > 0.999
```

如果误差过大，优先检查：

1. kernel 是否仍然是 `W * H_fast`；
2. 相位符号是否为 FFTW 对应的负号；
3. aperture 是否一致；
4. taper 是否一致；
5. direct 和 butterfly 是否使用了同一个 `KernelEvaluator`；
6. source/receiver 的索引是否颠倒；
7. `cmp=1` 是否误用了规则 butterfly。

---

## 29. 最终结论

按照本指南实现后，C++/SEBF 蝶形算法应当首先与当前频率域 direct summation 结果一致。

只要 direct summation 本身能与当前快速频率域外推代码一致，则 butterfly 接入后也应与其基本一致。

如果当前快速频率域代码和原始时间域外推代码已经基本一致，那么 butterfly 版本也可以获得与时间域外推基本一致的结果。

但是，这个一致性的前提是：

\[
\boxed{
\text{butterfly 使用的核函数必须是当前 C/C++ 外推核，而不是 MATLAB 示例核。}
}
\]

最终应使用：

\[
\boxed{
K_{\rm cpp}(r,s,\omega)
=
W(r,s)H_{\rm fast}(\omega,t(s,r)).
}
\]

不要使用：

\[
\boxed{
K_{\rm matlab}(r,s,\omega)
=
c_\theta(r,s)e^{i\omega t(s,r)}.
}
\]
