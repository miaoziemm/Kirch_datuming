# SEBF 相位残差蝶形算法与频率域波场外推完整实现指南（可独立实现版）

本文档用于指导开发者在不依赖外部资料的情况下，实现一套可接入当前频率域 Kirchhoff 型波场外推代码的 C++ 相位残差蝶形算法库。所有代码建议放在项目的 `SEBF/` 目录下。

本文档合并并修正了前面所有讨论中的要点，尤其包括：

1. 当前频率域波场外推代码的真实计算结构；
2. receiver-side、source-side `cmp=0`、source-side `cmp=1` 三类算子的正确定义；
3. 当前完整外推核函数 \(K=W\cdot AA\cdot H_{\rm fast}\)；
4. MATLAB 蝶形算法中必须保留的相位残差结构；
5. `cmp=1` 的正确理解：固定绝对检波点 \(r=x_s+h\) 后的 local common-receiver source-axis operator；
6. Chebyshev `order < leaf_size` 的必要约束；
7. 已完成的 `order < leaf_size` 验证结论：当前完整核在低频可用，中高频误差较大；
8. 推荐工程策略：低频 butterfly，高频 direct fallback；
9. 完整 SEBF 文件结构、核心接口、代码级伪代码和测试流程。

---

## 0. 当前验证结论与工程定位

### 0.1 已验证内容

已经用脚本测试了三类当前完整外推核：

\[
K_{\rm cpp}(r,s,\omega)
=
W(r,s)\operatorname{AA}(r,s,\omega)H_{\rm fast}(\omega,t(s,r)).
\]

测试对象包括：

1. receiver-side；
2. source-side `cmp=0`；
3. source-side `cmp=1` 固定绝对检波点 \(r=x_s+h\) 后的局部 source-axis gather。

测试使用合法参数：

\[
\boxed{order < leaf\_size}
\]

例如：

```text
(order, leaf_size) = (8,16), (16,32), (32,64), (60,64)
```

### 0.2 测试结论

低频误差可以较小，但中高频误差明显增大。例如 `order=60, leaf_size=64` 时：

```text
iw=1   error ≈ 1e-3
iw=5   error ≈ 1e-3
iw=10  error ≈ 1e-3
iw=20  error ≈ 3e-3
iw=40  error ≈ 1e-2
iw=80  error ≈ 6e-1
```

因此，当前指南不能宣称“完整核全频率 butterfly 已验证通过”。正确工程定位是：

\[
\boxed{
\text{低频使用 butterfly，高频使用 direct fallback。}
}
\]

初始建议：

```text
bf_order     = 60
bf_leaf_size = 64
bf_iw_max    = 20
```

后续再通过更多模型和误差测试调整 `bf_iw_max`。

---

## 1. 当前频率域外推主流程

当前代码的主流程可抽象为：

```text
read time-domain traces
        ↓
reverse trace in time
        ↓
FFT all traces → Uin(is, ih, iw)
        ↓
receiver-side frequency-domain summation
        ↓
source-side frequency-domain summation
        ↓
IFFT
        ↓
reverse trace in time
        ↓
write output
```

其中 butterfly 只能替换两类空间求和：

```text
receiver-side summation
source-side summation
```

不能改变：

- 文件输入输出；
- FFT/IFFT；
- trace reverse；
- `nfft` 策略；
- `stable/rtable` 读取；
- 走时索引；
- aperture；
- taper；
- antialias；
- 当前 \(H_{\rm fast}\)；
- `cmp=0/1` 的几何关系。

---

## 2. 当前快速版约束

### 2.1 `nfft`

当前快速版使用：

\[
N_{\rm fft}
=
\operatorname{nextpow2}(n_t+p_{\max}+2),
\]

而不是严格双程版本：

\[
\operatorname{nextpow2}(n_t+2p_{\max}+2).
\]

SEBF 默认必须和当前快速版一致。

### 2.2 无中间 IFFT + FFT

receiver pass 后，`Utmp` 直接进入 source pass。除非用户指定 `interm=` 输出中间结果，不应主动增加中间 IFFT/FFT。

### 2.3 无边界修正项

当前快速版没有严格时间域等价的边界项：

\[
E(\omega,\tau)x(0).
\]

SEBF 默认也不加入该项。

---

## 3. 三类外推算子

### 3.1 Receiver-side

固定频率 \(i_\omega\) 和固定 shot \(i_s\)，receiver-side 是一维 receiver-axis 算子：

\[
U_{\rm tmp}(i_s,i_h,\omega)
=
\sum_{i_c}
K_r(i_s,i_h,i_c,\omega)
U_{\rm in}(i_s,i_c,\omega).
\]

局部一维向量：

```text
input  : x[ic] = Uin(is, ic, iw)
output : y[ih] = Utmp(is, ih, iw)
```

### 3.2 Source-side `cmp=0`

固定频率 \(i_\omega\) 和固定 receiver index \(i_h\)，source-side 是一维 source-axis 算子：

\[
U_{\rm out}(i_s,i_h,\omega)
=
\sum_{i_c}
K_s(i_h,i_s,i_c,\omega)
U_{\rm tmp}(i_c,i_h,\omega).
\]

局部一维向量：

```text
input  : x[ic] = Utmp(ic, ih, iw)
output : y[is] = Uout(is, ih, iw)
```

### 3.3 Source-side `cmp=1`

`cmp=1` 表示 receiver 轴为相对 offset：

\[
x_r=x_s+h.
\]

source-side 计算时，按绝对检波点坐标：

\[
r=x_s+h
\]

重排成 common-receiver gather。

固定 \(r\) 后：

\[
i_h(i_s)
=
\operatorname{round}
\left(
\frac{r-i_s\Delta s}{\Delta h}
\right),
\]

\[
h_h(i_c)
=
\operatorname{round}
\left(
\frac{r-i_c\Delta s}{\Delta h}
\right).
\]

局部 source-axis 算子为：

\[
U_{\rm out}(i_s,i_h(i_s),\omega)
=
\sum_{i_c}
K_{cmp1}(r,i_s,i_c,\omega)
U_{\rm tmp}(i_c,h_h(i_c),\omega).
\]

因此 `cmp=1` 不是任意 irregular gather，而是：

\[
\boxed{
\text{固定绝对检波点 }r\text{ 后的 local common-receiver source-axis operator。}
}
\]

它也可以使用局部一维 butterfly。

---

## 4. 当前完整核函数

### 4.1 总体形式

\[
\boxed{
K(r,s,\omega)
=
W(r,s)\operatorname{AA}(r,s,\omega)H_{\rm fast}(\omega,t(s,r)).
}
\]

direct summation 和 butterfly summation 必须使用同一个 `KernelEvaluator`。

---

### 4.2 \(H_{\rm fast}\)

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
\omega_d=\frac{2\pi i_\omega}{N_{\rm fft}},
\]

\[
q=\frac{\tau}{\Delta t}, \qquad
m=\lceil q\rceil,\qquad
\delta=m-q.
\]

\[
G_k(\tau)=
\sqrt{
\left(
\frac{\tau+k\Delta t}{\tau}
\right)^2-1
}.
\]

\[
F_0=G_1-G_0,
\]

\[
F_k=G_{k+1}-2G_k+G_{k-1}.
\]

注意：FFTW 正变换对应负相位：

\[
e^{-i\omega t}.
\]

不能使用 MATLAB 示例中的：

\[
e^{+i\omega t}.
\]

---

### 4.3 Receiver-side 核

\[
K_r
=
w_{\rm taper}
\frac{\Delta h}{\pi}
\frac{rdatum\cdot\tau_r}
{rdatum^2+(i_c-i_h)^2\Delta h^2}
\operatorname{AA}_r
H_{\rm fast}(\omega,\tau_r).
\]

\[
\tau_r=rtable[cc][c].
\]

索引：

\[
c=
\operatorname{round}
\left(
\frac{
(cmp ? s_0+i_s\Delta s : 0)
+h_0+i_h\Delta h-rg_0
}{drg}
\right),
\]

\[
cc=
\operatorname{round}
\left(
\frac{
(cmp ? s_0+i_s\Delta s : 0)
+h_0+i_c\Delta h-rg_0
}{drg}
\right).
\]

---

### 4.4 Source `cmp=0` 核

\[
K_s
=
w_{\rm taper}
\frac{\Delta s}{\pi}
\frac{sdatum\cdot\tau_s}
{sdatum^2+(i_c-i_s)^2\Delta s^2}
\operatorname{AA}_s
H_{\rm fast}(\omega,\tau_s).
\]

\[
\tau_s=stable[cc][c].
\]

\[
c=
\operatorname{round}
\left(
\frac{s_0+i_s\Delta s-sg_0}{dsg}
\right),
\]

\[
cc=
\operatorname{round}
\left(
\frac{s_0+i_c\Delta s-sg_0}{dsg}
\right).
\]

---

### 4.5 Source `cmp=1` 核

固定 \(ir,r\) 后：

\[
i_s(q)=sleft+q\cdot jump,
\]

\[
i_c(p)=sleft+p\cdot jump.
\]

输出 offset：

\[
i_h(q)=
\operatorname{round}
\left(
\frac{r-i_s(q)\Delta s}{\Delta h}
\right).
\]

输入 offset：

\[
h_h(p)=
\operatorname{round}
\left(
\frac{r-i_c(p)\Delta s}{\Delta h}
\right).
\]

若 \(i_c\) 超出当前 \(i_s\) 的 aperture，则 \(K=0\)。否则：

\[
K_{cmp1}
=
w_{\rm taper}
\frac{\Delta s}{\pi}
\frac{sdatum\cdot\tau_s}
{sdatum^2+(i_c-i_s)^2\Delta s^2}
\operatorname{AA}_{cmp1}
H_{\rm fast}(\omega,\tau_s).
\]

---

## 5. 反假频权重

实现函数：

```cpp
float antialias_weight_iw(int iw, float dtau)
```

逻辑必须与当前代码一致：

```cpp
if (antialias <= 0) return 1;
if (iw <= 0) return 1;

dtau = fabs(dtau) * antialias;

freq = iw / (nfft * dt);
fnyq = 0.5 / dt;
fmax = 0.5 / dtau;

if (fmax >= fnyq) return 1;
if (fmax <= 0) return 0;

fpass = 0.8 * fmax;

if (freq <= fpass) return 1;
if (freq >= fmax) return 0;

x = (freq - fpass) / (fmax - fpass);
return 0.5 * (1 + cos(pi*x));
```

---

## 6. MATLAB 中必须保留的 butterfly 要素

用户提供的 MATLAB 代码中包含以下关键要素，SEBF 必须保留：

1. binary source/target trees；
2. Chebyshev nodes；
3. barycentric weights；
4. Lagrange basis；
5. Lagrange matrix；
6. 相位补偿初始化；
7. 第一段向上/互补递推；
8. switch；
9. 第二段向下/互补递推；
10. 相位补偿终止重构；
11. direct vs butterfly 对比测试。

但不能照搬 MATLAB 核：

\[
K_{\rm matlab}=c_\theta e^{+i\omega t}.
\]

当前必须使用：

\[
K_{\rm cpp}=W\cdot AA\cdot H_{\rm fast}.
\]

---

## 7. 相位残差定义

由于当前核不是单纯指数核，所以相位残差只从 \(H_{\rm fast}\) 中提取。

\[
H_{\rm fast}(\omega,\tau)=|H_{\rm fast}(\omega,\tau)|e^{i\Phi(\omega,\tau)}.
\]

相位差因子：

\[
P(a,b)=e^{i[\Phi(a)-\Phi(b)]}.
\]

代码中用复数比值，避免 unwrap：

```cpp
phase_ratio(Ha, Hb) = Ha * conj(Hb) / |Ha * conj(Hb)|
```

其中 `Ha`、`Hb` 是对应点对上的 \(H_{\rm fast}\)。

不要把 \(W\) 或 \(AA\) 放进相位残差，因为它们是实权重，且可能因符号或截断产生不必要跳变。

---

## 8. 相位残差 butterfly 公式

### 8.1 初始化

\[
\sigma_t^{A_0,B}
=
\sum_{s_j\in B}
L_t^B(s_j)
P\big((r_0^{A_0},s_j),(r_0^{A_0},s_t^B)\big)
x_j.
\]

### 8.2 第一段递推

\[
\sigma_t^{A,B}
=
\sum_{B_c\subset B}
\sum_{t'}
L_t^B(s_{t'}^{B_c})
P\big((r_0^A,s_{t'}^{B_c}),(r_0^A,s_t^B)\big)
\sigma_{t'}^{A_p,B_c}.
\]

### 8.3 Switch

\[
\gamma_t^{A,B}
=
\sum_{t'}
K(r_t^A,s_{t'}^B,\omega)
\sigma_{t'}^{A,B}.
\]

这里必须用完整核：

\[
K=W\cdot AA\cdot H_{\rm fast}.
\]

### 8.4 第二段递推

\[
\gamma_t^{A,B}
=
\sum_{B_c\subset B}
\sum_{t'}
L_{t'}^{A_p}(r_t^A)
P\big((r_t^A,s_0^{B_c}),(r_{t'}^{A_p},s_0^{B_c})\big)
\gamma_{t'}^{A_p,B_c}.
\]

### 8.5 终止重构

\[
y(r_i)
=
\sum_t
L_t^A(r_i)
P\big((r_i,s_0^{B_{\rm root}}),(r_t^A,s_0^{B_{\rm root}})\big)
\gamma_t^{A,B_{\rm root}}.
\]

---

## 9. `order < leaf_size` 约束

`order` 是每个 box 中使用的 Chebyshev 插值点数。`leaf_size` 是叶子 box 中包含的原始网格点数。

必须满足：

\[
\boxed{
order < leaf\_size
}
\]

否则插值压缩没有意义。

推荐：

```text
order = 16, leaf_size = 32
order = 24, leaf_size = 48
order = 32, leaf_size = 64
order = 60, leaf_size = 64  # 仅用于高精度验证，压缩率很低
```

禁止把：

```text
order = 64, leaf_size = 16
```

作为有效 butterfly 验证。

---

## 10. 推荐 SEBF 文件结构

```text
SEBF/
├── include/
│   ├── sebf_complex.h
│   ├── sebf_context.h
│   ├── sebf_kernel.h
│   ├── sebf_direct.h
│   ├── sebf_chebyshev.h
│   ├── sebf_tree.h
│   ├── sebf_butterfly.h
│   ├── sebf_wave_extrap.h
│   └── sebf_test_utils.h
├── src/
│   ├── sebf_kernel.cpp
│   ├── sebf_direct.cpp
│   ├── sebf_chebyshev.cpp
│   ├── sebf_tree.cpp
│   ├── sebf_butterfly.cpp
│   └── sebf_wave_extrap.cpp
└── tests/
    ├── test_matlab_reference.cpp
    ├── test_receiver_single_frequency.cpp
    ├── test_source_cmp0_single_frequency.cpp
    ├── test_source_cmp1_single_gather.cpp
    └── test_full_pipeline.cpp
```

---

## 11. `sebf_complex.h`

```cpp
#ifndef SEBF_COMPLEX_H
#define SEBF_COMPLEX_H

#include <cmath>

namespace sebf {

struct Complexf {
    float r;
    float i;

    Complexf() : r(0.0f), i(0.0f) {}
    Complexf(float rr, float ii) : r(rr), i(ii) {}

    Complexf& operator+=(const Complexf& b) {
        r += b.r;
        i += b.i;
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
    return Complexf(a.r*b.r - a.i*b.i,
                    a.r*b.i + a.i*b.r);
}

inline Complexf operator*(float a, const Complexf& b) {
    return Complexf(a*b.r, a*b.i);
}

inline Complexf operator*(const Complexf& a, float b) {
    return Complexf(a.r*b, a.i*b);
}

inline Complexf conj(const Complexf& a) {
    return Complexf(a.r, -a.i);
}

inline float abs2(const Complexf& a) {
    return a.r*a.r + a.i*a.i;
}

inline Complexf normalized_phase(const Complexf& z) {
    float mag = std::sqrt(z.r*z.r + z.i*z.i);
    if (mag < 1e-20f) return Complexf(1.0f, 0.0f);
    return Complexf(z.r/mag, z.i/mag);
}

inline Complexf phase_ratio(const Complexf& a, const Complexf& b) {
    return normalized_phase(a * conj(b));
}

} // namespace sebf

#endif
```

---

## 12. `sebf_context.h`

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
    int verb;

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

    float antialias;

    float** stable;
    float** rtable;
};

struct ButterflyOptions {
    int order;
    int leaf_size;
    int max_level;
    int use_butterfly;
    int verbose;
    int use_phase_residual;
    int bf_iw_max;
};

struct Cmp1GatherMap {
    int ir;
    float r;

    int sleft;
    int sright;
    int jump;
    int nloc;

    int* isrc;
    int* ih;
};

} // namespace sebf

#endif
```

---

## 13. 核函数接口 `sebf_kernel.h`

```cpp
#ifndef SEBF_KERNEL_H
#define SEBF_KERNEL_H

#include "sebf_complex.h"
#include "sebf_context.h"

namespace sebf {

enum class KernelType {
    Receiver,
    SourceCmp0,
    SourceCmp1Gather
};

class KernelEvaluator {
public:
    explicit KernelEvaluator(const FreqContextView* ctx);

    Complexf eval_receiver(int iw, int is, int ih, int ic) const;
    Complexf eval_source_cmp0(int iw, int ih, int is, int ic) const;
    Complexf eval_source_cmp1_gather(int iw,
                                     const Cmp1GatherMap* gm,
                                     int q_out,
                                     int p_in) const;

    Complexf eval_receiver_H(int iw, int is, int ih, int ic) const;
    Complexf eval_source_cmp0_H(int iw, int ih, int is, int ic) const;
    Complexf eval_source_cmp1_gather_H(int iw,
                                       const Cmp1GatherMap* gm,
                                       int q_out,
                                       int p_in) const;

    Complexf eval_regular(KernelType type,
                          int iw,
                          int gather_id,
                          int iout,
                          int iin) const;

    Complexf eval_regular_H(KernelType type,
                            int iw,
                            int gather_id,
                            int iout,
                            int iin) const;

    Complexf phase_factor_regular(KernelType type,
                                  int iw,
                                  int gather_id,
                                  int out_a,
                                  int in_a,
                                  int out_b,
                                  int in_b) const;

private:
    const FreqContextView* ctx_;

    float taper_weight(int left, int center, int right, int tap, int stride) const;
    float table_slope_first_index(float** table, int n, int row, int col, float dcoord) const;
    float antialias_weight_iw(int iw, float dtau) const;

    void build_filter_F(float tau, float* F) const;
    Complexf build_H_fast(float tau, int iw) const;
};

int build_cmp1_gather_map(const FreqContextView* ctx, int ir, Cmp1GatherMap* gm);
void free_cmp1_gather_map(Cmp1GatherMap* gm);

} // namespace sebf

#endif
```

Implementation notes:

- `eval_*()` returns full \(K=WAAH\).
- `eval_*_H()` returns only \(H_{\rm fast}\).
- `phase_factor_*()` uses `phase_ratio(Ha,Hb)`.

---

## 14. Direct summation

Implement direct reference functions:

```cpp
void direct_apply_regular_1d(
    const FreqContextView* ctx,
    KernelType type,
    int iw,
    int gather_id,
    int n_in,
    int n_out,
    const Complexf* x,
    Complexf* y);

void direct_apply_receiver_iw(...);
void direct_apply_source_cmp0_iw(...);
void direct_apply_source_cmp1_gather_iw(...);
```

`direct_apply_regular_1d`:

```cpp
for iout:
    y[iout] = 0
    for iin:
        K = kernel.eval_regular(type, iw, gather_id, iout, iin)
        y[iout] += K * x[iin]
```

---

## 15. Chebyshev utilities

Implement:

```cpp
void chebyshev_nodes_1d(int order, float a, float b, vector<float>& nodes);
void barycentric_weights_1d(const vector<float>& nodes, vector<float>& weights);
void lagrange_all_1d(const vector<float>& nodes,
                     const vector<float>& weights,
                     float x,
                     vector<float>& L);
void build_lagrange_matrix_1d(...);
```

Barycentric Lagrange basis:

\[
L_j(x)
=
\frac{w_j/(x-x_j)}
{\sum_k w_k/(x-x_k)}.
\]

If \(x=x_j\):

\[
L_j=1,\quad L_{k\ne j}=0.
\]

---

## 16. Tree structure

Use snapped Chebyshev nodes to call the existing discrete `KernelEvaluator`:

```cpp
struct Box1D {
    int id;
    int level;
    int begin;
    int end;
    int parent;
    int left;
    int right;

    float xmin;
    float xmax;
    float center_coord;
    int center_index;

    vector<int> cheb_index;
    vector<float> cheb_coord;
    vector<float> bary_w;
};
```

Rules:

1. generate continuous Chebyshev nodes;
2. snap each node to nearest grid index;
3. avoid duplicate indices;
4. set `cheb_coord` from snapped indices.

If `box_size <= order`, use direct fallback for that box or reduce local order. For the first implementation, enforce:

```text
order < leaf_size
```

globally.

---

## 17. Butterfly coefficient storage

```cpp
using CoeffTable = unordered_map<uint64_t, vector<Complexf>>;

uint64_t pair_key(int aid, int bid)
{
    return (uint64_t(uint32_t(aid)) << 32) | uint32_t(bid);
}
```

- `sigma[key]`: source-side phase-compensated coefficient;
- `gamma[key]`: target-side phase-compensated coefficient.

---

## 18. Butterfly algorithm steps

### 18.1 Build trees

```cpp
source_tree.build_uniform(n_in,  xin0,  dxin,  leaf_size, order);
target_tree.build_uniform(n_out, xout0, dxout, leaf_size, order);
```

If depths differ, first implementation should fallback to direct.

### 18.2 Initialize sigma

\[
\sigma_t^{A_0,B}
=
\sum_{s_j\in B}
L_t^B(s_j)
P((r_0^{A_0},s_j),(r_0^{A_0},s_t^B))
x_j.
\]

### 18.3 First recursion

\[
\sigma_t^{A,B}
=
\sum_{B_c}
\sum_{t'}
L_t^B(s_{t'}^{B_c})
P((r_0^A,s_{t'}^{B_c}),(r_0^A,s_t^B))
\sigma_{t'}^{A_p,B_c}.
\]

### 18.4 Switch

\[
\gamma_t^{A,B}
=
\sum_{t'}
K(r_t^A,s_{t'}^B,\omega)
\sigma_{t'}^{A,B}.
\]

### 18.5 Second recursion

\[
\gamma_t^{A,B}
=
\sum_{B_c}
\sum_{t'}
L_{t'}^{A_p}(r_t^A)
P((r_t^A,s_0^{B_c}),(r_{t'}^{A_p},s_0^{B_c}))
\gamma_{t'}^{A_p,B_c}.
\]

### 18.6 Termination

\[
y(r_i)
=
\sum_t
L_t^A(r_i)
P((r_i,s_0^{B_{\rm root}}),(r_t^A,s_0^{B_{\rm root}}))
\gamma_t^{A,B_{\rm root}}.
\]

---

## 19. `cmp=1` gather map

Implement:

```cpp
int build_cmp1_gather_map(const FreqContextView* ctx, int ir, Cmp1GatherMap* gm);
```

It must reproduce current code:

```cpp
s = fabs((ns-1)*ds)
h = fabs((nh-1)*dh)

jump = fabs(ds) >= fabs(dh) ? 1 : round(dh/ds)
dr   = fabs(ds) >= fabs(dh) ? fabs(dh) : fabs(ds)

r = ir*dr
  + (ds<=0 ? -1 : 0)*s
  + (dh<=0 ? -1 : 0)*h

sleft  = ...
sright = ...

isrc[p] = sleft + p*jump
ih[p]   = round((r - isrc[p]*ds)/dh)
```

For `cmp=1` local butterfly:

```cpp
x[p] = Utmp(isrc[p], ih[p], iw)
y[q] -> Uout(isrc[q], ih[q], iw)
```

---

## 20. Wave extrapolation wrappers

### 20.1 Receiver wrapper

For each `(iw,is)`:

```cpp
x[ic] = Uin(is, ic, iw)
butterfly_apply_regular_1d(...)
Utmp(is, ih, iw) += y[ih]
```

### 20.2 Source `cmp=0`

For each `(iw,ih)`:

```cpp
x[ic] = Utmp(ic, ih, iw)
butterfly_apply_regular_1d(...)
Uout(is, ih, iw) += y[is]
```

### 20.3 Source `cmp=1`

For each `(iw,ir)`:

```cpp
build_cmp1_gather_map(...)
x[p] = Utmp(isrc[p], ih[p], iw)
butterfly_apply_cmp1_regular_1d(...)
Uout(isrc[q], ih[q], iw) += y[q]
```

---

## 21. Frequency fallback strategy

Because current complete kernel is not validated at all frequencies under legal `order < leaf_size`, implement:

```cpp
if (use_butterfly && iw <= bf_iw_max) {
    use butterfly;
} else {
    use direct;
}
```

Initial recommendation:

```text
bf_order     = 60
bf_leaf_size = 64
bf_iw_max    = 20
```

Then adjust by tests.

---

## 22. Testing requirements

### 22.1 MATLAB reference kernel

First implement a test for MATLAB example kernel:

\[
K_{\rm matlab}=c_\theta e^{+i\omega t}.
\]

This verifies the butterfly recursion itself.

### 22.2 Current complete kernels

Test:

1. receiver single frequency;
2. source `cmp=0` single frequency;
3. source `cmp=1` single gather.

Use `order < leaf_size`.

Expected:

- low frequencies: small error;
- high frequencies: fallback direct if error is too large.

### 22.3 Full pipeline

Run:

```bash
use_butterfly=0
use_butterfly=1
```

Compare output only over frequencies where butterfly is enabled, or compare full output after high-frequency direct fallback.

---

## 23. Acceptance criteria

For frequencies using butterfly:

```text
relative L2 error < 1e-3 to 1e-2
```

depending on acceptable accuracy.

For full output:

```text
correlation > 0.999
```

and no visible artifacts.

If not satisfied, reduce `bf_iw_max`, increase `order`, decrease `leaf_size`, or use direct fallback for problematic frequencies.

---

## 24. Development order

```text
1. Implement Complexf
2. Implement KernelEvaluator full K and H-only
3. Implement direct summation
4. Implement MATLAB reference butterfly test
5. Implement Chebyshev utilities
6. Implement Tree1D
7. Implement phase-residual butterfly
8. Validate receiver/source0/source1 under order < leaf_size
9. Add bf_iw_max fallback
10. Connect receiver pass
11. Connect source cmp=0 pass
12. Connect source cmp=1 pass
13. Run full pipeline tests
```

---

## 25. Final warning

Do not claim that the current complete kernel is fully validated across all frequencies by butterfly. The correct statement is:

\[
\boxed{
\text{The MATLAB phase-residual butterfly structure is valid.}
}
\]

\[
\boxed{
\text{The current full extrapolation kernel is suitable for low-frequency butterfly under } order<leaf\_size.
}
\]

\[
\boxed{
\text{High frequencies currently require direct fallback or a more advanced } H_{\rm fast}\text{ decomposition.}
}
\]

This is the safest and most accurate implementation path.
