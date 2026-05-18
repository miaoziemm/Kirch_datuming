# SEBF 相位残差蝶形算法与频率域波场外推完整实现指南（已验证版）

本文档用于指导开发者在 `SEBF` 文件夹中实现一套 C++ 相位残差蝶形算法，并将其用于当前频率域 Kirchhoff 型波场外推代码。本文档基于用户提供的 MATLAB 蝶形算法代码重新整理，并在写指南前先用 Python 原型复现了 MATLAB 中的相位残差递推结构，验证 butterfly 结果与 direct summation 结果可以达到较低误差。

---

## 0. 验证结论

在生成本指南前，先写了一个简单 Python 原型，完整复现 MATLAB 代码中的以下要素：

1. `build_binary_tree_1d`
2. Chebyshev nodes
3. `bary_weights_1d`
4. `lagrange_basis_all_1d`
5. `build_lagrange_matrix_1d`
6. direct Kirchhoff summation
7. source-side phase-compensated initialization
8. first-half recursion
9. switch
10. second-half recursion
11. termination
12. final relative L2 error check

使用 MATLAB 示例核：

\[
K_{\rm matlab}(r,s)
=
c_\theta(r,s)\exp(i\omega t(r,s)),
\]

验证结果如下：

```text
n = 256, leaf_n = 16
p = 8   relative L2 error = 7.5346e-02
p = 12  relative L2 error = 5.7448e-03
p = 16  relative L2 error = 2.8014e-04
p = 20  relative L2 error = 1.3111e-05
p = 24  relative L2 error = 7.9568e-07
p = 32  relative L2 error = 4.0922e-09
```

进一步使用更接近 MATLAB 原始规模的测试：

```text
n = 1024, leaf_n = 64, p = 60
relative L2 error = 2.0045e-04
```

因此，MATLAB 中的相位残差 butterfly 递推结构本身是正确的，可以作为 SEBF C++ 实现的核心结构。

---

## 1. 本指南的核心原则

当前 C++ 波场外推核函数不是 MATLAB 示例核。MATLAB 的 butterfly 递推结构可以完整借鉴，但外推核函数必须替换为当前频率域代码中的核：

\[
\boxed{
K_{\rm cpp}(r,s,\omega)
=
W(r,s)
\operatorname{AA}(r,s,\omega)
H_{\rm fast}(\omega,t(s,r)).
}
\]

其中：

- \(W(r,s)\)：当前 Kirchhoff 几何权重与 taper；
- \(\operatorname{AA}\)：当前代码中的 `antialias_weight_iw()`；
- \(H_{\rm fast}\)：当前代码中的 `build_H_iw()`；
- \(t(s,r)=\tau\)：`stable` 或 `rtable` 中的走时。

MATLAB 中的：

\[
c_\theta(r,s)\exp(i\omega t(r,s))
\]

只能用于验证 butterfly 结构，不能用于当前 C++ 外推。

---

## 2. 当前外推代码的三个空间算子

### 2.1 Receiver-side

当前 receiver-side 是固定频率、固定 shot \(i_s\) 的一维 receiver-axis 算子：

\[
U_{\rm tmp}(i_s,i_h,\omega)
=
\sum_{i_c}
K_r(i_s,i_h,i_c,\omega)
U_{\rm in}(i_s,i_c,\omega).
\]

因此对每个 `(iw, is)`，可以调用一次一维 butterfly：

```text
input  axis: ic = 0 ... nh-1
output axis: ih = 0 ... nh-1
```

### 2.2 Source-side, `cmp=0`

`cmp=0` 时，source-side 是固定频率、固定 receiver index \(i_h\) 的一维 source-axis 算子：

\[
U_{\rm out}(i_s,i_h,\omega)
=
\sum_{i_c}
K_s(i_h,i_s,i_c,\omega)
U_{\rm tmp}(i_c,i_h,\omega).
\]

因此对每个 `(iw, ih)`，可以调用一次一维 butterfly：

```text
input  axis: ic = 0 ... ns-1
output axis: is = 0 ... ns-1
```

### 2.3 Source-side, `cmp=1`

`cmp=1` 表示 receiver 轴为相对 offset：

\[
x_r=x_s+h.
\]

source-side 外推时，当前代码按绝对检波点坐标：

\[
r=x_s+h
\]

重排为 common-receiver gather。

固定 \(r\) 后，定义：

\[
i_h(i_s)=
\operatorname{round}
\left(
\frac{r-i_s\Delta s}{\Delta h}
\right),
\]

\[
h_h(i_c)=
\operatorname{round}
\left(
\frac{r-i_c\Delta s}{\Delta h}
\right).
\]

于是：

\[
U_{\rm out}(i_s,i_h(i_s),\omega)
=
\sum_{i_c}
K_s(ir,i_s,i_c,\omega)
U_{\rm tmp}(i_c,h_h(i_c),\omega).
\]

因此 `cmp=1` 不是任意不规则 gather。它是固定绝对检波点 \(r\) 后的 local source-axis operator，可以构造局部一维 butterfly。

---

## 3. 当前 C++ 外推核函数

### 3.1 \(H_{\rm fast}\)

当前快速频率响应为：

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
q=\frac{\tau}{\Delta t},
\qquad
m=\lceil q\rceil,
\qquad
\delta=m-q.
\]

滤波器：

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
F_0=G_1-G_0,
\]

\[
F_k=G_{k+1}-2G_k+G_{k-1}.
\]

注意 FFTW 正变换下，时间延迟对应负相位：

\[
e^{-i\omega t}.
\]

所以当前 C++ 不能使用 MATLAB 中的正相位 \(e^{+i\omega t}\)。

---

### 3.2 Receiver-side 核

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

走时：

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

### 3.3 Source-side `cmp=0` 核

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

走时：

\[
\tau_s=stable[cc][c].
\]

索引：

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

### 3.4 Source-side `cmp=1` 核

固定 \(ir,r\) 后，局部 source-axis index 为 \(p,q\)。

\[
i_s(q)=sleft+q\cdot jump,
\]

\[
i_c(p)=sleft+p\cdot jump.
\]

\[
i_h(q)=
\operatorname{round}
\left(
\frac{r-i_s(q)\Delta s}{\Delta h}
\right),
\]

\[
h_h(p)=
\operatorname{round}
\left(
\frac{r-i_c(p)\Delta s}{\Delta h}
\right).
\]

若 \(i_c\) 超出当前 output source \(i_s\) 的 aperture，则：

\[
K=0.
\]

否则：

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

## 4. MATLAB 相位残差结构必须保留

用户提供的 MATLAB 代码中，相位残差出现在四个位置：

1. 初始化；
2. 第一段递推；
3. 第二段递推；
4. 终止重构。

Switch 阶段使用完整核函数矩阵。

C++ 版本也必须保留这些要素。

---

## 5. 当前 C++ 中的相位函数定义

由于当前核不是简单的：

\[
e^{i\omega t},
\]

而是：

\[
H_{\rm fast}(\omega,\tau),
\]

所以定义：

\[
H_{\rm fast}(\omega,\tau)
=
|H_{\rm fast}(\omega,\tau)|
e^{i\Phi(\omega,\tau)}.
\]

相位为：

\[
\Phi(\omega,\tau)=\arg H_{\rm fast}(\omega,\tau).
\]

代码中不要直接使用 `atan2` 和 unwrap，而使用复数比值：

\[
P(a,b)
=
e^{i[\Phi(a)-\Phi(b)]}
=
\frac{H(a)\overline{H(b)}}{|H(a)\overline{H(b)}|+\epsilon}.
\]

代码：

```cpp
inline Complexf phase_ratio(const Complexf& Ha, const Complexf& Hb)
{
    Complexf z(Ha.r * Hb.r + Ha.i * Hb.i,
               Ha.i * Hb.r - Ha.r * Hb.i);

    float mag = std::sqrt(z.r*z.r + z.i*z.i);

    if (mag < 1e-20f) {
        return Complexf(1.0f, 0.0f);
    }

    return Complexf(z.r/mag, z.i/mag);
}
```

注意：

- 相位残差只从 \(H_{\rm fast}\) 中提取；
- 完整核 \(K=W\cdot AA\cdot H_{\rm fast}\) 仍用于 direct 和 switch；
- \(W\) 与 \(AA\) 是实权重，不建议放进相位残差。

---

## 6. 相位残差 butterfly 公式

以下公式完全对应 MATLAB 代码中的结构，只是把 MATLAB 的 \(e^{i\omega t}\) 换成当前 C++ 的相位因子 \(P(a,b)\)。

---

### 6.1 初始化

MATLAB 形式：

\[
\sigma_t
=
e^{-i\omega t(r_0,s_t)}
\sum_{s_j\in B}
L_t^B(s_j)
e^{i\omega t(r_0,s_j)}
f(s_j).
\]

当前 C++ 应写成：

\[
\boxed{
\sigma_t^{A_0,B}
=
\sum_{s_j\in B}
L_t^B(s_j)
P\big((r_0^{A_0},s_j),(r_0^{A_0},s_t^B)\big)
x_j.
}
\]

---

### 6.2 第一段递推

MATLAB 形式：

\[
\sigma_t^{A,B}
=
e^{-i\omega t(r_0^A,s_t^B)}
\sum_{B_c}
\sum_{t'}
L_t^B(s_{t'}^{B_c})
e^{i\omega t(r_0^A,s_{t'}^{B_c})}
\sigma_{t'}^{A_p,B_c}.
\]

当前 C++：

\[
\boxed{
\sigma_t^{A,B}
=
\sum_{B_c\subset B}
\sum_{t'}
L_t^B(s_{t'}^{B_c})
P\big(
(r_0^A,s_{t'}^{B_c}),
(r_0^A,s_t^B)
\big)
\sigma_{t'}^{A_p,B_c}.
}
\]

---

### 6.3 Switch

MATLAB switch：

\[
\sigma_{\rm sw}(r_t^A)
=
\sum_{t'}
K(r_t^A,s_{t'}^B)
\sigma_{t'}^{A,B}.
\]

当前 C++ 完全相同，但 \(K\) 必须为当前完整核：

\[
\boxed{
\gamma_t^{A,B}
=
\sum_{t'}
K_{\rm cpp}(r_t^A,s_{t'}^B,\omega)
\sigma_{t'}^{A,B}.
}
\]

其中：

\[
K_{\rm cpp}=W\cdot AA\cdot H_{\rm fast}.
\]

---

### 6.4 第二段递推

MATLAB 形式：

\[
\gamma_t^{A,B}
=
\sum_{B_c}
e^{i\omega t(r_t^A,s_0^{B_c})}
\sum_{t'}
L_{t'}^{A_p}(r_t^A)
e^{-i\omega t(r_{t'}^{A_p},s_0^{B_c})}
\gamma_{t'}^{A_p,B_c}.
\]

当前 C++：

\[
\boxed{
\gamma_t^{A,B}
=
\sum_{B_c\subset B}
\sum_{t'}
L_{t'}^{A_p}(r_t^A)
P\big(
(r_t^A,s_0^{B_c}),
(r_{t'}^{A_p},s_0^{B_c})
\big)
\gamma_{t'}^{A_p,B_c}.
}
\]

---

### 6.5 终止重构

MATLAB 形式：

\[
u(r)
=
e^{i\omega t(r,s_0^B)}
\sum_t
L_t^A(r)
e^{-i\omega t(r_t^A,s_0^B)}
\gamma_t.
\]

当前 C++：

\[
\boxed{
y(r_i)
=
\sum_t
L_t^A(r_i)
P\big(
(r_i,s_0^{B_{\rm root}}),
(r_t^A,s_0^{B_{\rm root}})
\big)
\gamma_t^{A,B_{\rm root}}.
}
\]

---

## 7. SEBF 文件结构

建议：

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
    ├── test_matlab_kernel_reference.cpp
    ├── test_receiver_single_frequency.cpp
    ├── test_source_cmp0_single_frequency.cpp
    ├── test_source_cmp1_single_gather.cpp
    └── test_full_pipeline.cpp
```

---

## 8. 必须实现的 MATLAB 同等要素

### 8.1 barycentric 权重

对应 MATLAB：

```matlab
bary_weights_1d.m
```

C++：

```cpp
void barycentric_weights_1d(const std::vector<float>& nodes,
                            std::vector<float>& weights)
{
    int q = nodes.size();
    weights.assign(q, 1.0f);

    for (int j = 0; j < q; j++) {
        float w = 1.0f;
        for (int k = 0; k < q; k++) {
            if (k != j) {
                w /= (nodes[j] - nodes[k]);
            }
        }
        weights[j] = w;
    }
}
```

### 8.2 Lagrange basis

对应 MATLAB：

```matlab
lagrange_basis_all_1d.m
```

C++：

```cpp
void lagrange_all_1d(const std::vector<float>& nodes,
                     const std::vector<float>& weights,
                     float x,
                     std::vector<float>& L)
{
    int q = nodes.size();
    L.assign(q, 0.0f);

    float tol = 1e-7f;

    for (int k = 0; k < q; k++) {
        if (std::fabs(x - nodes[k]) < tol) {
            L[k] = 1.0f;
            return;
        }
    }

    float denom = 0.0f;
    for (int k = 0; k < q; k++) {
        denom += weights[k] / (x - nodes[k]);
    }

    for (int j = 0; j < q; j++) {
        L[j] = (weights[j] / (x - nodes[j])) / denom;
    }
}
```

### 8.3 Lagrange matrix

对应 MATLAB：

```matlab
build_lagrange_matrix_1d.m
```

C++：

```cpp
void build_lagrange_matrix_1d(const std::vector<float>& pvec,
                              const std::vector<float>& nodes,
                              const std::vector<float>& weights,
                              std::vector<float>& Lmat)
{
    int m = pvec.size();
    int q = nodes.size();

    Lmat.assign(m*q, 0.0f);

    std::vector<float> L;

    for (int i = 0; i < m; i++) {
        lagrange_all_1d(nodes, weights, pvec[i], L);
        for (int j = 0; j < q; j++) {
            Lmat[i*q + j] = L[j];
        }
    }
}
```

### 8.4 binary tree

对应 MATLAB：

```matlab
build_binary_tree_1d.m
```

C++ 需要保留：

- level；
- box index；
- i1/i2；
- center `x0`；
- width；
- Chebyshev nodes；
- parent/children。

---

## 9. Tree 中的 Chebyshev 节点

MATLAB 使用连续 Chebyshev 节点：

\[
x_t=x_0+\frac{w}{2}\xi_t.
\]

对于当前 C++ 外推，由于 `rtable/stable` 是离散表，建议第一版使用 snapped Chebyshev nodes：

1. 先生成连续 node；
2. 映射到最近 grid index；
3. 若重复，则选择最近未使用 index；
4. 保存：
   ```cpp
   cheb_index
   cheb_coord
   ```

这样可以直接调用当前离散 `KernelEvaluator`。

若以后希望更接近 MATLAB，可实现 `interp2_linear` 对 `rtable/stable` 做连续查询，但第一版不要求。

---

## 10. `KernelEvaluator` 必须提供的接口

```cpp
class KernelEvaluator {
public:
    Complexf eval_receiver(...);              // full K
    Complexf eval_source_cmp0(...);           // full K
    Complexf eval_source_cmp1_gather(...);    // full K

    Complexf eval_receiver_H(...);            // H_fast only
    Complexf eval_source_cmp0_H(...);         // H_fast only
    Complexf eval_source_cmp1_gather_H(...);  // H_fast only

    Complexf phase_factor_regular(...);
    Complexf phase_factor_cmp1(...);
};
```

完整核用于：

```text
direct summation
switch
```

H-only 用于：

```text
initialization phase residual
recursion-1 phase residual
recursion-2 phase residual
termination phase residual
```

---

## 11. `cmp=1` 的 `Cmp1GatherMap`

必须实现：

```cpp
struct Cmp1GatherMap {
    int ir;
    float r;

    int sleft;
    int sright;
    int jump;
    int nloc;

    std::vector<int> isrc;
    std::vector<int> ih;
};
```

构造方式必须与当前代码一致：

```cpp
build_cmp1_gather_map(ctx, ir, gm);
```

其中：

```cpp
isrc[p] = sleft + p * jump;
ih[p]   = round((r - isrc[p] * ds) / dh);
```

`cmp=1` 的 local butterfly 输入：

```cpp
x[p] = Utmp(isrc[p], ih[p], iw);
```

输出：

```cpp
Uout(isrc[q], ih[q], iw) += y[q];
```

注意：对输入向量而言，`ih[p]` 就是原代码里的 `hh`。

---

## 12. Butterfly 数据结构

使用：

```cpp
using CoeffTable = std::unordered_map<uint64_t, std::vector<Complexf>>;

uint64_t pair_key(int aid, int bid)
{
    return (uint64_t(uint32_t(aid)) << 32) | uint32_t(bid);
}
```

- `sigma[key]`：第一段 source-side phase-compensated coefficients；
- `gamma[key]`：switch 后 target-side phase-compensated coefficients。

---

## 13. 核心算法伪代码

```text
build source tree TS
build target tree TR

initialize sigma at:
    target root
    source leaves

for lv = 1 : floor(L/2):
    recursion-1 with phase residual

switch:
    gamma = full K * sigma

for lv = floor(L/2)+1 : L:
    recursion-2 with phase residual

termination:
    reconstruct y on target leaves
```

---

## 14. 初始化代码级公式

```cpp
sigma_t = sum_j L_t(s_j) * P((r0, s_j), (r0, s_t)) * x_j;
```

其中：

```cpp
P = phase_ratio(H(r0, s_j), H(r0, s_t));
```

---

## 15. 第一段递推代码级公式

```cpp
sigma_curr(t) +=
    L_t^B(s_child_node)
    * P((r0_A, s_child_node), (r0_A, s_parent_node))
    * sigma_prev(t_child);
```

---

## 16. Switch 代码级公式

```cpp
gamma(t_target) +=
    K_full(r_target_node, s_source_node)
    * sigma(t_source);
```

---

## 17. 第二段递推代码级公式

```cpp
gamma_curr(t_child_target) +=
    L_parent_target(t_parent at r_child_node)
    * P((r_child_node, s0_child_source),
        (r_parent_node, s0_child_source))
    * gamma_prev(t_parent);
```

---

## 18. 终止重构代码级公式

```cpp
y(real_target_point) +=
    L_t(real_target_point)
    * P((real_target_point, s0_source_root),
        (target_node, s0_source_root))
    * gamma(t);
```

---

## 19. 测试 1：MATLAB 示例核 reference

必须先实现一个 reference test，完全复现 MATLAB 核：

\[
K=c_\theta e^{i\omega t}.
\]

测试参数：

```text
v = 3000
omega = 15 * 2*pi
n = 256
dx = 5
dz = 100
leaf_n = 16
p = 20 or 24
```

期望：

```text
p = 20: relative L2 error about 1e-5
p = 24: relative L2 error about 1e-6
```

该测试用于证明 butterfly 结构和相位残差递推写对。

---

## 20. 测试 2：当前 C++ kernel direct vs butterfly

完成 reference test 后，再测试当前外推核：

\[
K=W\cdot AA\cdot H_{\rm fast}.
\]

测试：

1. receiver single frequency；
2. source `cmp=0` single frequency；
3. source `cmp=1` single gather；
4. full pipeline。

建议阈值：

```text
relative L2 error < 1e-3
```

如果当前核测试达不到 `1e-3`，先检查：

1. phase 是否来自 \(H_{\rm fast}\)；
2. switch 是否用了完整核；
3. 是否漏掉 antialias；
4. 是否漏掉 taper；
5. aperture 边界是否导致误差；
6. Chebyshev order 是否太低。

---

## 21. 推荐开发顺序

```text
1. 复现 MATLAB 示例核 butterfly，确保 reference test 通过
2. 实现 current KernelEvaluator full K
3. 实现 H-only phase factor
4. 实现 receiver direct vs butterfly
5. 实现 source cmp=0 direct vs butterfly
6. 实现 cmp=1 gather map
7. 实现 source cmp=1 local direct vs butterfly
8. 接入完整频率域外推代码
9. 比较 use_butterfly=0 和 use_butterfly=1
```

---

## 22. 最终验收标准

### MATLAB reference kernel

```text
relative L2 error < 1e-4
```

### Current receiver kernel

```text
relative L2 error < 1e-3
```

### Current source cmp=0 kernel

```text
relative L2 error < 1e-3
```

### Current source cmp=1 local gather

```text
relative L2 error < 1e-3
```

### Full wavefield extrapolation

```text
use_butterfly=0 vs use_butterfly=1
relative L2 error < 1e-3
correlation > 0.999
```

---

# 23. 最终提醒

必须完整保留 MATLAB 代码中的 butterfly 要素：

- binary tree；
- Chebyshev nodes；
- barycentric Lagrange interpolation；
- phase-compensated initialization；
- first-half recursion；
- switch；
- second-half recursion；
- phase-compensated termination；
- direct vs butterfly validation。

但必须替换 MATLAB 核函数：

\[
c_\theta e^{+i\omega t}
\]

为当前 C++ 外推核：

\[
W\cdot AA\cdot H_{\rm fast}.
\]

相位残差使用：

\[
P(a,b)=
\frac{
H_{\rm fast}(a)\overline{H_{\rm fast}(b)}
}{
|H_{\rm fast}(a)\overline{H_{\rm fast}(b)}|+\epsilon
}.
\]

这就是将 MATLAB 相位残差蝶形算法正确迁移到当前频率域波场外推代码中的关键。
