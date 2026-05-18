# SEBF 蝶形算法 C++ 接口实现指南

本文档用于指导另一个 AI 或程序员把已上传的 MATLAB 蝶形算法 `bf_kirch` 改写为 C++ 版本，并整理到项目的 `SEBF` 文件夹中。目标不是重写当前的 FFT、文件读写或 Kirchhoff 外推主程序，而是提供一个可以替换频率域空间求和的蝶形算法接口。

最终目标接口是前面讨论过的固定频率一维振荡积分算子：

```cpp
void butterfly_apply_1d(
    const std::vector<std::complex<double>>& x,
    std::vector<std::complex<double>>& y,
    const std::vector<double>& source_grid,
    const std::vector<double>& target_grid,
    const BFOptions1D& opt,
    const BFKernel1D& kernel);
```

它计算：

\[
y(r_i,\omega)=\sum_j K(r_i,s_j,\omega)x(s_j,\omega),
\]

其中核函数写成：

\[
K(r,s,\omega)=a(r,s,\omega)\exp\left(i\,\epsilon\,\omega\,\phi(r,s)\right),
\]

\(\epsilon\) 由 `opt.phase_sign` 给出。对 MATLAB 测试程序，`phase_sign=+1`。对 FFTW 正变换约定下的当前频率域 Kirchhoff 外推，通常使用 `phase_sign=-1`。

---

## 1. MATLAB 版本代码结构总结

上传的 MATLAB 目录包含以下文件：

```text
bf_kirch/
  bary_weights_1d.m
  bf_test.m
  build_binary_tree_1d.m
  build_lagrange_matrix_1d.m
  interp2_cubic.m
  interp2_linear.m
  lagrange_basis_all_1d.m
```

其中核心测试脚本是 `bf_test.m`。它完成了以下步骤。

### 1.1 构造直接 Kirchhoff 矩阵向量乘法

MATLAB 中直接求和为：

\[
u(r_j)=\sum_i c(s_i,r_j)\exp(i\omega \tau(s_i,r_j))f(s_i).
\]

脚本中使用：

```matlab
kirch_direct(j)=kirch_direct(j)+ctheta(i,j)*exp(1j*omega*tt(i,j))*f(i);
```

这里 `tt(i,j)` 是走时，`ctheta(i,j)` 是振幅权重。

### 1.2 建立 source 和 target 的一维二叉树

`build_binary_tree_1d.m` 对 source 和 target 网格分别建立二叉树。树的 level 为：

```text
root level = 0
leaf level = L
```

要求：

\[
\frac{N}{\texttt{leaf\_n}}=2^L.
\]

每个盒子保存：

```text
i1, i2     当前盒子覆盖的网格点范围
x0         几何中心
w          盒子宽度
cheb       映射到该盒子的 Chebyshev 节点
```

### 1.3 源叶子初始化

初始化阶段 target 为 root，source 为 leaf。系数含义为：

\[
\sigma_t^{A,B,-}\approx \sum_{s\in B} L_t^B(s)\exp\{i\omega[\tau(r_0^A,s)-\tau(r_0^A,s_t^B)]\}f(s).
\]

其中 \(A\) 是 target root box，\(B\) 是 source leaf box，\(s_t^B\) 是 source box 的第 \(t\) 个 Chebyshev 节点。

### 1.4 第一段递推

第一段递推同时让：

```text
target level: 0 -> mid
source level: L -> L-mid
```

即 target 向下细分，source 向上合并。

递推公式为：

\[
\sigma_t^{A,B,-}=\exp[-i\omega\tau(r_0^A,s_t^B)]
\sum_{B_c\subset B}\sum_{t'}L_t^B(s_{t'}^{B_c})
\exp[i\omega\tau(r_0^A,s_{t'}^{B_c})]
\sigma_{t'}^{A_p,B_c,-}.
\]

### 1.5 中间 switch

在中间层，将 source-side 系数转换为 target-side 系数：

\[
\sigma_t^{A,B,+}=\sum_{t'}a(r_t^A,s_{t'}^B)
\exp[i\omega\tau(r_t^A,s_{t'}^B)]
\sigma_{t'}^{A,B,-}.
\]

MATLAB 中的 `K_switch` 原来使用：

```matlab
K_switch(it,jt)=amp_val*exp(1j*omega*tau_val);
```

C++ 版本不能继续固定使用 MATLAB 的 `ctheta` 权重，而必须通过统一 kernel callback 计算新的 Kirchhoff 频率域核。

### 1.6 第二段递推

第二段递推继续让：

```text
target level: mid -> L
source level: L-mid -> 0
```

公式为：

\[
\sigma_t^{A,B,+}=\sum_{B_c\subset B}
\exp[i\omega\tau(r_t^A,s_0^{B_c})]
\sum_{t'}L_{t'}^{A_p}(r_t^A)
\exp[-i\omega\tau(r_{t'}^{A_p},s_0^{B_c})]
\sigma_{t'}^{A_p,B_c,+}.
\]

### 1.7 终止重构

target 到 leaf，source 到 root 后，对每个 target leaf box 重构：

\[
u(r)=\exp[i\omega\tau(r,s_0^B)]
\sum_t L_t^A(r)
\exp[-i\omega\tau(r_t^A,s_0^B)]
\sigma_t^{A,B,+}.
\]

最后与直接求和 `kirch_direct` 比较相对 L2 误差。

---

## 2. 与当前频率域 Kirchhoff 外推代码的核函数关系

当前频率域代码的空间求和形式是：

\[
U_{out}(r,\omega)=\sum_s K(r,s,\omega)U_{in}(s,\omega).
\]

快速版中，单端核函数可以写成：

\[
K(r,s,\omega)=W(r,s)H(\omega,\tau(r,s)).
\]

其中：

\[
W(r,s)=w_{tap}(r,s)\frac{\Delta x}{\pi}
\frac{z_d\tau(r,s)}{z_d^2+(s-r)^2}.
\]

快速版的 \(H\) 是：

\[
H(\omega,\tau)=
\left[(1-\delta)e^{-i\theta m}+\delta e^{-i\theta(m-1)}\right]
\sum_{k=0}^{N-2}\frac{F_k(\tau)}{\Delta t}e^{-i\theta k},
\]

其中：

\[
\theta=\omega\Delta t=\frac{2\pi i_w}{N_{fft}},
\]

\[
q=\frac{\tau}{\Delta t},\qquad m=\lceil q\rceil,\qquad \delta=m-q.
\]

若要放入蝶形算法的标准形式：

\[
K(r,s,\omega)=a(r,s,\omega)\exp(i\epsilon\omega\phi(r,s)),
\]

推荐设置：

```text
phase_sign = -1
phi(r,s)   = tau(r,s)
omega      = 2*pi*iw/(nfft*dt)
```

并令：

\[
a(r,s,\omega)=W(r,s)H(\omega,\tau(r,s))
\exp[-i\epsilon\omega\tau(r,s)].
\]

这样：

\[
a(r,s,\omega)\exp(i\epsilon\omega\tau)=W(r,s)H(\omega,\tau),
\]

与当前频率域快速版直接求和核一致。

---

## 3. SEBF 文件组织

建议在项目中新增：

```text
SEBF/
  include/
    sebf_butterfly1d.h
  src/
    sebf_butterfly1d.cpp
  tests/
    test_sebf_butterfly1d.cpp
  CMakeLists.txt
```

如果当前项目不用 CMake，也可以直接把 `src/sebf_butterfly1d.cpp` 加入 Makefile，并把 `SEBF/include` 加到 include path。

---

## 4. 文件：`SEBF/include/sebf_butterfly1d.h`

```cpp
#ifndef SEBF_BUTTERFLY1D_H
#define SEBF_BUTTERFLY1D_H

#include <complex>
#include <vector>

namespace sebf {

struct BFBox1D {
    int i1 = 0;
    int i2 = 0;
    double x0 = 0.0;
    double width = 0.0;
    std::vector<double> cheb;
};

using BFTree1D = std::vector<std::vector<BFBox1D>>;

struct BFOptions1D {
    int p = 32;
    int leaf_n = 64;
    double omega = 0.0;
    int phase_sign = +1;
    bool verbose = false;
};

using BFPhaseFunc = double (*)(double r, double s,
                               int iout, int iin,
                               const void *user);

using BFAmpFunc = std::complex<double> (*)(double r, double s,
                                           int iout, int iin,
                                           int iw,
                                           const void *user);

struct BFKernel1D {
    BFPhaseFunc phase = nullptr;
    BFAmpFunc amplitude = nullptr;
    const void *user = nullptr;
    int iw = 0;
};

std::vector<double> chebyshev_nodes_standard(int p);
std::vector<double> bary_weights_1d(const std::vector<double> &nodes);
std::vector<double> lagrange_basis_all_1d(double x,
                                          const std::vector<double> &nodes,
                                          const std::vector<double> &weights);
std::vector<double> build_lagrange_matrix_1d(const std::vector<double> &points,
                                             const std::vector<double> &nodes,
                                             const std::vector<double> &weights);
BFTree1D build_binary_tree_1d(const std::vector<double> &grid,
                              int leaf_n,
                              const std::vector<double> &cheb_std);

void butterfly_apply_1d(const std::vector<std::complex<double>> &x,
                        std::vector<std::complex<double>> &y,
                        const std::vector<double> &source_grid,
                        const std::vector<double> &target_grid,
                        const BFOptions1D &opt,
                        const BFKernel1D &kernel);

void direct_apply_1d(const std::vector<std::complex<double>> &x,
                     std::vector<std::complex<double>> &y,
                     const std::vector<double> &source_grid,
                     const std::vector<double> &target_grid,
                     const BFOptions1D &opt,
                     const BFKernel1D &kernel);

double relative_l2_error(const std::vector<std::complex<double>> &a,
                         const std::vector<std::complex<double>> &b);

} // namespace sebf

#endif

```

---

## 5. 文件：`SEBF/src/sebf_butterfly1d.cpp`

```cpp
#include "sebf_butterfly1d.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace sebf {

namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kNodeTol = 1e-13;

inline size_t idx3(int ibox_s, int ibox_r, int it, int nbox_r, int p)
{
    return (static_cast<size_t>(ibox_s) * nbox_r + ibox_r) * p + it;
}

bool is_power_of_two_int(int n)
{
    return n > 0 && (n & (n - 1)) == 0;
}

int integer_log2(int n)
{
    int L = 0;
    while (n > 1) {
        n >>= 1;
        ++L;
    }
    return L;
}

std::complex<double> exp_phase(int sign, double omega, double phase)
{
    const double a = static_cast<double>(sign) * omega * phase;
    return std::complex<double>(std::cos(a), std::sin(a));
}

} // namespace

std::vector<double> chebyshev_nodes_standard(int p)
{
    if (p <= 0) throw std::runtime_error("p must be positive");
    std::vector<double> x(p);
    for (int k = 0; k < p; ++k) {
        x[k] = std::cos((2.0 * (k + 1) - 1.0) * kPi / (2.0 * p));
    }
    return x;
}

std::vector<double> bary_weights_1d(const std::vector<double> &nodes)
{
    const int q = static_cast<int>(nodes.size());
    std::vector<double> w(q, 1.0);
    for (int j = 0; j < q; ++j) {
        for (int k = 0; k < q; ++k) {
            if (k != j) w[j] /= (nodes[j] - nodes[k]);
        }
    }
    return w;
}

std::vector<double> lagrange_basis_all_1d(double x,
                                          const std::vector<double> &nodes,
                                          const std::vector<double> &weights)
{
    const int q = static_cast<int>(nodes.size());
    if (static_cast<int>(weights.size()) != q) {
        throw std::runtime_error("nodes and weights size mismatch");
    }

    std::vector<double> L(q, 0.0);
    for (int j = 0; j < q; ++j) {
        if (std::abs(x - nodes[j]) < kNodeTol) {
            L[j] = 1.0;
            return L;
        }
    }

    double denom = 0.0;
    for (int j = 0; j < q; ++j) {
        L[j] = weights[j] / (x - nodes[j]);
        denom += L[j];
    }
    for (int j = 0; j < q; ++j) L[j] /= denom;
    return L;
}

std::vector<double> build_lagrange_matrix_1d(const std::vector<double> &points,
                                             const std::vector<double> &nodes,
                                             const std::vector<double> &weights)
{
    const int m = static_cast<int>(points.size());
    const int q = static_cast<int>(nodes.size());
    std::vector<double> M(static_cast<size_t>(m) * q, 0.0);
    for (int i = 0; i < m; ++i) {
        std::vector<double> L = lagrange_basis_all_1d(points[i], nodes, weights);
        for (int j = 0; j < q; ++j) M[static_cast<size_t>(i) * q + j] = L[j];
    }
    return M;
}

BFTree1D build_binary_tree_1d(const std::vector<double> &grid,
                              int leaf_n,
                              const std::vector<double> &cheb_std)
{
    const int n = static_cast<int>(grid.size());
    if (leaf_n <= 0 || n % leaf_n != 0 || !is_power_of_two_int(n / leaf_n)) {
        throw std::runtime_error("n/leaf_n must be a positive power of two");
    }

    const int L = integer_log2(n / leaf_n);
    BFTree1D tree(L + 1);

    for (int lev = 0; lev <= L; ++lev) {
        const int nbox = 1 << lev;
        const int box_n = n / nbox;
        tree[lev].resize(nbox);

        for (int b = 0; b < nbox; ++b) {
            const int i1 = b * box_n;
            const int i2 = (b + 1) * box_n - 1;
            const double xL = grid[i1];
            const double xR = grid[i2];
            const double x0 = 0.5 * (xL + xR);
            const double width = xR - xL;

            BFBox1D box;
            box.i1 = i1;
            box.i2 = i2;
            box.x0 = x0;
            box.width = width;
            box.cheb.resize(cheb_std.size());
            for (size_t k = 0; k < cheb_std.size(); ++k) {
                box.cheb[k] = x0 + 0.5 * width * cheb_std[k];
            }
            tree[lev][b] = std::move(box);
        }
    }
    return tree;
}

void direct_apply_1d(const std::vector<std::complex<double>> &x,
                     std::vector<std::complex<double>> &y,
                     const std::vector<double> &source_grid,
                     const std::vector<double> &target_grid,
                     const BFOptions1D &opt,
                     const BFKernel1D &kernel)
{
    if (!kernel.phase || !kernel.amplitude) {
        throw std::runtime_error("kernel callbacks must be provided");
    }
    if (x.size() != source_grid.size()) {
        throw std::runtime_error("x and source_grid size mismatch");
    }

    const int ns = static_cast<int>(source_grid.size());
    const int nr = static_cast<int>(target_grid.size());
    y.assign(nr, std::complex<double>(0.0, 0.0));

    for (int ir = 0; ir < nr; ++ir) {
        const double r = target_grid[ir];
        for (int is = 0; is < ns; ++is) {
            const double s = source_grid[is];
            const double ph = kernel.phase(r, s, ir, is, kernel.user);
            const std::complex<double> amp = kernel.amplitude(r, s, ir, is, kernel.iw, kernel.user);
            y[ir] += amp * exp_phase(opt.phase_sign, opt.omega, ph) * x[is];
        }
    }
}

void butterfly_apply_1d(const std::vector<std::complex<double>> &x,
                        std::vector<std::complex<double>> &y,
                        const std::vector<double> &source_grid,
                        const std::vector<double> &target_grid,
                        const BFOptions1D &opt,
                        const BFKernel1D &kernel)
{
    if (!kernel.phase || !kernel.amplitude) {
        throw std::runtime_error("kernel callbacks must be provided");
    }
    if (x.size() != source_grid.size()) {
        throw std::runtime_error("x and source_grid size mismatch");
    }
    if (source_grid.empty() || target_grid.empty()) {
        throw std::runtime_error("empty source or target grid");
    }
    if (opt.p <= 0 || opt.leaf_n <= 0) {
        throw std::runtime_error("invalid butterfly options");
    }

    const int ns = static_cast<int>(source_grid.size());
    const int nr = static_cast<int>(target_grid.size());
    if (ns != nr) {
        throw std::runtime_error("this compact implementation requires ns == nr; pad or use a rectangular implementation otherwise");
    }
    if (ns % opt.leaf_n != 0 || nr % opt.leaf_n != 0 ||
        !is_power_of_two_int(ns / opt.leaf_n) || !is_power_of_two_int(nr / opt.leaf_n)) {
        throw std::runtime_error("n/leaf_n must be a power of two for both source and target");
    }
    const int Ls = integer_log2(ns / opt.leaf_n);
    const int Lr = integer_log2(nr / opt.leaf_n);
    if (Ls != Lr) {
        throw std::runtime_error("source and target trees must have the same depth in this implementation");
    }
    const int L = Ls;
    const int p = opt.p;
    const int mid_lv = L / 2;

    const std::vector<double> cheb_std = chebyshev_nodes_standard(p);
    const BFTree1D TS = build_binary_tree_1d(source_grid, opt.leaf_n, cheb_std);
    const BFTree1D TR = build_binary_tree_1d(target_grid, opt.leaf_n, cheb_std);

    // ------------------------------------------------------------
    // Initialization: target root, source leaves.
    // sigma(B_leaf, A_root, t_s) are source-side coefficients.
    // ------------------------------------------------------------
    {
        const std::vector<BFBox1D> &boxes_s = TS[L];
        const std::vector<BFBox1D> &boxes_r = TR[0];
        const int nbox_s = static_cast<int>(boxes_s.size());
        const int nbox_r = static_cast<int>(boxes_r.size());
        std::vector<std::complex<double>> sigma(static_cast<size_t>(nbox_s) * nbox_r * p,
                                                std::complex<double>(0.0, 0.0));

        for (int ibs = 0; ibs < nbox_s; ++ibs) {
            const BFBox1D &B = boxes_s[ibs];
            const BFBox1D &A = boxes_r[0];
            const double r0 = A.x0;
            const std::vector<double> weightsB = bary_weights_1d(B.cheb);

            std::vector<double> src_points;
            src_points.reserve(B.i2 - B.i1 + 1);
            for (int i = B.i1; i <= B.i2; ++i) src_points.push_back(source_grid[i]);
            const std::vector<double> Lg = build_lagrange_matrix_1d(src_points, B.cheb, weightsB);

            for (int t = 0; t < p; ++t) {
                std::complex<double> acc(0.0, 0.0);
                const double phase_node = kernel.phase(r0, B.cheb[t], -1, -1, kernel.user);
                const std::complex<double> fac_node = exp_phase(-opt.phase_sign, opt.omega, phase_node);

                for (int ii = B.i1; ii <= B.i2; ++ii) {
                    const int local = ii - B.i1;
                    const double s = source_grid[ii];
                    const double phase_s = kernel.phase(r0, s, -1, ii, kernel.user);
                    const double Lval = Lg[static_cast<size_t>(local) * p + t];
                    acc += Lval * exp_phase(opt.phase_sign, opt.omega, phase_s) * x[ii];
                }
                sigma[idx3(ibs, 0, t, nbox_r, p)] = acc * fac_node;
            }
        }

        std::vector<std::complex<double>> sigma_prev = std::move(sigma);

        // ------------------------------------------------------------
        // First-half recursion: source boxes go upward, target boxes go downward.
        // ------------------------------------------------------------
        for (int lv = 1; lv <= mid_lv; ++lv) {
            const std::vector<BFBox1D> &boxes_r_now = TR[lv];
            const std::vector<BFBox1D> &boxes_s_now = TS[L - lv];
            const std::vector<BFBox1D> &boxes_s_prev = TS[L - lv + 1];

            const int nbox_r_now = static_cast<int>(boxes_r_now.size());
            const int nbox_s_now = static_cast<int>(boxes_s_now.size());
            const int nbox_r_prev = static_cast<int>(TR[lv - 1].size());

            std::vector<std::complex<double>> sigma_curr(static_cast<size_t>(nbox_s_now) * nbox_r_now * p,
                                                         std::complex<double>(0.0, 0.0));

            for (int ibr = 0; ibr < nbox_r_now; ++ibr) {
                const BFBox1D &A = boxes_r_now[ibr];
                const double r0 = A.x0;
                const int ibr_parent = ibr / 2;

                for (int ibs = 0; ibs < nbox_s_now; ++ibs) {
                    const BFBox1D &B = boxes_s_now[ibs];
                    const std::vector<double> weightsB = bary_weights_1d(B.cheb);
                    std::vector<std::complex<double>> acc(p, std::complex<double>(0.0, 0.0));

                    for (int child = 0; child < 2; ++child) {
                        const int ibs_child = 2 * ibs + child;
                        const BFBox1D &Bc = boxes_s_prev[ibs_child];
                        const std::vector<double> Lmat = build_lagrange_matrix_1d(Bc.cheb, B.cheb, weightsB);

                        std::vector<std::complex<double>> tmp(p);
                        for (int jt = 0; jt < p; ++jt) {
                            const double ph = kernel.phase(r0, Bc.cheb[jt], -1, -1, kernel.user);
                            tmp[jt] = exp_phase(opt.phase_sign, opt.omega, ph) *
                                      sigma_prev[idx3(ibs_child, ibr_parent, jt, nbox_r_prev, p)];
                        }

                        for (int t = 0; t < p; ++t) {
                            for (int jt = 0; jt < p; ++jt) {
                                acc[t] += Lmat[static_cast<size_t>(jt) * p + t] * tmp[jt];
                            }
                        }
                    }

                    for (int t = 0; t < p; ++t) {
                        const double ph = kernel.phase(r0, B.cheb[t], -1, -1, kernel.user);
                        sigma_curr[idx3(ibs, ibr, t, nbox_r_now, p)] =
                            acc[t] * exp_phase(-opt.phase_sign, opt.omega, ph);
                    }
                }
            }
            sigma_prev = std::move(sigma_curr);
        }

        // ------------------------------------------------------------
        // Switch: source-side nodal coefficients -> target-side nodal coefficients.
        // The full physical amplitude is applied here.
        // ------------------------------------------------------------
        {
            const std::vector<BFBox1D> &boxes_r_mid = TR[mid_lv];
            const std::vector<BFBox1D> &boxes_s_mid = TS[L - mid_lv];
            const int nbox_r_mid = static_cast<int>(boxes_r_mid.size());
            const int nbox_s_mid = static_cast<int>(boxes_s_mid.size());

            std::vector<std::complex<double>> sigma_sw(static_cast<size_t>(nbox_s_mid) * nbox_r_mid * p,
                                                       std::complex<double>(0.0, 0.0));

            for (int ibr = 0; ibr < nbox_r_mid; ++ibr) {
                const BFBox1D &A = boxes_r_mid[ibr];
                for (int ibs = 0; ibs < nbox_s_mid; ++ibs) {
                    const BFBox1D &B = boxes_s_mid[ibs];
                    for (int it = 0; it < p; ++it) {
                        const double r = A.cheb[it];
                        std::complex<double> acc(0.0, 0.0);
                        for (int jt = 0; jt < p; ++jt) {
                            const double s = B.cheb[jt];
                            const double ph = kernel.phase(r, s, -1, -1, kernel.user);
                            const std::complex<double> amp = kernel.amplitude(r, s, -1, -1, kernel.iw, kernel.user);
                            acc += amp * exp_phase(opt.phase_sign, opt.omega, ph) *
                                   sigma_prev[idx3(ibs, ibr, jt, nbox_r_mid, p)];
                        }
                        sigma_sw[idx3(ibs, ibr, it, nbox_r_mid, p)] = acc;
                    }
                }
            }
            sigma_prev = std::move(sigma_sw);
        }

        // ------------------------------------------------------------
        // Second-half recursion: target boxes go downward to leaves.
        // ------------------------------------------------------------
        for (int lv = mid_lv + 1; lv <= L; ++lv) {
            const std::vector<BFBox1D> &boxes_r_now = TR[lv];
            const std::vector<BFBox1D> &boxes_s_now = TS[L - lv];
            const std::vector<BFBox1D> &boxes_r_prev = TR[lv - 1];
            const std::vector<BFBox1D> &boxes_s_prev = TS[L - lv + 1];

            const int nbox_r_now = static_cast<int>(boxes_r_now.size());
            const int nbox_s_now = static_cast<int>(boxes_s_now.size());
            const int nbox_r_prev = static_cast<int>(boxes_r_prev.size());

            std::vector<std::complex<double>> sigma_curr(static_cast<size_t>(nbox_s_now) * nbox_r_now * p,
                                                         std::complex<double>(0.0, 0.0));

            for (int ibr = 0; ibr < nbox_r_now; ++ibr) {
                const BFBox1D &A = boxes_r_now[ibr];
                const int ibr_parent = ibr / 2;
                const BFBox1D &Ap = boxes_r_prev[ibr_parent];

                const std::vector<double> weightsAp = bary_weights_1d(Ap.cheb);
                const std::vector<double> Lmat = build_lagrange_matrix_1d(A.cheb, Ap.cheb, weightsAp);

                for (int ibs = 0; ibs < nbox_s_now; ++ibs) {
                    std::vector<std::complex<double>> acc(p, std::complex<double>(0.0, 0.0));
                    for (int child = 0; child < 2; ++child) {
                        const int ibs_child = 2 * ibs + child;
                        const BFBox1D &Bc = boxes_s_prev[ibs_child];
                        const double s0 = Bc.x0;

                        std::vector<std::complex<double>> tmp(p);
                        for (int jt = 0; jt < p; ++jt) {
                            const double ph = kernel.phase(Ap.cheb[jt], s0, -1, -1, kernel.user);
                            tmp[jt] = exp_phase(-opt.phase_sign, opt.omega, ph) *
                                      sigma_prev[idx3(ibs_child, ibr_parent, jt, nbox_r_prev, p)];
                        }

                        for (int it = 0; it < p; ++it) {
                            std::complex<double> sum(0.0, 0.0);
                            for (int jt = 0; jt < p; ++jt) {
                                sum += Lmat[static_cast<size_t>(it) * p + jt] * tmp[jt];
                            }
                            const double ph = kernel.phase(A.cheb[it], s0, -1, -1, kernel.user);
                            acc[it] += exp_phase(opt.phase_sign, opt.omega, ph) * sum;
                        }
                    }
                    for (int it = 0; it < p; ++it) {
                        sigma_curr[idx3(ibs, ibr, it, nbox_r_now, p)] = acc[it];
                    }
                }
            }
            sigma_prev = std::move(sigma_curr);
        }

        // ------------------------------------------------------------
        // Termination: reconstruct target leaf boxes.
        // Source is root, target is leaf.
        // ------------------------------------------------------------
        y.assign(nr, std::complex<double>(0.0, 0.0));
        const std::vector<BFBox1D> &boxes_r_leaf = TR[L];
        const BFBox1D &Broot = TS[0][0];
        const double s0 = Broot.x0;
        const int nbox_r_leaf = static_cast<int>(boxes_r_leaf.size());

        for (int ibr = 0; ibr < nbox_r_leaf; ++ibr) {
            const BFBox1D &A = boxes_r_leaf[ibr];
            std::vector<double> target_points;
            target_points.reserve(A.i2 - A.i1 + 1);
            for (int i = A.i1; i <= A.i2; ++i) target_points.push_back(target_grid[i]);

            const std::vector<double> weightsA = bary_weights_1d(A.cheb);
            const std::vector<double> LgA = build_lagrange_matrix_1d(target_points, A.cheb, weightsA);

            std::vector<std::complex<double>> tmp(p);
            for (int it = 0; it < p; ++it) {
                const double ph = kernel.phase(A.cheb[it], s0, -1, -1, kernel.user);
                tmp[it] = exp_phase(-opt.phase_sign, opt.omega, ph) *
                          sigma_prev[idx3(0, ibr, it, nbox_r_leaf, p)];
            }

            for (int local = 0; local <= A.i2 - A.i1; ++local) {
                const int ir = A.i1 + local;
                std::complex<double> sum(0.0, 0.0);
                for (int it = 0; it < p; ++it) {
                    sum += LgA[static_cast<size_t>(local) * p + it] * tmp[it];
                }
                const double ph = kernel.phase(target_grid[ir], s0, ir, -1, kernel.user);
                y[ir] = exp_phase(opt.phase_sign, opt.omega, ph) * sum;
            }
        }
    }
}

double relative_l2_error(const std::vector<std::complex<double>> &a,
                         const std::vector<std::complex<double>> &b)
{
    if (a.size() != b.size()) throw std::runtime_error("relative_l2_error size mismatch");
    double num = 0.0;
    double den = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const std::complex<double> d = a[i] - b[i];
        num += std::norm(d);
        den += std::norm(b[i]);
    }
    return std::sqrt(num / std::max(den, 1e-300));
}

} // namespace sebf

```

---

## 6. 文件：`SEBF/tests/test_sebf_butterfly1d.cpp`

这个测试不依赖 SE 库，只验证蝶形算法接口本身。它构造一个类似 MATLAB `bf_test.m` 的单频 Kirchhoff 核：

\[
K(r,s)=a(r,s)\exp(i\omega\tau(r,s)),
\]

并比较：

```text
direct_apply_1d vs butterfly_apply_1d
```

相对 L2 误差应为较低水平。当前设置下通常可达到约 \(10^{-4}\sim10^{-3}\)。如果误差偏大，可增大 `p` 或减小 `leaf_n`。

```cpp
#include "sebf_butterfly1d.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

struct TestKernelData {
    double v = 3000.0;
    double dz = 100.0;
};

static double test_phase(double r, double s, int, int, const void *user)
{
    const TestKernelData *d = static_cast<const TestKernelData *>(user);
    const double dx = s - r;
    return std::sqrt(dx * dx + d->dz * d->dz) / d->v;
}

static std::complex<double> test_amp(double r, double s, int, int, int, const void *user)
{
    const TestKernelData *d = static_cast<const TestKernelData *>(user);
    const double dx = s - r;
    const double R = std::sqrt(dx * dx + d->dz * d->dz);
    const double amp = d->dz / std::pow(R, 1.5);
    return std::complex<double>(amp, 0.0);
}

int main()
{
    const int n = 1024;
    const double dx = 5.0;
    std::vector<double> sgrid(n), rgrid(n);
    for (int i = 0; i < n; ++i) {
        sgrid[i] = (i + 1) * dx;
        rgrid[i] = (i + 1) * dx;
    }

    std::vector<std::complex<double>> x(n, std::complex<double>(0.0, 0.0));
    x[511] = std::complex<double>(1.0, 0.0);
    x[512] = std::complex<double>(1.0, 0.0);

    TestKernelData data;
    sebf::BFOptions1D opt;
    opt.p = 60;
    opt.leaf_n = 64;
    opt.omega = 15.0 * 2.0 * 3.14159265358979323846;
    opt.phase_sign = +1;

    sebf::BFKernel1D kernel;
    kernel.phase = test_phase;
    kernel.amplitude = test_amp;
    kernel.user = &data;
    kernel.iw = 0;

    std::vector<std::complex<double>> y_direct, y_bf;
    sebf::direct_apply_1d(x, y_direct, sgrid, rgrid, opt, kernel);
    sebf::butterfly_apply_1d(x, y_bf, sgrid, rgrid, opt, kernel);

    const double rel = sebf::relative_l2_error(y_bf, y_direct);
    std::cout << "relative L2 error = " << rel << "\n";
    return (rel < 5e-3) ? 0 : 1;
}

```

---

## 7. 文件：`SEBF/CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.10)
project(SEBF LANGUAGES CXX)

add_library(sebf STATIC
    src/sebf_butterfly1d.cpp
)

target_include_directories(sebf PUBLIC
    $${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(sebf PUBLIC cxx_std_11)

add_executable(test_sebf_butterfly1d
    tests/test_sebf_butterfly1d.cpp
)

target_link_libraries(test_sebf_butterfly1d PRIVATE sebf)

```

如果不用 CMake，可使用以下命令直接测试。

```bash
# Minimal standalone build command for the SEBF test
g++ -O2 -std=c++11 \
    -I SEBF/include \
    SEBF/src/sebf_butterfly1d.cpp \
    SEBF/tests/test_sebf_butterfly1d.cpp \
    -o test_sebf_butterfly1d

./test_sebf_butterfly1d

```

---

## 8. 如何把 SEBF 接入当前频率域 Kirchhoff 代码

当前快速频率域代码中，receiver-side 和 source-side 的直接求和分别是：

\[
U_{tmp}(i_s,i_h,\omega)=\sum_{i_c}K_r(i_s,i_h,i_c,\omega)U_{in}(i_s,i_c,\omega),
\]

\[
U_{out}(i_s,i_h,\omega)=\sum_{i_c}K_s(i_h,i_s,i_c,\omega)U_{tmp}(i_c,i_h,\omega).
\]

蝶形算法只替换这两个求和，不修改：

```text
FFT/IFFT
文件输入输出
Green 表读取
时间反转
nfft 设置
H(omega,tau) 的定义
```

### 8.1 receiver-side 调用方式

对每个固定的 `iw` 和 `is`，调用一次：

```cpp
std::vector<std::complex<double>> x(ctx->nh), y;
for (int ic = 0; ic < ctx->nh; ++ic) {
    size_t idx = IDX3(is, ic, iw, ctx->nh, ctx->nw);
    x[ic] = std::complex<double>(Uin[idx][0], Uin[idx][1]);
}

sebf::BFOptions1D opt;
opt.p = p;
opt.leaf_n = leaf_n;
opt.omega = 2.0 * M_PI * iw / (ctx->nfft * ctx->dt);
opt.phase_sign = -1;

sebf::BFKernel1D kernel;
kernel.phase = receiver_phase_callback;
kernel.amplitude = receiver_amplitude_callback;
kernel.user = &receiver_user_data;
kernel.iw = iw;

sebf::butterfly_apply_1d(x, y, source_grid_h, target_grid_h, opt, kernel);

for (int ih = 0; ih < ctx->nh; ++ih) {
    size_t idx = IDX3(is, ih, iw, ctx->nh, ctx->nw);
    Utmp[idx][0] += y[ih].real();
    Utmp[idx][1] += y[ih].imag();
}
```

这里 `source_grid_h` 和 `target_grid_h` 都是检波点坐标：

```cpp
source_grid_h[ic] = ctx->h0 + ic * ctx->dh;
target_grid_h[ih] = ctx->h0 + ih * ctx->dh;
```

### 8.2 source-side `cmp=0` 调用方式

对每个固定的 `iw` 和 `ih`，调用一次：

```cpp
std::vector<std::complex<double>> x(ctx->ns), y;
for (int ic = 0; ic < ctx->ns; ++ic) {
    size_t idx = IDX3(ic, ih, iw, ctx->nh, ctx->nw);
    x[ic] = std::complex<double>(Utmp[idx][0], Utmp[idx][1]);
}

sebf::BFOptions1D opt;
opt.p = p;
opt.leaf_n = leaf_n;
opt.omega = 2.0 * M_PI * iw / (ctx->nfft * ctx->dt);
opt.phase_sign = -1;

sebf::BFKernel1D kernel;
kernel.phase = source_phase_callback_cmp0;
kernel.amplitude = source_amplitude_callback_cmp0;
kernel.user = &source_user_data;
kernel.iw = iw;

sebf::butterfly_apply_1d(x, y, source_grid_s, target_grid_s, opt, kernel);

for (int is = 0; is < ctx->ns; ++is) {
    size_t idx = IDX3(is, ih, iw, ctx->nh, ctx->nw);
    Uout[idx][0] += y[is].real();
    Uout[idx][1] += y[is].imag();
}
```

这里：

```cpp
source_grid_s[ic] = ctx->s0 + ic * ctx->ds;
target_grid_s[is] = ctx->s0 + is * ctx->ds;
```

建议先只实现和验证 `cmp=0`，再处理 `cmp=1`。

---

## 9. Kirchhoff 快速频率域核函数 callback 写法

### 9.1 双线性插值走时表

MATLAB 版本用 `interp2_linear(tt,dx,dx,s,r,dx,dx)` 在二维走时表上插值。C++ 中也必须提供类似函数，用于在 Chebyshev 节点处计算 \(\tau(r,s)\)。示例：

```cpp
static double interp_table2_linear(float **tab,
                                   int n1, int n2,
                                   double o1, double d1,
                                   double o2, double d2,
                                   double x1, double x2)
{
    double u = (x1 - o1) / d1;
    double v = (x2 - o2) / d2;

    if (u < 0.0) u = 0.0;
    if (v < 0.0) v = 0.0;
    if (u > n1 - 1.0) u = n1 - 1.0;
    if (v > n2 - 1.0) v = n2 - 1.0;

    int i = (int)floor(u);
    int j = (int)floor(v);

    if (i >= n1 - 1) i = n1 - 2;
    if (j >= n2 - 1) j = n2 - 2;
    if (i < 0) i = 0;
    if (j < 0) j = 0;

    double a = u - i;
    double b = v - j;

    double q11 = tab[i][j];
    double q21 = tab[i + 1][j];
    double q12 = tab[i][j + 1];
    double q22 = tab[i + 1][j + 1];

    return (1.0 - a) * (1.0 - b) * q11
         + a         * (1.0 - b) * q21
         + (1.0 - a) * b         * q12
         + a         * b         * q22;
}
```

注意：上面的 `tab[i][j]` 对应当前 C 代码中的 `rtable[cc][c]` 或 `stable[cc][c]`。第一个坐标是输入点坐标，第二个坐标是输出点坐标。

### 9.2 快速版 H 函数

当前快速版的 `build_H_iw()` 可写为：

```cpp
static std::complex<double> build_H_fast(double tau,
                                         double dt,
                                         int nsam,
                                         int iw,
                                         int nfft)
{
    if (tau <= 0.0) return {0.0, 0.0};

    std::vector<double> G(nsam, 0.0);
    for (int k = 1; k < nsam; ++k) {
        double a = (tau + k * dt) / tau;
        G[k] = std::sqrt(a * a - 1.0);
    }
    for (int k = 0; k < nsam - 1; ++k) G[k] = G[k + 1] - G[k];
    for (int k = nsam - 2; k > 0; --k) G[k] = G[k] - G[k - 1];

    double q = tau / dt;
    int m = (int)std::ceil(q);
    double delta = (double)m - q;
    double theta = 2.0 * M_PI * iw / nfft;

    std::complex<double> e_step(std::cos(-theta), std::sin(-theta));
    std::complex<double> e_k(1.0, 0.0);
    std::complex<double> SF(0.0, 0.0);

    for (int k = 0; k <= nsam - 2; ++k) {
        SF += (G[k] / dt) * e_k;
        e_k *= e_step;
    }

    std::complex<double> e_m(std::cos(-theta * m), std::sin(-theta * m));
    std::complex<double> e_m1(std::cos(-theta * (m - 1)), std::sin(-theta * (m - 1)));
    std::complex<double> C = (1.0 - delta) * e_m + delta * e_m1;

    return C * SF;
}
```

### 9.3 receiver-side callback

receiver-side 物理核为：

\[
K_r(r,s,\omega)=W_r(r,s)H(\omega,\tau_r(s,r)).
\]

在蝶形算法接口中采用：

\[
K_r=a_r(r,s,\omega)\exp(-i\omega\tau_r(s,r)).
\]

因此：

\[
a_r=W_r H\exp(i\omega\tau_r).
\]

示例 callback：

```cpp
struct ReceiverBFUser {
    const FreqContext *ctx;
    int is;
};

static double receiver_phase_callback(double r, double s,
                                      int iout, int iin,
                                      const void *user)
{
    const ReceiverBFUser *u = static_cast<const ReceiverBFUser *>(user);
    const FreqContext *ctx = u->ctx;

    double abs_r = (ctx->cmp ? (ctx->s0 + u->is * ctx->ds) : 0.0) + r;
    double abs_s = (ctx->cmp ? (ctx->s0 + u->is * ctx->ds) : 0.0) + s;

    return interp_table2_linear(ctx->rtable,
                                ctx->nrg, ctx->nrg,
                                ctx->rg0, ctx->drg,
                                ctx->rg0, ctx->drg,
                                abs_s, abs_r);
}

static std::complex<double> receiver_amplitude_callback(double r, double s,
                                                        int iout, int iin,
                                                        int iw,
                                                        const void *user)
{
    const ReceiverBFUser *u = static_cast<const ReceiverBFUser *>(user);
    const FreqContext *ctx = u->ctx;

    double tau = receiver_phase_callback(r, s, iout, iin, user);
    if (tau <= 0.0) return {0.0, 0.0};

    // Aperture mask. If Chebyshev nodes are passed, iout/iin may be -1.
    // In that case use a coordinate-based aperture check.
    double dx = s - r;
    if (std::abs(dx) > ctx->aper * std::abs(ctx->dh)) return {0.0, 0.0};

    double dist = ctx->rdatum * ctx->rdatum + dx * dx;
    double W = ctx->dh / M_PI * ctx->rdatum * tau / dist;

    std::complex<double> H = build_H_fast(tau, ctx->dt, ctx->nsam, iw, ctx->nfft);
    double omega = 2.0 * M_PI * iw / (ctx->nfft * ctx->dt);

    // phase_sign = -1, so amplitude must multiply exp(+i omega tau)
    std::complex<double> correction(std::cos(omega * tau), std::sin(omega * tau));
    return W * H * correction;
}
```

### 9.4 source-side `cmp=0` callback

```cpp
struct SourceBFUser {
    const FreqContext *ctx;
    int ih;
};

static double source_phase_callback_cmp0(double r, double s,
                                         int iout, int iin,
                                         const void *user)
{
    const SourceBFUser *u = static_cast<const SourceBFUser *>(user);
    const FreqContext *ctx = u->ctx;

    return interp_table2_linear(ctx->stable,
                                ctx->nsg, ctx->nsg,
                                ctx->sg0, ctx->dsg,
                                ctx->sg0, ctx->dsg,
                                s, r);
}

static std::complex<double> source_amplitude_callback_cmp0(double r, double s,
                                                           int iout, int iin,
                                                           int iw,
                                                           const void *user)
{
    const SourceBFUser *u = static_cast<const SourceBFUser *>(user);
    const FreqContext *ctx = u->ctx;

    double tau = source_phase_callback_cmp0(r, s, iout, iin, user);
    if (tau <= 0.0) return {0.0, 0.0};

    double dx = s - r;
    if (std::abs(dx) > ctx->aper * std::abs(ctx->ds)) return {0.0, 0.0};

    double dist = ctx->sdatum * ctx->sdatum + dx * dx;
    double W = ctx->ds / M_PI * ctx->sdatum * tau / dist;

    std::complex<double> H = build_H_fast(tau, ctx->dt, ctx->nsam, iw, ctx->nfft);
    double omega = 2.0 * M_PI * iw / (ctx->nfft * ctx->dt);

    std::complex<double> correction(std::cos(omega * tau), std::sin(omega * tau));
    return W * H * correction;
}
```

---

## 10. 重要实现建议

### 10.1 先关闭 aperture 或使用大 aperture 测试

蝶形算法依赖核函数在块上具有低秩结构。如果 aperture mask 是硬截断，会在边界处引入不连续。建议调试时先设置：

```text
aperture >= max(ns, nh)
taper = 0
```

确认蝶形算法本身正确后，再加入 aperture 和 taper。

### 10.2 `cmp=1` 不要第一步实现

`cmp=1` 的 source-side 有：

```cpp
hh = (int)((r - ic * ctx->ds) / ctx->dh + 0.5f);
```

这是非规则共检波点重排。第一版 SEBF 只建议支持 receiver-side 和 source-side `cmp=0`。`cmp=1` 后续应实现 irregular butterfly task。

### 10.3 直接求和和蝶形算法必须使用同一个 callback

为了公平验证：

```text
direct_apply_1d
butterfly_apply_1d
```

必须使用同一组 `phase` 和 `amplitude` callback。否则误差无法归因。

---

## 11. 测试要求

### 11.1 SEBF 单元测试

编译并运行：

```bash
g++ -O2 -std=c++11 \
    -I SEBF/include \
    SEBF/src/sebf_butterfly1d.cpp \
    SEBF/tests/test_sebf_butterfly1d.cpp \
    -o test_sebf_butterfly1d

./test_sebf_butterfly1d
```

期望输出类似：

```text
relative L2 error = 2.0e-04
```

如果误差大于 `5e-3`，先检查：

1. `phase_sign` 是否正确；
2. `p` 是否太小；
3. `leaf_n` 是否太大；
4. amplitude 是否太不光滑；
5. phase 和 amplitude callback 是否与 direct 版本一致。

### 11.2 Kirchhoff 集成测试

在当前频率域代码中保留两个模式：

```text
use_butterfly=0  direct frequency-domain summation
use_butterfly=1  butterfly summation
```

对同一输入运行两次，比较输出：

\[
\epsilon=\frac{\|u_{bf}-u_{direct}\|_2}{\|u_{direct}\|_2}.
\]

建议先只测试：

```text
cmp=0
aperture >= max(ns,nh)
taper=0
单频或较少频率
```

再逐步恢复：

```text
完整频率范围
正常 aperture
正常 taper
cmp=1
```

---

## 12. 与 MATLAB 版本的主要差异

1. MATLAB 版本把 `ctheta` 作为 switch 阶段的权重；C++ 版本必须用当前频率域 Kirchhoff 核：

   \[
   K=W H.
   \]

2. MATLAB 版本 `interp2_linear(tt,...)` 直接插值走时矩阵；C++ 版本必须对 `stable` 和 `rtable` 做对应的二维线性插值。

3. MATLAB 版本默认正相位 `exp(+i omega tau)`；FFTW 频率域 Kirchhoff 外推通常应设置 `phase_sign=-1`。

4. MATLAB 版本只测试 `ns=nr` 的规则一维情况；当前 C++ 紧凑实现也要求 `n_source == n_target` 且 source/target 树深一致。如果实际 `ns != nh`，第一版可通过补零或单独实现 rectangular butterfly 扩展。

---

## 13. 最终落地步骤

建议按以下顺序实施：

```text
1. 新建 SEBF/include/sebf_butterfly1d.h
2. 新建 SEBF/src/sebf_butterfly1d.cpp
3. 新建 SEBF/tests/test_sebf_butterfly1d.cpp
4. 编译并运行 SEBF standalone test
5. 在当前频率域 Kirchhoff 代码中新增 use_butterfly 参数
6. 先接 receiver-side butterfly，source-side 仍用 direct
7. 验证 receiver intermediate 结果
8. 接 source-side cmp=0 butterfly
9. 验证最终结果
10. 最后再考虑 cmp=1 irregular butterfly
```

最重要的原则是：

\[
\boxed{\text{direct 和 butterfly 版本必须使用同一个 }K(r,s,\omega).}
\]

只允许改变空间求和方式，不允许改变核函数、FFT/IFFT、输入输出和走时表读取逻辑。
