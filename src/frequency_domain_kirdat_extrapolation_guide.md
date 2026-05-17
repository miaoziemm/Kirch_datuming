# 时间域 Kirchhoff 波场外推到频率域实现指南

本文档用于指导另一个 AI 或程序员，将当前时间域 `kirdat` 类波场外推代码改写为频率域实现。目标是：

1. 输入文件、输出文件、数据维度、道集组织方式、头信息保持与原程序完全一致；
2. 外推过程严格按照当前时间域代码对应的频率域公式实现；
3. 可以调用 FFTW / FFTW3 的函数进行时间轴 FFT；
4. 若要求与当前时间域程序逐样点一致，必须保留本文中给出的线性插值边界修正项。

---

## 1. 原时间域程序的核心结构

原程序做的是两次单端平面外推：

\[
D_{\rm out}=\mathcal{D}_s\mathcal{D}_rD_{\rm in},
\]

其中：

- \(\mathcal{D}_r\)：检波点端外推，即 common-shot gather 内沿接收点轴外推；
- \(\mathcal{D}_s\)：震源端外推，即 common-receiver gather 内沿震源轴外推；
- 当前代码在读入后先执行一次 `reverse_trace(nt, nh, ns, tr_in)`，最后写出前再执行一次 `reverse_trace(nt, nh, ns, tr_in)`。频率域版本必须保持这一点，否则输出时间方向会不一致。

单端外推的时间域结构为

\[
u_r(r,t_n)
=
\sum_{s\in\mathcal{A}(r)}
W(r,s)
\,\mathcal{T}_{\tau(s,r)}[u_s(s,t)]_n,
\]

其中：

- \(s\)：原始平面上的输入点；
- \(r\)：外推后平面上的输出点；
- \(\mathcal{A}(r)\)：以 \(r\) 为中心的空间孔径；
- \(\tau(s,r)=t(s,r)\)：从 \(s\) 到 \(r\) 的走时，由 Green 函数走时表给出；
- \(W(r,s)\)：Kirchhoff 空间权重；
- \(\mathcal{T}_{\tau}\)：当前代码中的时间延迟滤波算子。

---

## 2. Kirchhoff 空间权重

对任意单端外推，空间权重写为

\[
W(r,s)
=
w(r,s)
\frac{\Delta x}{\pi}
\frac{z_d\,\tau(s,r)}{z_d^2+|s-r|^2},
\]

其中：

- \(w(r,s)\)：代码中的 taper 系数 `coef`；
- \(\Delta x\)：该次外推方向上的空间采样间隔；
- \(z_d\)：两个平面之间的外推距离；
- \(\tau(s,r)\)：走时；
- \(|s-r|\)：输入点和输出点的水平距离。

### 2.1 检波点端外推

在 common-shot gather 中，固定震源 \(x_s\)，对接收点坐标外推。代码对应

\[
W_r(i_h,i_c)
=
\frac{\texttt{coef}}{\pi}
\frac{\Delta h\,z_r\,\tau_r(i_c,i_h)}
{z_r^2+(i_c-i_h)^2\Delta h^2},
\]

其中：

\[
z_r=\texttt{rdatum},\qquad \Delta h=\texttt{dh},
\]

\[
\tau_r(i_c,i_h)=\texttt{rtable[cc][c]}.
\]

这里 \(i_c\) 是输入接收点索引，\(i_h\) 是外推后接收点索引。

### 2.2 震源端外推

在 common-receiver gather 中，对震源坐标外推。代码对应

\[
W_s(i_s,i_c)
=
\frac{\texttt{coef}}{\pi}
\frac{\Delta s\,z_s\,\tau_s(i_c,i_s)}
{z_s^2+(i_c-i_s)^2\Delta s^2},
\]

其中：

\[
z_s=\texttt{sdatum},\qquad \Delta s=\texttt{ds},
\]

\[
\tau_s(i_c,i_s)=\texttt{stable[cc][c]}.
\]

这里 \(i_c\) 是输入震源索引，\(i_s\) 是外推后震源索引。

---

## 3. 时间域滤波核的离散定义

当前代码中，对每一个走时 \(\tau\)，先构造

\[
G_k(\tau)
=
\sqrt{
\left(
\frac{\tau+k\Delta t}{\tau}
\right)^2-1
},
\qquad k=0,1,\ldots,N-1,
\]

其中

\[
N=\texttt{nsam}=\left\lfloor \frac{\texttt{length}}{\Delta t}\right\rfloor+2.
\]

然后执行一次一阶差分和一次二阶差分。最终滤波系数为

\[
F_0(\tau)=G_1(\tau)-G_0(\tau),
\]

\[
F_k(\tau)=G_{k+1}(\tau)-2G_k(\tau)+G_{k-1}(\tau),
\qquad k=1,2,\ldots,N-2.
\]

因此有效滤波长度为

\[
L_f=N-1.
\]

---

## 4. 与当前 C 代码完全一致的时间域单对点公式

当前代码不是直接使用 `floor(tau/dt)`，而是从第一个满足

\[
t_n=n\Delta t\geq \tau
\]

的输出时间样点开始计算。因此定义

\[
m=\left\lceil \frac{\tau}{\Delta t}\right\rceil,
\]

\[
\delta=m-\frac{\tau}{\Delta t},
\qquad 0\leq\delta<1.
\]

代码中的 `shift` 等于

\[
\texttt{shift}=n-m.
\]

对于一个输入道 \(x_j\) 和一个输出道 \(y_n\)，单对点贡献为

\[
y_n
\;{+}{=}\;
W(r,s)
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}
F_k(\tau)
\left[(1-\delta)x_{n-m-k}+\delta x_{n-m-k+1}\right]
I(n\geq m+k),
\]

其中 \(I(\cdot)\) 是示性函数。这个示性函数来自代码中的条件

```c
(isam < nsam-1) && (shift-isam >= 0)
```

即

\[
k\leq n-m.
\]

这个公式是最重要的：如果需要逐样点复现时间域代码，频率域实现必须等价于这个公式。

---

## 5. 频率域公式：写成 \(\exp[i\omega t(s,r)]\) 的形式

为了写成用户需要的相位形式，采用如下傅里叶约定：

\[
U(\omega)=\sum_{n=0}^{N_t-1}u_n\exp(i\omega n\Delta t),
\]

\[
u_n=\frac{1}{N_\omega}\sum_{\ell=0}^{N_\omega-1}
U(\omega_\ell)\exp(-i\omega_\ell n\Delta t).
\]

在这个约定下，延迟 \(u(t-\tau)\) 对应频率域因子

\[
\exp(i\omega\tau).
\]

### 5.1 不考虑边界截断时的主频率响应

先定义线性卷积核

\[
h_{\tau,p}
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)
\left[(1-\delta)\,\mathbf{1}_{p=m+k}
+
\delta\,\mathbf{1}_{p=m+k-1}
\right].
\]

其频率响应为

\[
H^+(\omega,\tau)
=
\sum_p h_{\tau,p}\exp(i\omega p\Delta t).
\]

代入上式可得

\[
H^+(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)
\left[(1-\delta)\exp(i\omega(m+k)\Delta t)
+
\delta\exp(i\omega(m+k-1)\Delta t)
\right].
\]

因为

\[
\tau=(m-\delta)\Delta t,
\]

所以也可以分解为

\[
H^+(\omega,\tau)
=
\exp(i\omega\tau)
B^+(\omega,\tau)
C_{\rm lin}^+(\omega,\tau),
\]

其中

\[
B^+(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)\exp(i\omega k\Delta t),
\]

\[
C_{\rm lin}^+(\omega,\tau)
=
(1-\delta)\exp(i\omega\delta\Delta t)
+
\delta\exp[-i\omega(1-\delta)\Delta t].
\]

因此单端外推的频率域主公式为

\[
U_r(r,\omega)
=
\sum_{s\in\mathcal{A}(r)}
W(r,s)
H^+(\omega,\tau(s,r))
U_s(s,\omega),
\]

也就是

\[
U_r(r,\omega)
=
\sum_{s\in\mathcal{A}(r)}
W(r,s)
B^+(\omega,\tau)
C_{\rm lin}^+(\omega,\tau)
\exp[i\omega t(s,r)]
U_s(s,\omega).
\]

这里的主振荡相位就是

\[
\exp[i\omega t(s,r)].
\]

### 5.2 与当前时间域代码逐样点一致所需的边界修正

上面的 \(H^+\) 对应的是完整线性卷积。但是当前 C 代码中的 `kirdat_pick` 对每一个 \(k\) 都要求 `shift-k >= 0`。因此，完整线性卷积会比当前代码多出如下边界项：

\[
\frac{\delta}{\Delta t}F_k(\tau)x_0
\]

它出现在输出样点

\[
n=m+k-1,
\qquad k=0,1,\ldots,N-2.
\]

所以如果采用 FFT 实现线性卷积，必须在每个空间点对 \((s,r)\) 的卷积结果上执行修正：

\[
z_{m+k-1}
\leftarrow
z_{m+k-1}
-
\frac{\delta}{\Delta t}F_k(\tau)x_0,
\qquad k=0,1,\ldots,N-2,
\]

然后再乘以空间权重 \(W(r,s)\) 并累加到输出道中。若已经将空间权重乘入卷积结果，则修正应写为

\[
z_{m+k-1}
\leftarrow
z_{m+k-1}
-
W(r,s)\frac{\delta}{\Delta t}F_k(\tau)x_0.
\]

这一项不能省略。省略后，频率域结果会与当前时间域程序在每个走时对应的前沿附近相差若干样点。

---

## 6. FFTW 原生符号约定下的实现公式

FFTW 的常用正变换为

\[
X_\ell=\sum_{n=0}^{N_{\rm fft}-1}x_n\exp\left(-i\frac{2\pi n\ell}{N_{\rm fft}}\right),
\]

反变换为

\[
x_n=\frac{1}{N_{\rm fft}}\sum_{\ell=0}^{N_{\rm fft}-1}X_\ell
\exp\left(i\frac{2\pi n\ell}{N_{\rm fft}}\right).
\]

这与第 5 节的 \(\exp[i\omega t]\) 理论约定相反。因此，如果直接使用 FFTW 的 `fftwf_plan_dft_r2c_1d` 作为正变换，则应使用共轭形式：

\[
H^-(\omega,\tau)
=
\exp(-i\omega\tau)
B^-(\omega,\tau)
C_{\rm lin}^-(\omega,\tau),
\]

其中

\[
B^-(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)\exp(-i\omega k\Delta t),
\]

\[
C_{\rm lin}^-(\omega,\tau)
=
(1-\delta)\exp(-i\omega\delta\Delta t)
+
\delta\exp[i\omega(1-\delta)\Delta t].
\]

等价地，也可以直接用离散卷积核构造频率响应：

\[
H^-(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)
\left[(1-\delta)\exp[-i\omega(m+k)\Delta t]
+
\delta\exp[-i\omega(m+k-1)\Delta t]
\right].
\]

在实际 C/FFTW 实现中，推荐使用这个直接表达式，因为它最不容易写错。

---

## 7. 推荐的 FFTW 实现方式

为了与时间域代码完全一致，推荐使用“每个空间点对的线性卷积 + 边界修正 + 空间累加”的方式。

### 7.1 选择 FFT 长度

对给定走时 \(\tau\)，卷积核最大延迟为

\[
p_{\max}=m+N-2.
\]

为了避免 FFT 圆周卷积污染，应取

\[
N_{\rm fft}\geq N_t+p_{\max}.
\]

为了方便 FFTW 加速，可以取大于该值的 2 的整数次幂：

```c
int nfft = next_fast_size(nt + pmax);
```

若想简化实现，也可以对全程序取一个统一的

\[
N_{\rm fft}^{\rm global}
\geq
N_t+\max_{s,r}\left(m(s,r)+N-2\right),
\]

这样所有空间点对共用同一个 FFT 长度。

### 7.2 构造离散卷积核

对每一个 \(\tau\)，定义长度为 `nfft` 的实数组 `h`，初始化为 0。然后执行：

```c
q = tau / dt;
m = (int)ceilf(q);
delta = (float)m - q;

for (k = 0; k <= nsam-2; k++) {
    h[m+k]   += (1.0f-delta) * F[k] / dt;
    h[m+k-1] += delta        * F[k] / dt;
}
```

注意：

- 这里要求 `m >= 1`。正常物理走时 \(\tau>0\) 时一般成立；
- 如果可能出现 \(\tau=0\)，必须单独处理，否则 `m+k-1` 可能为负；
- `F[k]` 必须与原 `filt_set(tau)` 得到的 `filt[its][k]` 完全一致。

### 7.3 用 FFT 计算线性卷积

对输入道 `x[0:nt-1]` 做零填充到 `nfft`：

```c
for (i = 0; i < nfft; i++) {
    xpad[i] = (i < nt) ? x[i] : 0.0f;
    hpad[i] = h[i];
}
```

然后调用 FFTW：

```c
fftwf_plan px = fftwf_plan_dft_r2c_1d(nfft, xpad, X, FFTW_ESTIMATE);
fftwf_plan ph = fftwf_plan_dft_r2c_1d(nfft, hpad, H, FFTW_ESTIMATE);
fftwf_plan py = fftwf_plan_dft_c2r_1d(nfft, Y, ypad, FFTW_ESTIMATE);

fftwf_execute(px);
fftwf_execute(ph);

for (iw = 0; iw <= nfft/2; iw++) {
    float ar = X[iw][0];
    float ai = X[iw][1];
    float br = H[iw][0];
    float bi = H[iw][1];

    Y[iw][0] = ar*br - ai*bi;
    Y[iw][1] = ar*bi + ai*br;
}

fftwf_execute(py);

for (it = 0; it < nt; it++) {
    z[it] = ypad[it] / nfft;
}
```

这样得到的是完整线性卷积结果 \(z_n=(h_\tau*x)_n\)。

### 7.4 执行边界修正

为了与当前 C 代码严格一致，执行：

```c
x0 = x[0];
for (k = 0; k <= nsam-2; k++) {
    int nedge = m + k - 1;
    if (0 <= nedge && nedge < nt) {
        z[nedge] -= delta * F[k] * x0 / dt;
    }
}
```

然后将 `z[it]` 乘以空间权重并累加：

```c
for (it = 0; it < nt; it++) {
    out_trace[it] += W * z[it];
}
```

这样得到的 `out_trace` 才与原时间域 `kirdat_pick` 的离散行为一致。

---

## 8. 直接使用频率响应的实现方式

如果不显式构造 `h`，也可以直接构造频率响应。使用 FFTW 原生正变换时，令

\[
\omega_\ell=\frac{2\pi\ell}{N_{\rm fft}\Delta t},
\qquad \ell=0,1,\ldots,N_{\rm fft}/2.
\]

然后对每个频率计算

\[
H^-_\ell(\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)
\left[(1-\delta)\exp[-i\omega_\ell(m+k)\Delta t]
+
\delta\exp[-i\omega_\ell(m+k-1)\Delta t]
\right].
\]

再做

\[
Y_\ell=H^-_\ell X_\ell.
\]

最后 IFFT 得到 `z[it]`，再执行第 7.4 节的边界修正。这个方法与显式构造 `h` 再 FFT 是等价的。

如果要在理论说明中使用 \(\exp[i\omega t(s,r)]\)，则写成

\[
H^+(\omega,\tau)
=
B^+(\omega,\tau)C_{\rm lin}^+(\omega,\tau)
\exp[i\omega t(s,r)],
\]

但在 FFTW 原生代码中应使用符号相反的 \(H^-\)。

---

## 9. 双端外推的频率域实现流程

### 9.1 总体流程

频率域版本应保持原程序的外部行为不变：

```text
1. 读取 input_file，得到 tr_in[ns][nh][nt]
2. reverse_trace(nt, nh, ns, tr_in)
3. 读取 rgreen_file，得到 rtable
4. 读取 sgreen_file，得到 stable
5. filt_init(dt, length)
6. 检波点端外推：tr_in -> tr_out
7. 如果需要 interm，则写出 tr_out，保持与原程序一致
8. 将 tr_in 清零
9. 震源端外推：tr_out -> tr_in
10. 设置输出头信息
11. reverse_trace(nt, nh, ns, tr_in)
12. 写出 output_file
```

### 9.2 检波点端外推

对每个 `is` 和目标接收点 `ih`：

```text
for each input receiver ic in aperture:
    compute c, cc
    tau  = rtable[cc][c]
    dist = rdatum*rdatum + (ic-ih)*dh*(ic-ih)*dh
    W    = coef/M_PI * dh*rdatum*tau/dist
    build F[k] using filt_set-equivalent code
    compute z = exact_fft_time_operator(tr_in[is][ic], tau, F, dt, nt, nsam)
    tr_out[is][ih][:] += W * z[:]
```

数学形式为

\[
D_r(i_s,i_h,t_n)
=
\sum_{i_c\in\mathcal{A}(i_h)}
W_r(i_h,i_c)
\mathcal{T}_{\tau_r(i_c,i_h)}
\left[D_{\rm in}(i_s,i_c,t)\right]_n.
\]

### 9.3 震源端外推，`cmp=1`

`cmp=1` 时，接收点轴是相对震源的 offset。必须保留原代码中的 common-receiver 重排逻辑，尤其是

```c
hh = (r - ic*ds)/dh + 0.5;
```

外推形式为

\[
D_{\rm out}(i_s,i_h,t_n)
=
\sum_{i_c\in\mathcal{A}(i_s)}
W_s(i_s,i_c)
\mathcal{T}_{\tau_s(i_c,i_s)}
\left[D_r(i_c,hh,t)\right]_n.
\]

### 9.4 震源端外推，`cmp=0`

`cmp=0` 时，接收点轴已经是绝对坐标，不需要 `hh` 重排。外推形式为

\[
D_{\rm out}(i_s,i_h,t_n)
=
\sum_{i_c\in\mathcal{A}(i_s)}
W_s(i_s,i_c)
\mathcal{T}_{\tau_s(i_c,i_s)}
\left[D_r(i_c,i_h,t)\right]_n.
\]

---

## 10. 推荐封装的函数接口

建议将频率域时间算子封装成如下函数：

```c
void kirdat_pick_trace_freq(
    const float *trace_in,
    float *trace_out,
    int nt,
    float dt,
    float tau,
    float length,
    int nsam,
    int nfft,
    float *work_x,
    float *work_h,
    float *work_y,
    fftwf_complex *X,
    fftwf_complex *H,
    fftwf_complex *Y,
    fftwf_plan plan_x,
    fftwf_plan plan_h,
    fftwf_plan plan_y
);
```

该函数只完成

\[
\texttt{trace\_out}=\mathcal{T}_{\tau}[\texttt{trace\_in}],
\]

不乘空间权重。空间权重应在外层循环中乘入，以保持结构清晰。

函数内部必须执行：

```text
1. 根据 tau 构造 F[k]
2. 根据 m=ceil(tau/dt), delta=m-tau/dt 构造 h[p]
3. 对 trace_in 和 h 做零填充 FFT
4. 频率域相乘
5. IFFT，除以 nfft
6. 执行边界修正
7. 返回 trace_out[0:nt-1]
```

---

## 11. 精度和一致性检查

为了验证频率域版本是否正确，必须做以下测试。

### 11.1 单道单走时测试

固定一个输入道 `x` 和一个走时 `tau`，比较：

```text
time-domain: kirdat_pick loop
frequency-domain: kirdat_pick_trace_freq
```

误差应满足

\[
\frac{\|y_{\rm freq}-y_{\rm time}\|_2}{\|y_{\rm time}\|_2}<10^{-5}
\]

或达到单精度 FFT 可接受的误差水平。

### 11.2 单端外推测试

只做检波点端外推，比较 `tr_out`。注意原程序中 `tr_in` 已经经过一次 `reverse_trace`，测试时必须保持相同的数据方向。

### 11.3 双端完整测试

完整执行

\[
D_{\rm out}=\mathcal{D}_s\mathcal{D}_rD_{\rm in}.
\]

比较最终写出的 `output_file`。频率域版本必须保持：

- `nt`, `nh`, `ns` 不变；
- `dt`, `dh`, `ds` 不变；
- `o1`, `o2`, `o3` 与原程序一致；
- `cmp=0` 和 `cmp=1` 两个分支的道集组织方式不变；
- 最终输出前执行第二次 `reverse_trace`。

---

## 12. 不能省略或随意修改的部分

改写代码时，不要做以下改动：

1. 不要把 \(m\) 写成 `floor(tau/dt)`。当前代码等价于

   \[
   m=\left\lceil\frac{\tau}{\Delta t}\right\rceil.
   \]

2. 不要把 `delta` 写成 `tau/dt - floor(tau/dt)`。当前代码等价于

   \[
   \delta=m-\frac{\tau}{\Delta t}.
   \]

3. 不要省略线性插值响应 \(C_{\rm lin}\)。否则结果不等价。

4. 不要省略边界修正项。否则前沿样点与时间域代码不一致。

5. 不要改变 `reverse_trace` 的位置。

6. 不要改变 `cmp=1` 时的 common-receiver 重排逻辑。

7. 不要使用未零填充的 FFT 直接相乘。那会得到圆周卷积，而不是时间域代码对应的线性卷积。

8. 不要把 FFTW 原生负指数正变换与理论中的 \(\exp[i\omega t]\) 直接混用。若使用 FFTW 原生正变换，核函数应使用 \(H^-\)。

---

## 13. 最终应实现的单端频率域公式

理论写法采用正号相位约定：

\[
U_r(r,\omega)
=
\sum_{s\in\mathcal{A}(r)}
W(r,s)
B^+(\omega,\tau)
C_{\rm lin}^+(\omega,\tau)
\exp[i\omega t(s,r)]
U_s(s,\omega),
\]

其中

\[
B^+(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)\exp(i\omega k\Delta t),
\]

\[
C_{\rm lin}^+(\omega,\tau)
=
(1-\delta)\exp(i\omega\delta\Delta t)
+
\delta\exp[-i\omega(1-\delta)\Delta t],
\]

\[
\tau=t(s,r),
\qquad
m=\left\lceil\frac{\tau}{\Delta t}\right\rceil,
\qquad
\delta=m-\frac{\tau}{\Delta t}.
\]

实现时若采用 FFTW 原生正变换，则使用

\[
U_r(r,\omega)
=
\sum_{s\in\mathcal{A}(r)}
W(r,s)
H^-(\omega,\tau)
U_s(s,\omega),
\]

其中

\[
H^-(\omega,\tau)
=
\frac{1}{\Delta t}
\sum_{k=0}^{N-2}F_k(\tau)
\left[(1-\delta)\exp[-i\omega(m+k)\Delta t]
+
\delta\exp[-i\omega(m+k-1)\Delta t]
\right].
\]

IFFT 后必须执行边界修正

\[
z_{m+k-1}
\leftarrow
z_{m+k-1}
-
\frac{\delta}{\Delta t}F_k(\tau)x_0,
\qquad k=0,1,\ldots,N-2.
\]

然后再进行空间权重累加。

---

## 14. 最终应实现的双端公式

完整程序应实现

\[
D_{\rm out}(\omega)=\mathbf{K}_s(\omega)\mathbf{K}_r(\omega)D_{\rm in}(\omega),
\]

其中检波点端为

\[
D_r(i_s,i_h,\omega)
=
\sum_{i_c\in\mathcal{A}(i_h)}
K_r(i_h,i_c,\omega)D_{\rm in}(i_s,i_c,\omega),
\]

\[
K_r(i_h,i_c,\omega)
=
\frac{\texttt{coef}}{\pi}
\frac{\Delta h\,z_r\,\tau_r(i_c,i_h)}
{z_r^2+(i_c-i_h)^2\Delta h^2}
H(\omega,\tau_r(i_c,i_h)).
\]

震源端为

\[
D_{\rm out}(i_s,i_h,\omega)
=
\sum_{i_c\in\mathcal{A}(i_s)}
K_s(i_s,i_c,\omega)D_r(i_c,i_h,\omega),
\]

\[
K_s(i_s,i_c,\omega)
=
\frac{\texttt{coef}}{\pi}
\frac{\Delta s\,z_s\,\tau_s(i_c,i_s)}
{z_s^2+(i_c-i_s)^2\Delta s^2}
H(\omega,\tau_s(i_c,i_s)).
\]

这里 \(H\) 在理论说明中可以用 \(H^+\)，在 FFTW 原生实现中应使用 \(H^-\)。

---

## 15. 编译链接建议

如果使用单精度 FFTW：

```bash
gcc -O3 -fopenmp your_code.c -lfftw3f -lm -o your_program
```

如果使用 FFTW 线程接口：

```bash
gcc -O3 -fopenmp your_code.c -lfftw3f_threads -lfftw3f -lm -o your_program
```

初始化时需要：

```c
fftwf_init_threads();
fftwf_plan_with_nthreads(omp_get_max_threads());
```

但更安全的 OpenMP 写法是：每个线程使用自己的 FFTW 输入输出缓冲区和 plan，避免多个线程同时执行同一个 plan 并写同一缓冲区。

---

## 16. 简要结论

可以将当前时间域外推公式完全转换为频率域形式。理论上可以写成

\[
U_r(r,\omega)
=
\sum_s
A(r,s,\omega)
\exp[i\omega t(s,r)]
U_s(s,\omega),
\]

但为了与当前代码完全一致，必须令

\[
A(r,s,\omega)=W(r,s)B^+(\omega,\tau)C_{\rm lin}^+(\omega,\tau),
\]

并在 FFT 实现中保留当前 `kirdat_pick` 由 `shift-isam >= 0` 引入的边界修正项。否则只能得到近似等价的频率域外推，而不是与原时间域代码逐样点一致的实现。
