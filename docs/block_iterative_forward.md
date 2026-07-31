# FMM 初场分块迭代正演

`block_iterative_forward` 是一个串行二维声学正演程序。程序沿深度方向建立等厚核心区和等宽重叠区，调用 `reckirch` 中已有的 `efmm_*` 接口计算首波走时，以 FMM 射线场初始化下部区域，再用 FD8 局部求解和 Robin 特征条件执行乘法型向下—向上扫描。

## 算法流程

1. 使用 `sewave::read_rsf2d` 读取 RSF 速度模型。
2. 将深度网格等分为 `blocks` 个核心区，各区域上下分别扩展 `overlap` 个网格。`overlap` 是单侧宽度。
3. 对每炮调用 `efmm_init`、`efmm_set_vel` 和 `efmm_solver`，获得全模型首波走时 \(T(\boldsymbol{x})\)。
4. 采用二维几何扩散近似构造初始射线场：

   \[
   u^{\mathrm{ray}}(\boldsymbol{x},t)
   =
   a_{\mathrm{ray}}
   \sqrt{\frac{h}{\max(v_sT(\boldsymbol{x}),h)}}
   w\!\left(t-T(\boldsymbol{x})\right),
   \]

   其中 \(h=\min(\Delta x,\Delta z)\)，\(v_s\) 为震源点速度，默认 `ray_scale=1`。
5. 第一区域用真实震源进行一次局部 FD8 正演；下部区域保留 FMM 射线界面历史作为迭代初值。
6. 向下扫描使用

   \[
   \mathcal{B}^-u=u_t-\beta c u_z
   \]

   的相邻区域差值驱动下部区域。每个下部区域只读取已经更新的上部区域，因此此阶段不会提前引入深部反射。
7. 向上扫描使用

   \[
   \mathcal{B}^+u=u_t+\beta c u_z
   \]

   的相邻区域差值，将深部反射和散射逐块返回到地表区域。
8. 从第一区域的检波深度提取地震记录，使用 `sewave::write_rsf3d` 输出。

局部区域内部采用时间二阶、空间八阶有限差分。人工界面附近的三个网格使用空间二阶模板，使单行 Robin 数据能够进入 FD8 内部；区域核心区仍采用 FD8。左右边界以及物理顶、底边界的吸收公式与 `sewave2d` 保持一致。基准脚本固定 `flag_smooth=0` 和 `type_compute_Laplace=0`，以避免速度光滑或伪谱核造成不一致。

## 构建

```bash
cmake -S . -B build -DSE_USE_MPI=OFF
cmake --build build --target block_iterative_forward sewave2d rsf_l2cmp -j
```

## 单次运行

```bash
build/bin/block_iterative_forward \
  velocity=vel.rsf output=block_record.rsf \
  nt=3210 dt=0.0005 fdom=15 \
  sx=2500 sz=190 rz=190 \
  blocks=4 overlap=22 cycles=1 \
  beta=1.0 gate=0.003 \
  nbc=40 L=30 alpha=1 \
  ray_output=fmm_ray_record.rsf \
  traveltime=fmm_time.rsf
```

`traveltime=` 当前只支持 `ns=1`；多炮正演可以省略该参数。多炮在炮循环中串行执行。

## 与 sewave2d 串行比较

```bash
NT=3210 DT=0.0005 FDOM=15 \
SX=2500 SZ=190 RZ=190 \
REPEATS=3 \
scripts/benchmark_block_iterative.sh vel.rsf
```

脚本执行以下控制：

- `OMP_NUM_THREADS=1`，并将常见数学库线程数固定为 1；
- `sewave2d` 使用声学 FD8、`flag_smooth=0`；
- 两个程序采用相同的震源、时间采样、吸收边界和检波几何；
- 使用现有 `rsf_l2cmp` 输出 `rel_l2` 与 `max_abs`；
- GNU `time` 可用时，同时记录最大 RSS；否则 RSS 写为 `-1`。

输出目录中包含 `benchmark_summary.txt`、`benchmark_metrics.csv`、两种地震记录、FMM 射线初始记录和逐次运行日志。

## 已验证配置

独立 C++ 同核测试采用 174×500 Marmousi 模型、\(\Delta x=\Delta z=10\) m、`nt=3210`、`dt=0.0005` s、`fdom=15` Hz、单炮 500 道记录。默认四块配置的相对 L2 记录误差为 0.947%，串行时间约为同核全域 FD8 的 2.15 倍。

一次四块双扫包含第一区域震源传播、3 次向下校正和 3 次向上校正，共 7 次局部传播。因此该方法当前用于验证分块射线—波动校正和后续并行/多级算法研究；在串行条件下，不应预期它快于一次全域显式 FD8 正演。
