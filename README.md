# Kirch_datuming

## `kirchmig2d_cig` 两侧直达波/递归噪声切除

在生成角道集（`adj=1`）前，可以对每个炮集启用关于震源位置左右对称的平滑时变切除：

```text
direct_mute=1
direct_mute_slope=0.0005
direct_mute_intercept=0.02
direct_mute_taper=0.05
```

默认切除边界为
`t = direct_mute_intercept + direct_mute_slope * abs(receiver_x - source_x)`；
设置 `direct_mute_invert=1` 后改为减号，从而形成尖点位于震源位置的倒 V 形。
边界之前的数据置零，边界之后在 `direct_mute_taper` 秒内使用升余弦从 0
平滑恢复到 1，因此不会产生硬截断。`cmp=1` 时接收点位置按 `source + offset`
计算；`cmp=0` 时第二维坐标直接作为接收点位置。该功能默认关闭，且
`direct_mute_taper` 必须大于零。

## `kirchmig2d_cig` 内存策略

源、接收点走时表及其导数现在与 `kirchmig2d` 一样按 `(z,x)` 切片读取，
每个表只缓存插值所需的两个相邻切片，不再把完整的三维走时体载入内存。
可选的 `sgradx/sgradz/rgradx/rgradz` 也使用一对共享切片缓冲区流式读取。
这一修改不改变走时插值、角度计算或 Kirchhoff 累加顺序；主要常驻内存由
与炮点/接收点数量成正比降为与 `nz*nx` 成正比。完整角道集
`nz*nx*angle_n` 仍保留在内存中。

## 五层盐丘 ADCIG 脚本

`run_salt_adcig.sh` 为五个递归层分别设置平滑直达波切除参数。当前使用一套
以震源为尖点的 V 形切除预设：五层均设置 `direct_mute_invert=0`、截距为
`0 s`，斜率为 `50,40,30,20,10 s/km`，缓冲长度为
`0.05,0.06,0.07,0.08,0.10 s`。由于模型横向坐标使用 km，第一层在距离震源
`0.01 km`（10 m）处的切除时间正好为 `0.5 s`；距离越远，切除时间越大，
并将小于该边界时间的全部能量去除。斜率随深度减小使深层切除线更平；
五组数组 `direct_mute_enable/`
`direct_mute_invert/direct_mute_slope/`
`direct_mute_intercept/direct_mute_taper` 均可在脚本开头独立调整。运行方式：

```text
bash run_salt_adcig.sh [横向采样序号]
```

每层 ADCIG 计算完成后，同一个 `run_salt_adcig.sh` 会继续检查五个输入角度轴
是否一致，自动保留可用的 0--60 度范围，合成完整三维 ADCIG，并输出指定位置
的原始深度-角度道集和仅沿角度方向轻微平滑的显示版本；不再另外启动第二个
脚本。
