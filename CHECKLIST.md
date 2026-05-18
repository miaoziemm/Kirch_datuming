# Kirchdat2d 频率域重构验收清单（对照 `src/butterfly_ready_frequency_refactor_guide.md`）

> 目的：把指导文档中的“硬约束”映射到当前实现中的具体代码位置，便于后续改动快速回归验收。

## 1) 外部参数与使用方式

- [x] 保持参数名与主语义：`input_file/output_file/sgreen_file/rgreen_file/model_file/verb/cmp/aperture/taper/length/interm`。  
  代码：`src/kirchdat2d_frequency_cs.cpp` 参数读取段。  
  位置：`L282-L287`, `L293-L303`, `L344-L347`
- [x] `should_datum` 仅从 `model_file` header 读取（无命令行兜底扩展）。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L297-L299`

## 2) 输入/输出维度与头信息

- [x] 输入保持 3D `(t,h,s)`。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L305-L313`
- [x] 输出保持 3D，且头信息与原实现一致（`n/d/o`）。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L382-L391`

## 3) 时间反转与 interm 行为

- [x] 首次反转：读入后立即 `reverse_trace`。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L320-L322`
- [x] 二次反转：最终写出前 `reverse_trace`。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L393-L394`
- [x] `interm`：检波点端外推后 IFFT 写出，不额外反转。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L371-L374`

## 4) 计算结构（两阶段外推）

- [x] 检波点端外推（receiver-side）先执行。  
  代码：`apply_receiver_operator_freq_iw` 调用。  
  位置：`L367-L369`
- [x] 震源端外推（source-side）后执行。  
  代码：`apply_source_operator_freq_iw` 调用。  
  位置：`L376-L378`
- [x] `cmp=1` 保留 `ir/r/sleft/sright/ih/hh` 映射逻辑（未简化为固定 `ih`）。  
  代码：`apply_source_operator_freq_iw` 的 `cmp==1` 分支。  
  位置：`L183-L226`

## 5) 与时域一致性的关键点

- [x] `taper_weight` 与时域整数除法语义对齐（包括 `cmp=1` 场景阈值行为）。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L103-L114`
- [x] `tap<=0` 时等价“无 taper”（返回 1）。  
  代码：`src/kirchdat2d_frequency_cs.cpp`。  
  位置：`L105`

## 6) 辅助工具（测试链路）

- [x] 合成数据生成工具 `rsf_synth_gen` 已加入构建。  
  代码：`src/CMakeLists.txt`, `src/rsf_synth_gen.cpp`。  
  位置：`src/CMakeLists.txt L56-L57, L87-L88, L106-L107, L124-L125`
- [x] 误差比较工具 `rsf_l2cmp` 已加入构建。  
  代码：`src/CMakeLists.txt`, `src/rsf_l2cmp.cpp`。  
  位置：`src/CMakeLists.txt L56-L57, L87-L88, L106-L107, L124-L125`

---

## 推荐回归命令（最小集合）

1. 编译：
   ```bash
   cmake --build build -j4
   ```
2. 端到端对比（示例）：
   - 生成：`rsf_synth_gen ...`
   - 时域：`kirchdat2d_auto_cs ...`
   - 频域：`kirchdat2d_freq ...`
   - 误差：`rsf_l2cmp a=... b=...`

