# Seismic_Engine
A Seismic Toolkit for forward or migration

TODO LIST
-  Add a common shot data convert tool for different seismic data format
-  IO optimization
-  完善文件系统，支持多种数据格式，考虑是否需要将文件系统与SEDATA合并，在文件系统中直接将基本数据进行读写
-  SEDATA中只需要保存数据处理系统即可
-  文件系统中函数进行分类，与json相关的放在json解析中，与数据解析相关的（比如su数据，segy数据解析）放在数据解析中
-  思考文件头中还需要什么信息（文件名，大小，路径，时间戳，类型）
-  道集转换，比如共炮点，共检波点，共中心点，共偏移距
-  

基本函数库
内存分配，json解析，通用结构体

数学函数库
包含eigen，以及打包为通用的seimicmath库，包含常用的数学函数，矩阵运算，向量运算等，以及fft运算

文件系统架构
架构
- 将数据分解成基本的处理元素（主动源：共炮，共中心，共检波，共偏移距，超道集，叠后数据，被动源）database类
- 数据格式 SU，SEGY，SE
- 从文件中加载seismicdata类 seismicdata.read() seismicdata.segyread() 写入文件seismicdata.write()。。
seimicdata.get() 获取道集数据 seismicdata.get()，获取道集数据
seismicdata.set() 设置道集数据 seismicdata.set()，设置道集数据
转换databse类型 seismicdata.convert()，转换为不同的道集类型
生成简单道集功能 seismicdata.generate()，根据输入的参数生成一个简单的道集

- 基本的模型元素（速度，品质因子，各向异性参数）model类
从文件中加载model类 model.read() model.segyread() 写入文件model.write()。。
生成一个初始模型 model.generate()，根据输入的参数生成一个初始模型

- 结果/中间元素（射线路径，走时表，结果数据  时间、深度域）
processdata类


- 更高级的处理工具
- 正演工具（forward类）继承model类和database类，可以进行dataset的生成与保存
正演工具分为两类 raybased和wavebased，raybased使用射线追踪方法进行正演，wavebased使用波动方程进行正演
不同的正演方式可能输出的内容不一致，这个思考一下是否需要统一，是否需要按照不同的正演方式进行类的定义
- 反演工具（inversion类）继承model类和database类，可以进行model的生成与保存

- 第二种定义方法
- 射线类 ray namespace 包含正演偏移和反演类 SERAY::forward 类 SERAY::inversion 类 SERAY::migration
- 波动方程类 wave namespace 包含正演偏移和反演类 SEWAVE::forward 类 SEWAVE::inversion 类 SEWAVE::migration
- 走时类 traveltime namespace 包含正演偏移和反演类 SETRAVELTIME::forward 类 SETRAVELTIME::inversion 类 SETRAVELTIME::migration

- 文件类型
- 基本的处理元素，每次保存时生成两个文件，第一是se的头文件，第二是se@的数据文件夹
- 基本的模型元素，每次保存时生成两个文件，第一是se的头文件，第二是se@的数据文件夹
  模型的data数据文件中，按顺序保存各类模型，地表，速度，品质因子，各向异性参数等
  数据保存格式类似rsf的格式，一个是数据描述文件，一个保存纯数据文件