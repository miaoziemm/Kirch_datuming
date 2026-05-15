#将marmsmooth.rsf分层
nlayers=5
overlap=10
model_file=marmsmooth.rsf
data_file=data.rsf

nx=768
nz=261
dx=0.01
dz=0.01


../build/bin/layervel input_file=${model_file} nlayers=${nlayers} overlap=${overlap} output_file_base=marmlayer.rsf

# 对第一层进行走时计算和成像
../build/bin/shotgen shotfile=sht_1.rsf nshot=${nx} oy=0 dy=${dx} 
../build/bin/shotgen shotfile=rcv_1.rsf nshot=${nx} oy=0 dy=${dx} 

../build/bin/eikods2d in=marmlayer_1.rsf out=time1s.rsf shotfile=sht_1.rsf tdl1=tdl1s.rsf tds1=tds1s.rsf b1=2 b2=2
# ../build/bin/eikods2d in=marmlayer_1.rsf out=time1r.rsf shotfile=rcv_1.rsf tdl1=tdl1r.rsf tds1=tds1r.rsf b1=2 b2=2
cp time1s.rsf time1r.rsf
cp tds1s.rsf tds1r.rsf
cp tdl1s.rsf tdl1r.rsf

../build/bin/kirchmig2d seismic_data=${data_file} migration=mig_1.rsf stable=time1s.rsf sderiv=tds1s.rsf rtable=time1r.rsf rderiv=tds1r.rsf

../build/bin/greensol_auto ttabel_file=time1s.rsf tgreen_file=time1s_green.rsf model_file=marmlayer_1.rsf
../build/bin/greensol_auto ttabel_file=time1r.rsf tgreen_file=time1r_green.rsf model_file=marmlayer_1.rsf

../build/bin/kirchdat2d_auto input_file=${data_file} output_file=rdata_1.rsf aperture=300 taper=0 length=0.05 sgreen_file=time1s_green.rsf rgreen_file=time1s_green.rsf model_file=marmlayer_1.rsf


# 从第二层开始循环
for i in $(seq 2 $nlayers)
do
    # 打印一下当前处理的层数
    echo "Processing layer $i"
    # 生成当前层的shot和rcv文件
    ../build/bin/shotgen shotfile=sht_${i}.rsf nshot=${nx} oy=0 dy=${dx} 
    ../build/bin/shotgen shotfile=rcv_${i}.rsf nshot=${nx} oy=0 dy=${dx} 
    # 计算走时并成像
    ../build/bin/eikods2d in=marmlayer_${i}.rsf out=time${i}s.rsf shotfile=sht_${i}.rsf tdl1=tdl${i}s.rsf tds1=tds${i}s.rsf b1=2 b2=2
    # ../build/bin/eikods2d in=marmlayer_${i}.rsf out=time${i}r.rsf shotfile=rcv_${i}.rsf tdl1=tdl${i}r.rsf tds1=tds${i}r.rsf b1=2 b2=2
    cp time${i}s.rsf time${i}r.rsf
    cp tds${i}s.rsf tds${i}r.rsf
    cp tdl${i}s.rsf tdl${i}r.rsf

    ../build/bin/kirchmig2d seismic_data=rdata_$(($i-1)).rsf migration=mig_${i}.rsf stable=time${i}s.rsf sderiv=tds${i}s.rsf rtable=time${i}r.rsf rderiv=tds${i}r.rsf
    ../build/bin/greensol_auto ttabel_file=time${i}s.rsf tgreen_file=time${i}s_green.rsf model_file=marmlayer_${i}.rsf
    ../build/bin/greensol_auto ttabel_file=time${i}r.rsf tgreen_file=time${i}r_green.rsf model_file=marmlayer_${i}.rsf

    # 外推到下一层
    ../build/bin/kirchdat2d_auto input_file=rdata_$(($i-1)).rsf output_file=rdata_${i}.rsf aperture=300 taper=0 length=0.05 sgreen_file=time${i}s_green.rsf rgreen_file=time${i}r_green.rsf model_file=marmlayer_${i}.rsf
done

# 最后一层走时计算加成像
../build/bin/shotgen shotfile=sht_${nlayers}.rsf nshot=${nx} oy=0 dy=${dx} 
../build/bin/shotgen shotfile=rcv_${nlayers}.rsf nshot=${nx} oy=0 dy=${dx} 

../build/bin/eikods2d in=marmlayer_${nlayers}.rsf out=time${nlayers}s.rsf shotfile=sht_${nlayers}.rsf tdl1=tdl${nlayers}s.rsf tds1=tds${nlayers}s.rsf b1=2 b2=2
    # ../build/bin/eikods2d in=marmlayer_${nlayers}.rsf out=time${nlayers}r.rsf shotfile=rcv_${nlayers}.rsf tdl1=tdl${nlayers}r.rsf tds1=tds${nlayers}r.rsf b1=2 b2=2
cp time${nlayers}s.rsf time${nlayers}r.rsf
cp tds${nlayers}s.rsf tds${nlayers}r.rsf
cp tdl${nlayers}s.rsf tdl${nlayers}r.rsf

../build/bin/kirchmig2d seismic_data=rdata_$(($nlayers-1)).rsf migration=mig_${nlayers}.rsf stable=time${nlayers}s.rsf sderiv=tds${nlayers}s.rsf rtable=time${nlayers}r.rsf rderiv=tds${nlayers}r.rsf


../build/bin/layercom full_model=marmsmooth.rsf layer_model_base=marmlayer.rsf layer_image_base=mig.rsf nlayers=${nlayers} output_file=mig_full.rsf
