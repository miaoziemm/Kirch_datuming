cd cyq1

# migration reference
ns=768
dx=0.01
os=0

# travel time
../build/bin/shotgen shotfile=sht.rsf nshot=$ns oy=$os dy=$dx
../build/bin/shotgen shotfile=rcv.rsf nshot=$ns oy=$os dy=$dx

../build/bin/eikods2d in=marmsmooth.rsf out=tts.rsf shotfile=sht.rsf tds1=tds0s.rsf tdl1=tdl0s.rsf
../build/bin/eikods2d in=marmsmooth.rsf out=ttr.rsf shotfile=rcv.rsf  tds1=tds0r.rsf tdl1=tdl0r.rsf

# migration
../build/bin/kirchmig2d seismic_data=data.rsf migration=mig_ref.rsf stable=tts.rsf sderiv=tds0s.rsf rtable=ttr.rsf rderiv=tds0r.rsf




