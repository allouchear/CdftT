#!/bin/bash
cd /home/lgardre/Softwares/cdftt/04-enePC_20_512_0
export OMP_NUM_THREADS=128
export OUT=/home/lgardre/Softwares/cdftt/04-enePC_20_512_0/xcdftt_enePC_20_512_0.out
echo "DIR=/home/lgardre/Softwares/cdftt/04-enePC_20_512_0" > $OUT

chmod u+x xcdftt_enePC_20_512_0

time ./xcdftt_enePC_20_512_0 128 >  $OUT
