#!/bin/bash
cd /home/lgardre/Softwares/cdftt/03-enePC_20_256_0
export OMP_NUM_THREADS=64
export OUT=/home/lgardre/Softwares/cdftt/03-enePC_20_256_0/xcdftt_enePC_20_256_0.out
echo "DIR=/home/lgardre/Softwares/cdftt/03-enePC_20_256_0" > $OUT

chmod u+x xcdftt_enePC_20_256_0

time ./xcdftt_enePC_20_256_0 64 >  $OUT
