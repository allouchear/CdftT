#!/bin/bash
cd /home/lgardre/Softwares/cdftt/06-enePC_16_256_0
export OMP_NUM_THREADS=32
export OUT=/home/lgardre/Softwares/cdftt/06-enePC_16_256_0/xcdftt_enePC_16_256_0.out
echo "DIR=/home/lgardre/Softwares/cdftt/06-enePC_16_256_0" > $OUT

chmod u+x xcdftt_enePC_16_256_0

time ./xcdftt_enePC_16_256_0 32 >  $OUT
