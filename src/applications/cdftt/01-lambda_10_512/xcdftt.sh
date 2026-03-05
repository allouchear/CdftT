#!/bin/bash
cd /home/lgardre/Softwares/cdftt
export OMP_NUM_THREADS=64
export OUT=/home/lgardre/Softwares/cdftt/xcdftt.out
echo "DIR=/home/lgardre/Softwares/cdftt" > $OUT

chmod u+x xcdftt

time ./xcdftt 64 >  $OUT
