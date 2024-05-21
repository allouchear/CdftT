#!/bin/bash
cd ../../..
source env.sh
cd -
cd ../../lib
./compcdfttlib.sh
cd -
rm cdftt
make
