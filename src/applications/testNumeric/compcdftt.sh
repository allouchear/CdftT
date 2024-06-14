#!/bin/bash
cd ../../..
source env.sh
cd -
cd ../../lib
#./cleancdfttlib.sh
./compcdfttlib.sh
cd -
./cleancdftt.sh
make
