#!/bin/bash
cd ../..
source env.sh
cd -
make dep
make -j 10 $1