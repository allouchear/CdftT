#!/bin/bash
cd ../../..
source env.sh
cd -
cd ../../lib
./cleancdfttlib.sh
cd -
make clean
