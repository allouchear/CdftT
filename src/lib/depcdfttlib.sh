#!/bin/bash
cd ../..
source env.sh
cd -
rm */Dep.mk
make dep
