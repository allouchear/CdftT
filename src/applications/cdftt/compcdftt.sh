#!/bin/bash

# Help message
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo "Usage: $0 [target]"
    echo "Compile the CDFTT library for the specified target."
    echo "Targets:"
    echo "  all     Compile for all targets (default)"
    echo "  profile Compile for profile target"
    exit 0
fi

cd ../../..
source env.sh
cd -
cd ../../lib
./cleancdfttlib.sh
./compcdfttlib.sh $1
cd -
./cleancdftt.sh
make $1
