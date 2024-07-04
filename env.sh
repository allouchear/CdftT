#!/bin/bash
# SCRIPT=$(readlink -f $0) only for linux
# Works on both linux and Mac OS X.

function read_link() {
    local path=$1
        if [ -d $path ] ; then
        local abspath=$(cd $path; pwd)
        else
            local prefix=$(cd $(dirname -- $path) ; pwd)
            local suffix=$(basename $path)
            local abspath="$prefix/$suffix"
        fi
    echo $abspath
}
SCRIPT=$(read_link $0)
export  WORKDIR=`dirname $SCRIPT`
echo WORKDIR=$WORKDIR
export LIBCDFTTDIR=$WORKDIR
export LIBCDFTTSRC=$LIBCDFTTDIR/src/applications
export LIBCDFTT=$LIBCDFTTDIR/lib/libcdftt.so
#export PYTHONPATH=$LIBCDFTTDIR/python/src/Utils:$LIBCDFTTDIR/python/src/Molecule
echo $LIBCDFTTSRC
#source /home/allouche/Softwares/nvhpc/env.sh
