#!/bin/bash

id=1
dir=$(pwd)
cd cxSMThermalMassOnly
make clean; make
cd $dir
datadir=$dir/Data
mkdir -p $datadir
mkdir -p $datadir/cxx
mkdir -p $datadir/cosmo
while [ $id -le 10 ]
  do
    printf -v ID "%02d" ${id}
    screen -S CXSM0828"$ID" -d -m bash -c "./cxSMSingleRun.sh $datadir "$id"; exec bash;"
    id=$[$id+1]
  done
