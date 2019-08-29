#!/bin/bash

runs=$1
dir=$(pwd)
cd cxSMThermalMassOnly
make clean; make cxSMCritical_Scan_HT.x
cd $dir

sed 's/__TIMES__/'$runs'/g' HTjob.sub.tmp > HTjob.sub
condor_submit HTjob.sub