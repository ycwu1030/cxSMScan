#!/bin/bash

runs=$1
dir=$(pwd)
cd cxSMThermalMassOnly
make clean; make cxSMCritical_Scan_Z2_HT.x
cd $dir

sed 's/__TIMES__/'$runs'/g' HTjob_Z2.sub.tmp > HTjob_Z2.sub
condor_submit HTjob_Z2.sub