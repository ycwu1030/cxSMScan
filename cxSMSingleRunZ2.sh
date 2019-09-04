#!/bin/bash

dir=$1
rid=$2
id=1
while [ $id -le 10 ]
  do
     kid=$[$rid*100+$id]
     ./cxSMThermalMassOnly/cxSMCritical_Scan_Z2.x $dir/cxx $kid
     infile=$dir/cxx/cxSMPointsZ2_$kid.dat
     outfile=$dir/cosmo/cxSMPointsZ2_$kid.dat
     python cxSMScannerFullZ2.py -i $infile -o $outfile
     id=$[$id+1]
  done
