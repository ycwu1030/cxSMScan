#!/bin/bash

dir=$1
rid=$2
id=1
while [ $id -le 10 ]
  do
     kid=$[$rid*100+$id]
     ./cxSMThermalMassOnly/cxSMCritical_Scan.x $dir/cxx $kid
     infile=$dir/cxx/cxSMPoints_$kid.dat
     outfile=$dir/cosmo/cxSMPoints_$kid.dat
     python cxSMScannerFull.py -i $infile -o $outfile
     id=$[$id+1]
  done
