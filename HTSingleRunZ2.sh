#!/bin/bash

cxxfile=$1
pyfile=$2
id=$3

./cxSMThermalMassOnly/cxSMCritical_Scan_Z2_HT.x $cxxfile $id
python cxSMScannerFull.py -i $cxxfile -o $pyfile