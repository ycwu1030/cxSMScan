#!/bin/bash

id=1
while [ $id -le 10 ]
  do
    printf -v ID "%02d" ${id}
    screen -X -S CXSM0828"$ID" quit
    id=$[$id+1]
  done
