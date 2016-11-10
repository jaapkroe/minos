#!/bin/bash
#for n in 1000 2000 4000 8000 16000 32000 64000 128000 ; do
for m in `seq 60`; do
for i in `seq 1 3 90` 100 115 130 150 170 190 216; do 
  n=$((i**3)); 
  printf "%-12d" $n >&2
  ../write-config -n $n 
  minos -q config.xyz; 
  #read x
done 2>&1 | tee -a scaling.dat
done
