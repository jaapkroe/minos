#!/bin/bash
#for n in 1000 2000 4000 8000 16000 32000 64000 128000 ; do
for s in `seq 0 1 24` ; do 
  n=$(echo "2^($s)" | bc ) ; 
  printf "%-12d" $n >&2
  ./write-config -n $n 
  minos -q config.xyz; 
  #read x
done 2>&1 | tee scaling.dat

