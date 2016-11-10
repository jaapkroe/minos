#!/bin/bash

awk '{ if(NF==7) { if(t[$1]=="" || $3<t[$1]) t[$1]=$3 } }
  END {
    for (N in t) {
      printf "%-9d %10.4f\n",N,t[N]
    }
  }' scaling.dat | sort -n
