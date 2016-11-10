#!/bin/bash

awk '{ if(NF==7) { if(M[$1]=="")M[$1]=1; t[$1][M[$1]]=$3; M[$1]+=1; } }
  END {
    for (n in t) {
      for (m in t[n]) {
        tmp[m] = t[n][m]
      }
      asort(tmp)
      for (m=1;m<=length(tmp)*0.20;m++) {
        print n,tmp[m]
      }
      delete tmp
    }
  }' scaling.dat | sort -nk1 -k2.2
