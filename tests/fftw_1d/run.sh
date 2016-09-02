#!/bin/bash

# compile
make
# generate data (see script for differn
./fft_1d-data_generatation.sh
# make the fourier transform
./fft_1d > fft_1d.out

# generate real space data from fourier coefficients
awk 'BEGIN{N=1000; pi2=2*atan2(0,-1); n=0} {
  a[n]=$2; b[n]=$3; n++
  } END{
    for(i=0;i<N;i++) { 
      x=pi2*i/N; 
      f=a[0]; 
# real part of
#       sum_q a_q exp(-i q x) 
#                 = q_0 + a_1 * cos(x) + a_1 cos(-x) + b_1 sin(x) + (-b_1) sin(-x) 
#                       + ...
#                 = q_0 + 2 a_1 cos(x) + 2 b_1 sin(x) + ...
      for(j=1;j<n;j++) f += 2* ( a[j]*cos(j*x) - b[j]*sin(j*x) ); 
      print x,f 
    }
  }' fft_1d.out > fft_1d.gen

# plot results
gnuplot5 fft_1d.gpl
