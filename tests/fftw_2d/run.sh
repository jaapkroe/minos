#!/bin/bash

# compile
make
# generate data (see script for differn
./fft_2d-data_generatation.sh
# make the fourier transform
./fft_2d > fft_2d.out
# plot results
gnuplot fft_2d.gpl
