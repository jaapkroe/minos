# NOTE : requirements: 
# 	o lemon graph library
# 	o a recent compiler (c++11) support

CPPFLAGS=-O3 -g -std=c++11 -Wall #-fdiagnostics-color=auto 
#CPPFLAGS+=-DLEMON
CPPFLAGS+=-I/vol/tcm00/jaap/prog/opt/include -L/vol/tcm00/jaap/prog/opt/lib

#CPPFLAGS+=-fopenmp 
LDLIBS=-lm -lnfft3 -lfftw3

nfft_1d: nfft_1d.cpp
example: example.cpp

clean:
	rm -f nfft_1d example
