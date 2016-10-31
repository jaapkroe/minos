# NOTE : requirements: 
# 	o lemon graph library
# 	o a recent compiler (c++11) support

CPPFLAGS=-O3 -std=c++11 -Wall #-fdiagnostics-color=auto 
#CPPFLAGS+=-DHASLEMON
#CPPFLAGS+=-fopenmp 
LDLIBS=-lm -lfftw3

minos: minos.cpp

clean:
	rm -f minos
