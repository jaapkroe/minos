# NOTE : requirements: 
# 	o lemon graph library
# 	o a recent compiler (c++11) support

CPPFLAGS=-O3 -g -std=c++11 -Wall #-fdiagnostics-color=auto 
#CPPFLAGS+=-DHASLEMON
#CPPFLAGS+=-fopenmp 
LDLIBS=-lm -lrt

minos: minos.cpp

clean:
	rm -f minos
