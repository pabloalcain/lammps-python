LDFLAGS=-L/home/pablo/lammps/src
LDLIBS=-llammps -lfftw3 -lmpi -lpthread
INFLAGS=-I/home/pablo/lammps/src
CXXFLAGS=-Wno-write-strings -g -O0 -fPIC

all: simple.o gofr.o minkowski.o
	mpic++ $(LDFLAGS) simple.o gofr.o minkowski.o $(LDLIBS) -o simple.e

library: gofr.o minkowski.o
	mpic++ $(LDFLAGS) gofr.o minkowski.o $(LDLIBS) -shared -o ../libanalysis.so 

simple.o: simple.cpp gofr.h minkowski.h
	mpic++ $(INFLAGS) $(CXXFLAGS) -c simple.cpp

gofr.o: gofr.cpp gofr.h
	mpic++ $(INFLAGS) $(CXXFLAGS) -Wall -c gofr.cpp

minkowski.o: minkowski.cpp minkowski.h
	mpic++ $(INFLAGS) $(CXXFLAGS) -Wall -c minkowski.cpp

clean:
	rm *.o simple.e *.so
