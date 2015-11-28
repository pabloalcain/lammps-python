CC=gcc
CFLAGS=-O3 -std=gnu99 -fPIC
LDFLAGS=-L$(LAMMPS)/src
LDLIBS=-llammps -lmpi -lpthread
INFLAGS=-I$(LAMMPS)/src
CXXFLAGS=-Wno-write-strings -g -O3 -fPIC


library: minkowski.o rdf.o
	mpic++ $(LDFLAGS) $^ $(LDLIBS) -shared -o libanalysis.so

install: library
	cp libanalysis.so /usr/local/lib

simple: simple.o rdf.o minkowski.o
	mpic++ $(LDFLAGS) $^ $(LDLIBS) -o simple.e

simple.o: simple.cpp simple.h
	mpic++ $(INFLAGS) $(CXXFLAGS) -c $<

rdf.o: rdf.c rdf.h
	$(CC) $(CFLAGS) -c $< -o $@

minkowski.o: minkowski.cpp minkowski.h
	mpic++ $(INFLAGS) $(CXXFLAGS) -Wall -c minkowski.cpp

clean:
	rm -f *.o simple.e *.so
