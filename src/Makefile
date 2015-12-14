CSRC=$(wildcard *.c)
CXXSRC=$(wildcard *.cpp)
COBJ=$(patsubst %.c,%_c.o,$(CSRC))
CXXOBJ=$(patsubst %.cpp,%_cpp.o,$(CXXSRC))

LIB=libanalysis.so
EXE=simple.e


CC=gcc
CXX=mpic++
LD=mpic++
CFLAGS=-g -O3 -fPIC -std=gnu99 -fopenmp
LDFLAGS=-L$(LAMMPS)/src -lgomp
LDLIBS=#-llammps
CXXFLAGS=-Wno-write-strings -g -O3 -fPIC -I$(LAMMPS)/src

default: $(LIB)

executable: $(EXE)

library: $(LIB)

all: $(LIB) $(EXE)

install: $(LIB)
	cp $< /usr/local/lib

$(LIB): $(COBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -shared -o $@

$(EXE): $(CXXOBJ) $(COBJ)
	$(LD) $(LDFLAGS) $^ $(LDLIBS) -o $@

%_cpp.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%_c.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rfv $(COBJ) $(CXXOBJ) $(EXE) $(LIB)
