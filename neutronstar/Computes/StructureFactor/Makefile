CSRC=$(wildcard *.c)
COBJ=$(patsubst %.c,%_c.o,$(CSRC))
LIB=libssf.so

CC=gcc
LD=gcc

LDLIBS=-lgomp
CFLAGS=-g -O3 -fPIC -std=gnu99 -fopenmp -Wall

default: library

library: $(LIB)

$(LIB): $(COBJ)
	$(LD) $^ $(LDLIBS) $(LDFLAGS) -shared -o $@

%_c.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rfv $(COBJ) $(LIB)
