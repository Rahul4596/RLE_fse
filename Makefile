GPP=gcc
CCFLAGS=-msse -Wall
CCFLAGS+=-g
CCFLAGS+=-O3
LDFLAGS=-lm 

OBJECTS=src/bfio.o \
	src/image_io.o \
	src/alloc.o 

all: compress

compress: $(OBJECTS) src/RLE_copy_2.c Makefile
	$(GPP) $(CCFLAGS) $(OBJECTS) src/RLE_FSE_2.c -o RLE $(LDFLAGS)

%.o : %.c
	$(GPP) $(CCFLAGS) -o $@ -c $<

clean:
	rm -f compress
	rm -f src/*.o src/*~
	rm -f *.o *~
