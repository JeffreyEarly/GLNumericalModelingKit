#
# Makefile for building MatrixFFT examples.
#
CC = gcc

LIBS=    -lMatrixFFT -lm -lstdc++
FRAMES=  -framework Accelerate -framework CoreFoundation
LIBPATH= 
LDFLAGS= $(LIBS) $(LIBPATH) $(FRAMES)
CFLAGS= -g -Wno-four-char-constants -Wall -Werror

EXECUTABLES= complex1Dfft complex2Dfft real1Dfft real2Dfft

all:	$(EXECUTABLES)

complex1Dfft: complex1Dfft.o
	$(CC) -o $@ $@.o $(LDFLAGS) -framework Accelerate -framework CoreFoundation

complex2Dfft: complex2Dfft.o
	$(CC) -o $@ $@.o $(LDFLAGS) -framework Accelerate -framework CoreFoundation

real1DFft: real1Dfft.o
	$(CC) -o $@ $@.o $(LDFLAGS) -framework Accelerate -framework CoreFoundation

real2DFft: real2Dfft.o
	$(CC) -o $@ $@.o $(LDFLAGS) -framework Accelerate -framework CoreFoundation

%.o: %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -f $(EXECUTABLES) *.o
