#
# Makefile for build Fortran/MatrixFFT examples.
#
CC = g++
FC = gfortran
LD = gfortran

LDLIBS  = -lMatrixFFT -lm -lstdc++-static

#
# Currently gfortran does not build a 64-bit binary on Mac OS X; by default,
# as of Snow Leopard (Mac OS 10.6), 64-bit build is the default for C code
# of the current host is 64-bit capable. Thus (for now) we have to explcitly
# force a 32-bit build for compilation of the wrapper; you also have to build
# the 32-bit (i386) version of MatrixFFT to work with this.
#
ARCH_CFLAGS = -m32

all:	callFFTComplex callFFTComplex2D callFFTReal callFFTReal2D

callFFTComplex: fortranFFTWrapper.o callFFTComplex.o
	$(LD) -o callFFTComplex callFFTComplex.o fortranFFTWrapper.o $(LDLIBS) -framework Accelerate -framework CoreFoundation

callFFTComplex2D: fortranFFTWrapper.o callFFTComplex2D.o
	$(LD) -o callFFTComplex2D callFFTComplex2D.o fortranFFTWrapper.o $(LDLIBS) -framework Accelerate -framework CoreFoundation

callFFTReal: fortranFFTWrapper.o callFFTReal.o
	$(LD) -o callFFTReal callFFTReal.o fortranFFTWrapper.o $(LDLIBS) -framework Accelerate -framework CoreFoundation

callFFTReal2D: fortranFFTWrapper.o callFFTReal2D.o
	$(LD) -o callFFTReal2D callFFTReal2D.o fortranFFTWrapper.o $(LDLIBS) -framework Accelerate -framework CoreFoundation

%.o: %.f90
	$(FC) $< -c

fortranFFTWrapper.o: fortranFFTWrapper.cpp
	$(CC) fortranFFTWrapper.cpp  -g $(ARCH_CFLAGS) -c 

clean:
	rm -f callFFTReal callFFTComplex callFFTComplex2D callFFTReal2D *.o
