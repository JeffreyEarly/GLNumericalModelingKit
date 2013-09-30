/*	File: complexBufUtils.cpp 
	
	Description:
		FFTComplex buffer routines.
	
	Copyright:
		Copyright (C) 2009 Apple Inc.  All rights reserved.
	
	Disclaimer:
		IMPORTANT:  This Apple software is supplied to you by Apple
		Computer, Inc. ("Apple") in consideration of your agreement to
		the following terms, and your use, installation, modification
		or redistribution of this Apple software constitutes acceptance
		of these terms.  If you do not agree with these terms, please
		do not use, install, modify or redistribute this Apple
		software.

		In consideration of your agreement to abide by the following
		terms, and subject to these terms, Apple grants you a personal,
		non-exclusive license, under Apple’s copyrights in this
		original Apple software (the "Apple Software"), to use,
		reproduce, modify and redistribute the Apple Software, with or
		without modifications, in source and/or binary forms; provided
		that if you redistribute the Apple Software in its entirety and
		without modifications, you must retain this notice and the
		following text and disclaimers in all such redistributions of
		the Apple Software.  Neither the name, trademarks, service
		marks or logos of Apple Computer, Inc. may be used to endorse
		or promote products derived from the Apple Software without
		specific prior written permission from Apple.  Except as
		expressly stated in this notice, no other rights or licenses,
		express or implied, are granted by Apple herein, including but
		not limited to any patent rights that may be infringed by your
		derivative works or by other works in which the Apple Software
		may be incorporated.

		The Apple Software is provided by Apple on an "AS IS" basis.
		APPLE MAKES NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING
		WITHOUT LIMITATION THE IMPLIED WARRANTIES OF NON-INFRINGEMENT,
		MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, REGARDING
		THE APPLE SOFTWARE OR ITS USE AND OPERATION ALONE OR IN
		COMBINATION WITH YOUR PRODUCTS.

		IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT,
		INCIDENTAL OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
		TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
		DATA, OR PROFITS; OR BUSINESS INTERRUPTION) ARISING IN ANY WAY
		OUT OF THE USE, REPRODUCTION, MODIFICATION AND/OR DISTRIBUTION
		OF THE APPLE SOFTWARE, HOWEVER CAUSED AND WHETHER UNDER THEORY
		OF CONTRACT, TORT (INCLUDING NEGLIGENCE), STRICT LIABILITY OR
		OTHERWISE, EVEN IF APPLE HAS BEEN ADVISED OF THE POSSIBILITY OF
		SUCH DAMAGE.
*/
/*
 * complexBufUtils.cpp - FFTComplex buffer routines.
 *			  
 * Created 01/07/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#include <libMatrixFFT/complexBufUtils.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libMatrixFFT/devRandom.h>
#include <libMatrixFFT/fftUtils.h>
#include "PolyComplex.h"

/*
 * Allocate an array of well-aligned memory. 
 * Free the returned freePtr free().
 */
void *fftAllocAlign(
    size_t                  allocSize,      // to be allocated
    void                    **freePtr)      // to be freed
{
    void *ptr = malloc(allocSize + FFT_MEM_ALIGNMENT);
    if(ptr == NULL) {
        return NULL;
    }
    *freePtr = ptr;
    return (void *)FFT_ALIGN(ptr, FFT_MEM_ALIGNMENT);
}

#pragma mark --- Signal initialization ---

#define RRAND_BUFSIZE	2000

void genRandComplex(
	FFTComplex *buf,
	size_t numSamples)		/* in each of {real,imag} */
{
	unsigned bufSize = RRAND_BUFSIZE;
	
	if(bufSize > (numSamples * 2)) {
		bufSize = numSamples * 2;
	}
	
	/* fill up our buffer */
	FFTFloat *inBuf = (FFTFloat *)malloc(bufSize * sizeof(FFTFloat));
	FFTFloat *inp = inBuf;
	for(size_t dex=0; dex<bufSize; dex++) {
		*inp++ = getRandomDouble();
	}

	/* copy to caller's buffer */
	inp = inBuf;
	FFTFloat *inEnd = inBuf + bufSize;
	PolyComplex pc(buf);
	
	for(size_t dex=0; dex<numSamples; dex++) {
		pc.real(*inp++);
		pc.imag(*inp++);
		++pc;
		if(inp == inEnd) {
			inp = inBuf;
		}
	}
	free(inBuf);
}

void genIncrComplex(
	FFTComplex *buf,
	size_t numSamples)			/* in each of {real,imag} */
{
	size_t n = numSamples * 2;
	PolyComplex pc(buf);
	for(size_t dex=0; dex<n; ) {
		pc.real((FFTFloat)(dex++));
		pc.imag((FFTFloat)(dex++));
		++pc;
	}
}

/* generate random FFTFloat signal */
void genRandFloat(
	FFTFloat *buf,
	size_t numSamples)
{
	unsigned bufSize = RRAND_BUFSIZE;
	
	if(bufSize > numSamples) {
		bufSize = numSamples;
	}
	
	/* fill up our buffer */
	FFTFloat *inBuf = (FFTFloat *)malloc(bufSize * sizeof(FFTFloat));
	FFTFloat *inp = inBuf;
	for(size_t dex=0; dex<bufSize; dex++) {
		*inp++ = getRandomDouble();
	}

	/* copy to caller's buffer */
	FFTFloat *outp = buf; 
	inp = inBuf;
	FFTFloat *inEnd = inBuf + bufSize;
	
	for(size_t dex=0; dex<numSamples; dex++) {
		*outp++ = *inp++;
		if(inp == inEnd) {
			inp = inBuf;
		}
	}
	free(inBuf);
}

/* generate incrementing FFTFloat signal */
void genIncrFloat(
	FFTFloat *buf,
	size_t numSamples)
{
	for(size_t dex=0; dex<numSamples; ) {
		*buf++ = dex++;
	}
}

/* copy from float array (i.e. vImage format) to FFTComplex */
void fftCopyFloatToSplit(
	const float *src,
	FFTComplex *dst,
	size_t numSamples)			/* total real-signal samples */
{
	PolyComplex pc(dst);
	
	for(size_t dex=0; dex<numSamples; dex+=2) {
		pc.real(*src++);
		if(dex == (numSamples - 1)) {
			/* odd number of samples, done */
			break;
		}
		pc.imag(*src++);
		++pc;
	}
}

/* 
 * Return random double such that  
 *      minVal < rtn < maxVal
 */
#define RDR_RES 0xffff     /* should NOT be pwr of two to ensure non-integral values */

static FFTFloat randDoubleRange(
    FFTFloat minVal,
    FFTFloat maxVal)
{
    /* 
     * First a random double 0 < rd < 1.0 
     */
    unsigned randInt = genRandInt(1, RDR_RES-1);
    FFTFloat rd = (FFTFloat)randInt / (FFTFloat)(RDR_RES);
    
    /* 
     * Now get the value that far into the range of minVal to maxVal
     */
    return minVal + ((maxVal - minVal) * rd);
}

/*
 * Generate random tau array. 
 */
#define TUA_ALWAYS_NONNEG   0

void fftGenTau(
    FFTFloat *tau,
    size_t N)
{
    FFTFloat dN = (double)(N<<1);
    FFTFloat mdN = TUA_ALWAYS_NONNEG ? 0.0 : -dN;
    for(size_t j=0; j<N; j++) {
        tau[j] = randDoubleRange(mdN, dN);
    }
}

#pragma mark --- Buffer comparison and analysis ---

/* 
 * Analyze amplitudes of two complex signals. Return max delta and RMSE between
 * the two.  
 */
void fftCompareAmplitudes(
	const FFTComplex *buf1,
	const FFTComplex *buf2,
	size_t numSamples,			/* In each of real & imag */
	double *maxDelta,			/* RETURNED */
	double *rmse)				/* RETURNED */
{
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	PolyComplex pc1((FFTComplex *)buf1);
	PolyComplex pc2((FFTComplex *)buf2);
	
	for(size_t dex=0; dex<numSamples; dex++) {
		double d = pc1.real() - pc2.real();
		double deltaSquare = d * d;
		d = pc1.imag() - pc2.imag();
		deltaSquare += (d * d);
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		++pc1;
		++pc2;
	}
	*maxDelta = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);
}


/* 
 * Analyze amplitudes of two real signals. Return max delta and RMSE between
 * the two.  
 */
void fftCompareRealAmplitudes(
	const FFTComplex *buf1,
	const FFTComplex *buf2,
	size_t numSamples,			/* total real samples */
	double *maxDelta,			/* RETURNED */
	double *rmse)				/* RETURNED */
{
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	PolyComplex pc1((FFTComplex *)buf1);
	PolyComplex pc2((FFTComplex *)buf2);
	
	for(size_t dex=0; dex<numSamples/2; dex++) {
		double d = pc1.real() - pc2.real();
		double deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		
		d = pc1.imag() - pc2.imag();
		deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		++pc1;
		++pc2;
	}
	*maxDelta = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);
}

/* 
 * Find min and max values of amplitudes of all the elements of an FFT output.
 */ 
void fftAnalyzeAmplitudes(
	const FFTComplex *zBuf,
	size_t numSamples,
	double *minAmpl,		/* RETURNED */
	double *maxAmpl)		/* RETURNED */
{
	double minA = 100000000.0;
	double maxA = 0.0;
	PolyComplex pc((FFTComplex *)zBuf);
	
	for(size_t dex=0; dex<numSamples; dex++) {
		double re = pc.real();
		double im = pc.imag();
		++pc;
		double ampSq = (re * re) + (im * im);
		double amp = sqrt(ampSq);
		if(amp < minA) {
			minA = amp;
		}
		if(amp > maxA) {
			maxA = amp;
		}
	}
	*minAmpl = minA;
	*maxAmpl = maxA;
}

#pragma mark --- Buffer dump/display ---

#define MAX_FLT_TO_PRINT	16

void fftDump1DFloat(	
    const char			*title,
    const FFTFloat		*buf,
    size_t				numSamples)	/* max to print */
{
	if(title) {
		printf("%s:\n", title);
	}
	size_t toPrint = MAX_FLT_TO_PRINT * MAX_FLT_TO_PRINT;
	if(toPrint > numSamples) {
		toPrint = numSamples;
	}
	for(size_t off=0; off<toPrint; off++) {
		printf("%.4f, ", buf[off]);
		if(((off % MAX_FLT_TO_PRINT) == (MAX_FLT_TO_PRINT - 1)) && (off != toPrint)) {
			printf("\n");
		}
	}
	printf("\n");
}

void fftDumpBuf(	
	const char			*title,
	const FFTFloat		*buf,
	size_t				numRows,
	size_t				numCols)
{
	if(title) {
		printf("%s:\n", title);
	}
	
	const char *eolStr = "\n";
	const char *endStr = NULL;
	size_t rowsToPrint = numRows;
	if(rowsToPrint > MAX_FLT_TO_PRINT) {
		rowsToPrint = MAX_FLT_TO_PRINT;
		endStr = " ...\n";
	}
	size_t colsToPrint = numCols;
	if(colsToPrint > MAX_FLT_TO_PRINT) {
		colsToPrint = MAX_FLT_TO_PRINT;
		eolStr = " ...\n";
	}
	
	for(size_t row=0; row<rowsToPrint; row++) {
		for(size_t col=0; col<colsToPrint; col++) {
			size_t off = (row * numCols) + col;
			printf("%.2f, ", buf[off]);
		}
		printf("%s", eolStr);
	}
	if(endStr) {
		printf("%s", endStr);
	}
}

/* Dump an array of size_t's to stdout */
void fftDump1DSize(
    const char                  *title,
    const size_t                *buf,
    size_t                      numSamples)
{
	if(title) {
		printf("%s:\n", title);
	}
	size_t toPrint = MAX_FLT_TO_PRINT * MAX_FLT_TO_PRINT;
	if(toPrint > numSamples) {
		toPrint = numSamples;
	}
	for(size_t off=0; off<toPrint; off++) {
		printf("%6lu  ", (unsigned long)(buf[off]));
		if(((off % MAX_FLT_TO_PRINT) == (MAX_FLT_TO_PRINT - 1)) && (off != toPrint)) {
			printf("\n");
		}
	}
	printf("\n");
}
