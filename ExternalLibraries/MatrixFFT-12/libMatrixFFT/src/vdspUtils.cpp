/*	File: vdspUtils.cpp 
	
	Description:
		Common vDSP utility functions for FFT tests.
	
	Copyright:
		Copyright (C) 2008 Apple Inc.  All rights reserved.
	
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
 * vdspUtils.cpp - common vDSP utility functions for FFT tests.
 *
 * Created 7/17/2007. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <libMatrixFFT/vdspUtils.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/devRandom.h>
#include "fftPriv.h"
#include "fftIntel.h"
#include "fftDebug.h"

#pragma mark --- Buffer allocation ---

/*
 * Allocate an vDSPComplex with specified number of complex elements.
 * Returns nonzero on malloc failure. Caller must free the returned
 * vDSPComplex pointers.
 */
int fftAllocDSPComplex(
	vDSPComplex				*buf,
	size_t					numComplex)
{
	size_t bufSize = numComplex * sizeof(FFTFloat);
	buf->realp = (FFTFloat *)malloc(bufSize);
	buf->imagp = (FFTFloat *)malloc(bufSize);
	if((buf->realp == NULL) || (buf->imagp == NULL)) {
		COND_FREE(buf->realp);
		COND_FREE(buf->imagp);
		printf("***fftAllocSplitComplex: malloc error\n");
		return -1;
	}
	return 0;
}

/*
 * Allocate an vDSPComplex with specified number of complex elements,
 * aligned to specified alignSize, which must be a power of 2. 
 * Returns nonzero on malloc failure. 
 * Caller must free the contents of the returned freeBuf. 
 */
int fftAllocDSPComplexAlign(
	vDSPComplex				*buf,
	size_t					numComplex,
	unsigned				alignSize,
	vDSPComplex				*freeBuf)
{
	size_t bufSize = (numComplex * sizeof(FFTFloat)) + alignSize;
	freeBuf->realp = (FFTFloat *)malloc(bufSize);
	freeBuf->imagp = (FFTFloat *)malloc(bufSize);
	if((freeBuf->realp == NULL) || (freeBuf->imagp == NULL)) {
		COND_FREE(freeBuf->realp);
		COND_FREE(freeBuf->imagp);
		printf("***fftAllocDSPComplexAlign: malloc error\n");
		return -1;
	}
	buf->realp = (FFTFloat *)FFT_ALIGN(freeBuf->realp, alignSize);
	buf->imagp = (FFTFloat *)FFT_ALIGN(freeBuf->imagp, alignSize);
	#ifdef	DEBUG
	if((freeBuf->realp != buf->realp) || (freeBuf->imagp != buf->imagp)) {
		printf("...caught unaligned DSP pointer\n");
	}
	#endif
	return 0;
}


void fftFreeDSPComplex(
	vDSPComplex				*buf)
{
	COND_FREE(buf->realp);
	COND_FREE(buf->imagp);
}

/* copy from one vDSPComplex to another */
void fftCopyDSPComplex(
	const vDSPComplex		*src,
	vDSPComplex				*dst,
	size_t					numSamples)
{
	size_t bufSize = numSamples * sizeof(FFTFloat);
	memmove(dst->realp, src->realp, bufSize);
	memmove(dst->imagp, src->imagp, bufSize);
}

/* Multiply each element in an vDSPComplex by specified scale factor */
void fftScaleDSPComplex(
	vDSPComplex *buf,
	FFTFloat scaleFact,
	size_t numSamples)		/* in each of realp, imag */
{
	FFTvScale(buf->realp, scaleFact, numSamples);
	FFTvScale(buf->imagp, scaleFact, numSamples);
}

/* copy from float array (i.e. vImage format) to vDSPComplex */
void fftCopyFloatToDSPComplex(
	const float *src,
	vDSPComplex *dst,
	size_t numSamples)			/* total real-signal samples */
{
	FFTFloat *dstReal = dst->realp;
	FFTFloat *dstImag = dst->imagp;
	
	for(size_t dex=0; dex<numSamples; dex+=2) {
		*dstReal++ = *src++;
		if(dex == (numSamples - 1)) {
			/* odd number of samples, done */
			break;
		}
		*dstImag++ = *src++;
	}
}

#pragma mark --- Buffer initialization ---

/* 
 * Generate random real data (-1 <= x <= 1) signal into a caller-allocated 
 * DSPSplitComplex or DSPDoubleSplitComplex.
 * It takes a LONG time to generate a lot of random floats, so we just generate
 * a small(er) number and reuse them.
 */
#define RAND_BUFSIZE	2000

void genRandComplexDSP(
	vDSPComplex *zBuf,
	size_t numSamples)
{
	unsigned bufSize = RAND_BUFSIZE;
	
	if(bufSize > numSamples * 2) {
		bufSize = numSamples * 2;
	}
	
	/* fill up our buffer */
	FFTFloat *inBuf = (FFTFloat *)malloc(bufSize * sizeof(FFTFloat));
	FFTFloat *inp = inBuf;
	for(unsigned dex=0; dex<bufSize; dex++) {
		*inp++ = getRandomDouble();
	}

	/* copy to caller's buffer */
	FFTFloat *realp = zBuf->realp; 
	FFTFloat *imagp = zBuf->imagp;
	inp = inBuf;
	FFTFloat *inEnd = inBuf + bufSize;
	
	for(size_t dex=0; dex<numSamples; dex++) {
		*realp++ = *inp++;
		*imagp++ = *inp++;
		if(inp == inEnd) {
			inp = inBuf;
		}
	}
	free(inBuf);
}

void genIncrComplexDSP(
	vDSPComplex *buf,
	size_t numSamples)			/* in each of {real,imag} */
{
	FFTFloat *realp = buf->realp; 
	FFTFloat *imagp = buf->imagp;
	
	size_t n = numSamples * 2;
	for(size_t dex=0; dex<n; ) {
		*realp++ = dex++;
		*imagp++ = dex++;
	}
}

/* Generate 1-dimension "chirp" signal */
void fftGenChirp1DDSP(
	vDSPComplex *buf,
	size_t numSamples)			/* in each of {real,imag} */
{
	FFTFloat *realp = buf->realp; 
	FFTFloat *imagp = buf->imagp;
	double pi = M_PI;
	double oneOverN = 1.0 / (FFTFloat)numSamples;
	
	for(size_t j=0; j<numSamples; j++) {
		/* arg := pi * j^2 / N */
		double arg = j * j;
		arg *= pi;
		arg *= oneOverN;
		*realp++ = cos(arg);
		*imagp++ = sin(arg);
	}
}

/* Generate 2-dimension "chirp" signal */
void fftGenChirp2DDSP(
	vDSPComplex *buf,
	size_t width,
	size_t height)
{
	FFTFloat *realp = buf->realp; 
	FFTFloat *imagp = buf->imagp;
	double pi = M_PI;
	
	for(size_t k=0; k<height; k++) {
		for(size_t j=0; j<width; j++) {
			/* arg := pi * (j^2/W + k^2/H) */
			double arg = ((double)(j * j)) / (double)width;
			arg += ((double)(k * k)) / (double)height;
			arg *= pi;
			*realp++ = cos(arg);
			*imagp++ = sin(arg);
		}
	}
}

/* flush contents of an vDSPComplex from cache */
void fftFlushComplexDSP(
	vDSPComplex *buf,
	size_t numComplex)
{
    size_t numBytes = numComplex * sizeof(FFTFloat);
    fftFlushCache(buf->realp, numBytes);
    fftFlushCache(buf->imagp, numBytes);
    
}

#pragma mark --- Buffer comparison and analysis ---

/* 
 * Analyze amplitudes of two complex signals. Return max delta and RMSE between
 * the two.  
 */
void fftCompareAmplitudesDSP(
	const vDSPComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* In each of real & imag */
	double *maxDelta,			/* RETURNED */
	double *rmse)				/* RETURNED */
{
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	const FFTFloat *real1 = buf1->realp;
	const FFTFloat *imag1 = buf1->imagp;
	const FFTFloat *real2 = buf2->realp;
	const FFTFloat *imag2 = buf2->imagp;
	
	for(size_t dex=0; dex<numSamples; dex++) {
		double d = *real1++ - *real2++;
		double deltaSquare = d * d;
		d = *imag2++ - *imag1++;
		deltaSquare += (d * d);
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
	}
	*maxDelta = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);
}

/* 
 * Analyze amplitudes of two real signals. Return max delta and RMSE between
 * the two.  
 */
void fftCompareRealAmplitudesDSP(
	const vDSPComplex *buf1,
	const vDSPComplex *buf2,
	size_t numSamples,			/* total real samples */
	double *maxDelta,			/* RETURNED */
	double *rmse)				/* RETURNED */
{
	double totalErr = 0.0;		/* running total of delta squared */
	double currMax = 0.0;		/* squared! */
	
	const FFTFloat *real1 = buf1->realp;
	const FFTFloat *imag1 = buf1->imagp;
	const FFTFloat *real2 = buf2->realp;
	const FFTFloat *imag2 = buf2->imagp;
	
	for(size_t dex=0; dex<numSamples/2; dex++) {
		double d = *real1++ - *real2++;
		double deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
		
		d = *imag1++ - *imag2++;
		deltaSquare = d * d;
		if(deltaSquare > currMax) {
			currMax = deltaSquare;
		}
		totalErr += deltaSquare;
	}
	*maxDelta = sqrt(currMax);
	*rmse = sqrt(totalErr / numSamples);
}


/* 
 * Find min and max values of amplitudes of all the elements of an FFT output.
 */ 
extern void fftAnalyzeAmplitudesDSP(
	const vDSPComplex *zBuf,
	size_t numSamples,
	double *minAmpl,		/* RETURNED */
	double *maxAmpl)		/* RETURNED */
{
	double minA = 100000000.0;
	double maxA = 0.0;
	const FFTFloat *rep = zBuf->realp;
	const FFTFloat *imp = zBuf->imagp;
	
	for(size_t dex=0; dex<numSamples; dex++) {
		double re = *rep++;
		double im = *imp++;
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

#define MAX_TO_PRINT		8
#define MAX_ROWS_TO_PRINT	16

/* 
 * Dump up to 8 rows and columns of a square vDSPComplex data to stdout.
 */
void fftDumpDSPMatrix(	
	const char				*title,
	const vDSPComplex		*buf,
	unsigned				n)				/* total size = 2 ^ 2n */
{
	if(title) {
		printf("%s:\n", title);
	}
	
	size_t x = (size_t)1 << n;
	size_t toPrint = x;
	const char *eolStr = "\n";
	if(toPrint > MAX_TO_PRINT) {
		toPrint = MAX_TO_PRINT;
		eolStr = " ...\n";
	}
	for(size_t row=0; row<toPrint; row++) {
		for(size_t col=0; col<toPrint; col++) {
			size_t off = (row * x) + col;
			printf("{%.2f,%.2f}, ", buf->realp[off], buf->imagp[off]);
		}
		printf("%s", eolStr);
	}
	if(x > MAX_TO_PRINT) {
		printf("%s", eolStr);
	}
}

/*
 * Dump up to 8 rows and columns of vDSPComplex data to stdout.
 */
void fftDumpDSPMatrixRect(	
	const char				*title,
	const vDSPComplex		*buf,
	size_t					numRows,
	size_t					numCols)
{
	if(title) {
		printf("%s:\n", title);
	}
	
	const char *eolStr = "\n";
	const char *endStr = NULL;
	size_t rowsToPrint = numRows;
	if(rowsToPrint > MAX_ROWS_TO_PRINT) {
		rowsToPrint = MAX_ROWS_TO_PRINT;
		endStr = " ...\n";
	}
	size_t colsToPrint = numCols;
	if(colsToPrint > MAX_TO_PRINT) {
		colsToPrint = MAX_TO_PRINT;
		eolStr = " ...\n";
	}
	
	for(size_t row=0; row<rowsToPrint; row++) {
		for(size_t col=0; col<colsToPrint; col++) {
			size_t off = (row * numCols) + col;
			printf("{%.2f,%.2f}, ", buf->realp[off], buf->imagp[off]);
		}
		printf("%s", eolStr);
	}
	if(endStr) {
		printf("%s", endStr);
	}
}

/* 
 * Dump 1-dimension DSPSplit*Complex data to stdout.
 */
template <class SplitComplex>
static void fftDump1DDSPComplexT(	
	const char				*title,
	const SplitComplex		*buf,
	size_t					numSamples)	/* max to print */
{
	if(title) {
		printf("%s:\n", title);
	}
	size_t toPrint = MAX_TO_PRINT * MAX_TO_PRINT;
	if(toPrint > numSamples) {
		toPrint = numSamples;
	}
	for(size_t off=0; off<toPrint; off++) {
		printf("{%.2f,%.2f}, ", buf->realp[off], buf->imagp[off]);
		if(((off % MAX_TO_PRINT) == (MAX_TO_PRINT - 1)) && (off != (toPrint-1))) {
			printf("\n");
		}
	}
    if(toPrint != numSamples) {
		printf("\n...");
	}
	else {
        printf("\n");	
	}
}

/* 
 * Dump 2-dimension DSPSplit*Complex data to stdout.
 */

template <class SplitComplex>
static void fftDump2DDSPComplexT(	
	const char				*title,
	const SplitComplex		*buf,
	size_t					numRows,	/* max to print */
	size_t					numCols)
{
	if(title) {
		printf("%s:\n", title);
	}
	size_t rowsToPrint = MAX_TO_PRINT;
	if(rowsToPrint > numRows) {
		rowsToPrint = numRows;
	}
	size_t colsToPrint = MAX_TO_PRINT;
	if(colsToPrint > numCols) {
		colsToPrint = numCols;
	}
	for(size_t row=0; row<rowsToPrint; row++) {
		for(size_t col=0; col<colsToPrint; col++) {
			size_t off = (row * numCols) + col;
			printf("{%.2f,%.2f}, ", buf->realp[off], buf->imagp[off]);
		}
		printf("\n");
	}
	printf("\n");
}


/* 
 * Dump 1- or 2-dimension DSPSplit*Complex data to stdout.
 */
void fftDump1DDSPComplexFloat(	
	const char					*title,
	const DSPSplitComplex		*buf,
	size_t						numSamples)		/* max to print */
{
	fftDump1DDSPComplexT<DSPSplitComplex>(title, buf, numSamples);
}

void fftDump1DDSPComplexDouble(	
	const char					*title,
	const DSPDoubleSplitComplex	*buf,
	size_t						numSamples)		/* max to print */
{
	fftDump1DDSPComplexT<DSPDoubleSplitComplex>(title, buf, numSamples);
}

void fftDump1DDSPComplex(	
	const char					*title,
	const vDSPComplex			*buf,
	size_t						numSamples)	/* max to print */
{
	#if		FFT_DOUBLE_PREC
	fftDump1DDSPComplexT<DSPDoubleSplitComplex>(title, buf, numSamples);
	#else
	fftDump1DDSPComplexT<DSPSplitComplex>(title, buf, numSamples);
	#endif
}

void fftDump2DDSPComplexFloat(	
	const char					*title,
	const DSPSplitComplex		*buf,
	size_t						numRows,	/* max to print */
	size_t						numCols)
{
	fftDump2DDSPComplexT<DSPSplitComplex>(title, buf, numRows, numCols);
}
	
void fftDump2DDSPComplexDouble(	
	const char					*title,
	const DSPDoubleSplitComplex	*buf,
	size_t						numRows,	/* max to print */
	size_t						numCols)
{
	fftDump2DDSPComplexT<DSPDoubleSplitComplex>(title, buf, numRows, numCols);
}

extern void fftDump2DDSPComplex(	
	const char					*title,
	const vDSPComplex			*buf,
	size_t						numRows,		/* max to print */
	size_t						numCols)
{
	#if		FFT_DOUBLE_PREC
	fftDump2DDSPComplexT<DSPDoubleSplitComplex>(title, buf, numRows, numCols);
	#else
	fftDump2DDSPComplexT<DSPSplitComplex>(title, buf, numRows, numCols);
	#endif
}

#pragma mark --- Convert between interleaved complex and vDSP split ---

#if		!FFT_SPLIT_COMPLEX

/* 
 * Routines to convert between interleaved complex format and vDSP's native
 * split format. 
 * First the brute force versions, for small signals and non-Intel.
 */
static void fftVDSPToIntSimple(
	const vDSPComplex *inBuf,
	FFTComplex *outBuf,
	size_t num)
{
	const FFTFloat *realp = inBuf->realp;
	const FFTFloat *imagp = inBuf->imagp;
	
	for(size_t dex=0; dex<num; dex++) {
		outBuf->real = *realp++;
		outBuf->imag = *imagp++;
		outBuf++;
	}
}

static void fftIntToVDSPSimple(
	const FFTComplex *inBuf,
	vDSPComplex *outBuf,
	size_t num)
{
	FFTFloat *realp = outBuf->realp;
	FFTFloat *imagp = outBuf->imagp;
	
	for(size_t dex=0; dex<num; dex++) {
		*realp++ = inBuf->real;
		*imagp++ = inBuf->imag;
		inBuf++;
	}
}

#if		FFT_INTEL

/* 
 * Intel vectorized versions.
 * USE_VDSP_CTOZ enables the use of vDSP_ctoz and its associates; this results
 * in a very a slight performance degradation when compared to our own vectorized 
 * version.
 */
#define USE_VDSP_CTOZ	0

#if		USE_VDSP_CTOZ

static void fftVDSPToIntVector(
	const vDSPComplex *inBuf,
	FFTComplex *outBuf,
	size_t num)
{
	vDSPComplexInt *dspInt = (vDSPComplexInt *)outBuf;
	FFT_vDSP_ztoc(inBuf, dspInt, num);
}

static void fftIntToVDSPVector(
	const FFTComplex *inBuf,
	vDSPComplex *outBuf,
	size_t num)
{
	vDSPComplexInt *dspInt = (vDSPComplexInt *)inBuf;
	FFT_vDSP_ctoz(dspInt, outBuf, num);
}

#else /* !USE_VDSP_CTOZ */

/* Custom hand-rolled vectorized versions */
static void fftVDSPToIntVector(
	const vDSPComplex *inBuf,
	FFTComplex *outBuf,
	size_t num)
{
	const FFTFloat *realp = inBuf->realp;
	const FFTFloat *imagp = inBuf->imagp;
	FFTVector rv;
	FFTVector iv;

	RFASSERT(FFT_IS_ALIGNED(realp, VECTOR_SIZE));
	RFASSERT(FFT_IS_ALIGNED(imagp, VECTOR_SIZE));
	
	for(size_t dex=0; dex<num; dex+=FFT_FLOATS_PER_VECTOR) {
		rv = FFTVectLoad(realp);
		iv = FFTVectLoad(imagp);
		fftStoreComplex(outBuf, iv, rv);
		realp  += FFT_FLOATS_PER_VECTOR;
		imagp  += FFT_FLOATS_PER_VECTOR;
		outBuf += FFT_FLOATS_PER_VECTOR;
	}
}

static void fftIntToVDSPVector(
	const FFTComplex *inBuf,
	vDSPComplex *outBuf,
	size_t num)
{
	FFTFloat *realp = outBuf->realp;
	FFTFloat *imagp = outBuf->imagp;
	FFTVector rv;
	FFTVector iv;

	RFASSERT(FFT_IS_ALIGNED(realp, VECTOR_SIZE));
	RFASSERT(FFT_IS_ALIGNED(imagp, VECTOR_SIZE));

	for(size_t dex=0; dex<num; dex+=FFT_FLOATS_PER_VECTOR) {
		fftLoadComplex(inBuf, iv, rv);
		FFTVectStore(realp, rv);
		FFTVectStore(imagp, iv);
		realp  += FFT_FLOATS_PER_VECTOR;
		imagp  += FFT_FLOATS_PER_VECTOR;
		inBuf += FFT_FLOATS_PER_VECTOR;
	}
}

#endif	/* USE_VDSP_CTOZ */
#endif	/* FFT_INTEL */

void fftVDSPToInt(
	const vDSPComplex *inBuf,
	FFTComplex *outBuf,
	size_t num)
{
	#if		FFT_INTEL
		if(FFT_IS_ALIGNED(outBuf, VECTOR_SIZE) &&
		   FFT_IS_ALIGNED(num, FFT_FLOATS_PER_VECTOR)) {
		   fftVDSPToIntVector(inBuf, outBuf, num);
		   return;
		}
	#endif	/* FFT_INTEL */
	
	/* else drop thru and take brute-force case */
   fftVDSPToIntSimple(inBuf, outBuf, num);
}

void fftIntToVDSP(
	const FFTComplex *inBuf,
	vDSPComplex *outBuf,
	size_t num)
{
	#if		FFT_INTEL
		if(FFT_IS_ALIGNED(inBuf, VECTOR_SIZE) &&
		   FFT_IS_ALIGNED(num, FFT_FLOATS_PER_VECTOR)) {
		   fftIntToVDSPVector(inBuf, outBuf, num);
		   return;
		}
	#endif	/* FFT_INTEL */
	
	/* else drop thru and take brute-force case */
   fftIntToVDSPSimple(inBuf, outBuf, num);
}

#endif	/* !FFT_SPLIT_COMPLEX */

