/*	File: fftUtils.h 
	
	Description:
		Common public utility functions for FFT tests.
	
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
 * fftUtils.h - common public utility functions for FFT tests.
 *
 * Created 9/23/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_FFT_UTILS_H_
#define _FFT_UTILS_H_

#include <libMatrixFFT/fftPrecision.h>
#include <libMatrixFFT/vdspUtils.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#pragma mark --- Test setup and logging ---

/* 
 * Get current OS-measured "user time", or absolute (wall) time, in seconds, 
 * as a double.
 */
double fftGetTime(bool wallTime);		/* true --> CFAbsoluteTimeGetCurrent() 
										 * false --> user time via getrusage() */

/* print standard test banner */
void fftPrintTestBanner(
	const char *fftType,		/* e.g. "One-dimension real" */
	const char *impl,			/* e.g. "Accelerate", "FFTW" */
	bool doublePrec,
	const char *signalType,		/* e.g. "random", "chirp" */
	const char *options,		/* optional */
	unsigned loops,
	unsigned numThreads);

/* append a string suitable for 'options' arg to fftPrintTestBanner() */
extern void appendOptStr(
	char *optStr,
	const char *newOpt);

/* 
 * Obtain a compact string representation, in K or M as appropriate,
 * of specified input size_t. Caller must free() the result.
 */
extern char *fftStringRep(
	size_t inSize);

/* 
 * Obtain compact string representation as power of 2 if possible, else
 * as an unsigned long long. Caller must free() the result. 
 */
extern char *fftStringRepPow2(
	size_t inSize);

/* 
 * Parse an input string, possibly with a 'K', 'k', 'M', or 'm' suffix,
 * or as 2^xxx, or as plain decimal, into a size_t.
 */
extern size_t fftParseStringRep(
	const char *inStr);

/* 
 * Round up size to next power of 2.
 *
 * Returns 1, 2, 4, 8....
 * Returns 0, 1, 2, 3... in *log2NumSamples
 */
size_t fftRoundNumSamples(
	size_t numSamples,
	unsigned *log2NumSamples);	/* RETURNED */

#define COND_FREE(p)	if(p != NULL) { free(p); p = NULL; }

#define	max(a, b)	((a) < (b) ? (b) : (a))

#pragma mark --- Misc. functions ---

/* 
 * Determine if a number is a power of two; if so, return true and  
 * (optionally) the exponent in exp. Else return false and return the
 * largest power of 2 less than num.
 */
extern bool fftIsPowerOfTwo(
	size_t num,
	unsigned *exp);		/* optionally RETURNED */

/* 
 * Determine the number of active virtual CPU cores. 
 */
extern unsigned numCpuCores();

/*
 * Determine L2+L3+... cache size per CPU core, i.e. total cache 
 * size - other than L1 - per core.
 */
extern uint64_t l2CachePerCore();

/*
 * Flush all cache associated with specified memory range.
 */
extern void fftFlushCache(
    void *addr,
    size_t numBytes);

/*
 * Determine number of iterations to measure, given signal size.
 */
extern unsigned fftIterationsForSize(
    size_t numComplex);
    
/*
 * Max vDSP sizes on a 32 GB machine with 64 bit binaries. 
 */
#if		FFT_DOUBLE_PREC

#define VDSP_MAX_1D_COMPLEX		((size_t)1 << 28)
#define VDSP_MAX_2D_COMPLEX		((size_t)1 << 28)
#define VDSP_MAX_1D_REAL		((size_t)1 << 28)
#define VDSP_MAX_2D_REAL		((size_t)1 << 28)

#else	/* single precision */

#define VDSP_MAX_1D_COMPLEX		((size_t)1 << 28)
#define VDSP_MAX_2D_COMPLEX		((size_t)1 << 28)
#define VDSP_MAX_1D_REAL		((size_t)1 << 28)
#define VDSP_MAX_2D_REAL		((size_t)1 << 28)


#endif	/* FFT_DOUBLE_PREC */

#ifdef __cplusplus
}
#endif

#endif	/* _FFT_UTILS_H_ */

