/*	File: fftwUtils.h
	
	Description:
		Common utility functions for FFTW tests.
	
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
 * fftwUtils.h - common utility functions for FFTW tests.
 *
 * Created 01/08/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_FFTW_UTILS_H_
#define _FFTW_UTILS_H_

#include <fftw3.h>
#include <libMatrixFFT/devRandom.h>
#include <math.h>
#include <stdlib.h>

/* default plan flag */
#if 0
/* This can take 12 hours or more to plan a single FFT...only use this
 * when you have a LOT of time.
 */
#define FFTW_PLAN_FLAG_DEF		FFTW_PATIENT
#define FFTW_PLAN_FLAG_DEF_STR	"PATIENT"
#else
#define FFTW_PLAN_FLAG_DEF		FFTW_MEASURE
#define FFTW_PLAN_FLAG_DEF_STR	"MEASURE"
#endif

/* 
 * Generate random real data (-1 <= x <= 1) signal into a caller-allocated 
 * 2-D array.
 * It takes a LONG time to generate a lot of random floats, so we just generate
 * a small(er) number and reuse them.
 */
#define RAND_BUFSIZE	2000

#define FFT_CONST_DATA	0

#if		FFT_CONST_DATA

template <class elt>
static void genRandSignal(
	elt *buf,
	size_t numSamples)
{
	for(size_t dex=0; dex<numSamples; dex++) {
		*buf++ = 1.0;
	}
}

#else	/* !FFT_CONST_DATA */

template <class elt>
static void genRandSignal(
	elt *buf,
	size_t numSamples)
{
	unsigned bufSize = RAND_BUFSIZE;
	
	if(bufSize > numSamples) {
		bufSize = numSamples;
	}
	
	/* fill up our buffer */
	elt *inBuf = (elt *)malloc(bufSize * sizeof(elt));
	elt *inp = inBuf;
	for(unsigned dex=0; dex<bufSize; dex++) {
		*inp++ = getRandomDouble();
	}
	
	/* copy to caller's buffer */
	elt *outp = buf;
	inp = inBuf;
	elt *inEnd = inBuf + bufSize;
	for(size_t dex=0; dex<numSamples; dex++) {
		*outp++ = *inp++;
		if(inp == inEnd) {
			inp = inBuf;
		}
	}
	free(inBuf);
}

#endif	/* FFT_CONST_DATA */

/* 
 * Find min and max values of amplitudes of all the elements of an FFT output.
 */ 
template <class FFTWComplex, class elt>
static void analyzeAmplitudes(
	FFTWComplex *zBuf,
	unsigned numSamples,
	double *minAmpl,		/* RETURNED */
	double *maxAmpl)		/* RETURNED */
{
	double minA = 100000000.0;
	double maxA = 0.0;
	
	for(unsigned dex=0; dex<numSamples; dex++) {
		double re = (*zBuf)[0];
		double im = (*zBuf)[1];
		zBuf++;
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

#define COND_FFTW_FREE(p)	if(p) { fftw_free(p); }

#ifdef	__cplusplus
extern "C" {
#endif

/* 
 * Obtain and import wisdom
 * Returns nonzero on error (which should generally not be a big deal; our wisdom files
 * are in /tmp/, so they get nuked on each boot).
 */
extern int fftwGetWisdom(bool doublePrec);

/*
 * Save current wisdom to file
 */
extern int fftwSaveWisdom(bool doublePrec);

#ifdef	__cplusplus
}
#endif

#endif	/* _FFT_UTILS_H_ */

