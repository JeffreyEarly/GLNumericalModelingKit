/*	File: fftPriv.cpp 
	
	Description:
		Common private utility functions for MatrixFFT library.
	
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
 * fftPriv.cpp - common private utility functions for MatrixFFT library.
 *			  
 * Note: the precision of the floats used in this module are determined at 
 * compile time via the MatrixFFTConfig.h module. 
 *
 * Created 9/25/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "fftPriv.h"
#include "fftDebug.h"
#include <stdio.h>


/* dump a submatrix to stdout */

#if		FFT_SPLIT_COMPLEX
void fftDumpSub(	
	const char			*title,
	const FFTFloat		*buf,
	size_t				rowSize)	/* in FFTFloats */
{
	if(title) {
		printf("%s:\n", title);
	}
	
	for(size_t row=0; row<FFT_FLOATS_PER_SUBMATRIX; row++) {
		for(size_t col=0; col<FFT_FLOATS_PER_SUBMATRIX; col++) {
			size_t off = (row * rowSize) + col;
			printf("%.2f, ", buf[off]);
		}
		printf("\n");
	}
}

#else	/* FFT_SPLIT_COMPLEX */

void fftDumpSub(	
	const char			*title,
	const FFTComplex	*buf,
	size_t				rowSize)	/* in FFTFloats */
{
	if(title) {
		printf("%s:\n", title);
	}
	
	for(size_t row=0; row<FFT_COMPLEX_PER_SUBMATRIX; row++) {
		for(size_t col=0; col<FFT_COMPLEX_PER_SUBMATRIX; col++) {
			size_t off = (row * rowSize) + col;
			printf("{%.2f,%.2f}, ", buf[off].real, buf[off].imag);
		}
		printf("\n");
	}
}

#endif	/* FFT_SPLIT_COMPLEX */

/* Enforce alignment of FFTComplexs */

#ifdef	DEBUG
#if		FFT_SPLIT_COMPLEX

void fftAssertBufAlign(
	FFTComplex	*buf,
	size_t		alignSize)
{
	RFASSERT(FFT_IS_ALIGNED(buf->real, alignSize));
	RFASSERT(FFT_IS_ALIGNED(buf->imag, alignSize));
}

#else	/* !FFT_SPLIT_COMPLEX */

void fftAssertBufAlign(
	FFTComplex	*buf,
	size_t		alignSize)
{
	RFASSERT(FFT_IS_ALIGNED(buf, alignSize));
}

#endif	/* FFT_SPLIT_COMPLEX */
#endif	/* DEBUG */


/* 
 * Alloc or realloc mfftPlan->auxBuf. 
 */
MFFTReturn fftAllocAuxBuf(
	MatrixFFTPlan		mfftPlan,
	size_t				newSize)		/* in in FFTComplex elements */
{
	if(mfftPlan->auxBufSize >= newSize) {
		return MR_Success;
	}
	if(mfftPlan->auxBufSize) {
		fftFreeComplexArrayAlign(mfftPlan->auxBuf, mfftPlan->auxBufFree);
		mfftPlan->auxBuf = NULL;
		mfftPlan->auxBufFree = NULL;
		mfftPlan->auxBufSize = 0;
	}
	mdprint("allocing auxBuf, newSize %llu\n", (unsigned long long)newSize);
	mfftPlan->auxBuf = fftAllocComplexArrayAlign(newSize, CACHE_LINE_SIZE, &mfftPlan->auxBufFree);
	if(mfftPlan->auxBuf == NULL) {
		printf("***malloc error in fftAllocAuxBuf\n");
		return MR_Memory;
	}
	fftAssertBufAlign(mfftPlan->auxBuf, CACHE_LINE_SIZE);
	mfftPlan->auxBufSize = newSize;
	return MR_Success;
}

#pragma mark --- Threaded FFTvScale() ---

/***
 *** This is not currently used in the library; it was written to optimize normalization, 
 *** which is done much more efficiently in other ways now.
 ***/
 
/* Called out from thread module */
static MFFTReturn fftScaleThrCallback(TP_TaskUnion *u)
{
	TP_ComplexScale *tcs = &u->scale;
	
	tpThreadDebug("fftScaleThrCallback top: startElt %lu eltsToScale %lu\n", 
		(unsigned long)tcs->startElt, (unsigned long)tcs->eltsToScale);

	#if		FFT_SPLIT_COMPLEX
	
	FFTComplex buf = *tcs->buf;
	buf.real += tcs->startElt;
	buf.imag += tcs->startElt;
	fftScaleComplex(&buf, tcs->scaleFact, tcs->eltsToScale);
	
	#else	/* interleaved */
	
	FFTComplex *buf = tcs->buf + tcs->startElt;
	fftScaleComplex(buf, tcs->scaleFact, tcs->eltsToScale);
	
	#endif	/* FFT_SPLIT_COMPLEX*/
	
	return MR_Success;
}

#define MIN_SCALE_ELTS_PER_THREAD	(256 * 1024)

MFFTReturn fftScaleThr(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	FFTFloat			scaleFact,
	
	/* total size of buf */
	size_t				totalElts)
{
	tpThreadDebug("fftScaleThr top: totalElts %lu\n", (unsigned long)totalElts);

	unsigned numThreads = mfftPlan->threadPoolP->numThreads;
	if(numThreads <= 1) {
		/* No threading */
		fftScaleComplex(buf, scaleFact, totalElts);
		return MR_Success;
	}
	
	/* To ensure good alignment we're going to force numThreads to be a power of 2 */
	while(!fftIsPowerOfTwo(numThreads, NULL)) {
		numThreads--;
		RFASSERT(numThreads >= 2);
	}
	RFASSERT(numThreads >= 2);

	size_t eltsPerThread = totalElts / numThreads;
	while(eltsPerThread < MIN_SCALE_ELTS_PER_THREAD) {
		numThreads--;
		if(numThreads <= 1) {
			break;
		}
		eltsPerThread = totalElts / numThreads;
	}
	
	if(numThreads <= 1) {
		/* No threading */
		fftScaleComplex(buf, scaleFact, totalElts);
		return MR_Success;
	}
		
	/* set up a TP_ComplexScale for each thread */	
	for(unsigned dex=0; dex<numThreads; dex++) {
		TP_PerThread *pt     = &mfftPlan->threadPoolP->perThread[dex];
		TP_Task *task        = &pt->task;
		TP_ComplexScale *tcs = &task->u->scale;

		task->op         = TPO_App;
		task->threadFcn  = fftScaleThrCallback;
		tcs->buf         = buf;
		tcs->scaleFact   = scaleFact;
		tcs->totalElts   = totalElts;
		tcs->startElt    = dex * eltsPerThread;
		tcs->eltsToScale = eltsPerThread;
	}
	
	/* GO */
	MFFTReturn ourRtn = tpThreadDispatch(mfftPlan->threadPoolP, numThreads);
	if(ourRtn) {
		return ourRtn;
	}
	ourRtn = tpThreadFinish(mfftPlan->threadPoolP, numThreads);
	return ourRtn;

}
