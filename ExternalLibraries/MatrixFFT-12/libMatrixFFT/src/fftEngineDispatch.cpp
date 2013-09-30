/*	File: fftEngineDispatch.cpp  
	
	Description:
		Interface between MatrixFFT and FFTEngine.
	
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
 * fftEngineDispatch.cpp - Interface between MatrixFFT and FFTEngine. 
 *			  
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include "fftEngineDispatch.h"
#include "MatrixFFTPlan.h"
#include "ThreadPool.h"
#include "fftThreadOps.h"
#include "FFTEngine.h"
#include "fftPriv.h"
#include "fftDebug.h"
#include <math.h>

#define CE_LOG_SUMMARY	0
#if		CE_LOG_SUMMARY
#define csprintf(args...)	printf(args)
#else
#define csprintf(args...)
#endif

#pragma mark --- private routines ---

/*
 * ThreadPool callout. This calls a specific engine with a subset of work (relative to 
 * the task presented to our SPI) to perform. 
 * *ALL* work performed thru this module goes thru this single function. 
 */
static MFFTReturn fftEngineThreadCall(
	TP_TaskUnion *u)
{
	TP_Engine *ce = &u->engine;
	FFTEngine *ffte = ce->engine;
	FFTEngineOp eop = ce->op;
	
	RFASSERT(ffte != NULL);
	RFASSERT(ce->inBuf != NULL);
	RFASSERT(ce->outBuf != NULL);
	RFASSERT(ce->eltsToProcess != 0);
	ceprintf("%s %s start %llu eltsToProc %llu\n",
		ffte->mEngineName, fftEngineOps[eop], 
		(unsigned long long)ce->startElt, (unsigned long long)ce->eltsToProcess);
	
	MFFTReturn mrtn = MR_Internal;
	ffte->startTimer();
	switch(eop) {
		case FEO_FftRows:
			mrtn = ffte->fftRows(ce->forward, ce->inBuf, ce->outBuf, ce->normFact, ce->startElt, 
				ce->eltsToProcess);
			break;
		case FEO_TwistFft:
			mrtn = ffte->twistPlusFftRow(ce->inBuf, ce->outBuf, ce->startElt, ce->eltsToProcess);
			break;
		case FEO_FftInvTwist:
			mrtn = ffte->fftRowPlusInvTwist(ce->inBuf, ce->outBuf, ce->startElt, ce->eltsToProcess);
			break;
		case FEO_FftCols:
			mrtn = ffte->fftColumns(ce->forward, ce->inBuf, ce->outBuf, ce->normFact, ce->startElt, 
				ce->eltsToProcess);
			break;
		case FEO_FftRowsReal:
			mrtn = ffte->fftRowsReal(ce->forward, ce->inBuf, ce->outBuf, ce->normFact, ce->startElt, 
				ce->eltsToProcess);
			break;
		case FEO_RealTwist:
			mrtn = ffte->realTwist(ce->forward, ce->inBuf, ce->outBuf, ce->startElt, ce->eltsToProcess);
			break;
		default:
			RFASSERT(0);
			break;
	}
	ffte->stopTimer();
	return mrtn;
}	

#define LB_LOG_OPS		0
#if		LB_LOG_OPS
#define lbprintf(args...)	printf(args)
#else
#define lbprintf(args...)
#endif

#if		LB_LOG_OPS
static void fftEngineLogTimes(
	MatrixFFTPlan		mfftPlan,	
	FFTEngineOp			op,
	unsigned			numEngines)
{
	for(unsigned dex=0; dex<numEngines; dex++) {
		FFTEngine *engine = mfftPlan->fftEngines[dex];
		/* timers == 0 means this engine wasn't used */
		if(engine->startTime() != 0.0) {
			printf("%s: %-12s load %.2f elapsed %.2f ms\n", 
				engine->mEngineName, fftEngineOps[op],
				engine->mOpRatios[op],
				engine->elapsed() * 1000.0);
		}
	}
}
#else	/* LB_LOG_OPS */
#define fftEngineLogTimes(mfftPlan, op, numEngines)
#endif	/* LB_LOG_OPS */

/*
 * Adjust mOpRatios for the engines we just performed work with. 
 * Currently a nop in this library's configuration.
 */
static void fftEngineBalance(
	MatrixFFTPlan		mfftPlan,	
	FFTEngineOp			op,
	unsigned			numEngines)
{
}

/* 
 * Core function to divvy up one of the tasks expressed in our public SPI
 * between the various FFTEngines we know about. We fire off a worker thread for
 * each engine and wait for completion. 
 */
static MFFTReturn fftEngineDispatch(
	MatrixFFTPlan		mfftPlan,	
	FFTEngineOp			op,
	bool				forward,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact,
	size_t				startElt = 0)
{	
	#if		CE_LOG_SUMMARY
	double startTime = CFAbsoluteTimeGetCurrent();
	double endTime;
	#endif	/* CE_LOG_SUMMARY */
	
	/* 
	 * Ask engine 0 for the number of elements to process for this op. 
	 * When we have a heterogenous collection of FFTEngine types, this
	 * will be much more complicated. For now, they are all the same. 
	 */
	size_t totalElts = mfftPlan->fftEngines[0]->totalEltsForOp(op);
	RFASSERT(totalElts > startElt);
	totalElts -= startElt;
	
	/*
	 * First handle unaligned startElt case to allow subsequent ops to perform 
	 * at max performance.
	 */
	size_t opAlign = mfftPlan->fftEngines[0]->alignmentForOp(op);
	size_t startEltMod = startElt % opAlign;
	if((startEltMod != 0) && (totalElts > opAlign)) {
		size_t unalignedElts = opAlign - startEltMod;

		/* one engine, low performance, but small size */
		TP_TaskUnion tu;
		TP_Engine *ce	  = &tu.engine;
		ce->engine		  = mfftPlan->fftEngines[0];
		ce->mfftPlan	  = mfftPlan;
		ce->op			  = op;
		ce->forward		  = forward;
		ce->inBuf		  = inBuf;
		ce->outBuf		  = outBuf;
		ce->startElt	  = startElt;
		ce->eltsToProcess = unalignedElts;
		ce->normFact      = normFact;
		MFFTReturn mrtn   = fftEngineThreadCall(&tu);
		fftEngineLogTimes(mfftPlan, op, 1);
		if(mrtn) {
			return mrtn;
		}
		
		startElt  += unalignedElts;
		totalElts -= unalignedElts;
		if(totalElts == 0) {
			#if CE_LOG_SUMMARY
			endTime = CFAbsoluteTimeGetCurrent();
			csprintf("%s: %.2fs\n", fftEngineOps[op], endTime - startTime);
			#endif
			return MR_Success;
		}
	}
	
	unsigned numEngines = mfftPlan->numHostEngines;
	if(numEngines == 1) {
		/* Trivial case - no divvying up, no thread dispatch */
		TP_TaskUnion tu;
		TP_Engine *ce	  = &tu.engine;
		ce->engine		  = mfftPlan->fftEngines[0];
		ce->mfftPlan	  = mfftPlan;
		ce->op			  = op;
		ce->forward		  = forward;
		ce->inBuf		  = inBuf;
		ce->outBuf		  = outBuf;
		ce->startElt	  = startElt;
		ce->eltsToProcess = totalElts;
		ce->normFact      = normFact;
		MFFTReturn mrtn   = fftEngineThreadCall(&tu);
		fftEngineLogTimes(mfftPlan, op, 1);
		#if CE_LOG_SUMMARY
		endTime = CFAbsoluteTimeGetCurrent();
		csprintf("%s: %.2fs\n", fftEngineOps[op], endTime - startTime);
		#endif
		return mrtn;
	}
	
	RFASSERT(mfftPlan->threadPoolP->numThreads >= numEngines);
	size_t eltsToGo = totalElts;
	TP_PerThread *perThread = mfftPlan->threadPoolP->perThread;
	unsigned numThreads = 0;
	
	size_t currStart = startElt;
	
	for(unsigned dex=0; dex<numEngines; dex++) {
		FFTEngine *engine = mfftPlan->fftEngines[dex];
		RFASSERT(engine != NULL);
		size_t thisElts;
		if(dex == (numEngines - 1)) {
			/* last engine, take the remainder regardless */
			thisElts = eltsToGo;
		}
		else {
			/* Dole out according to ability per opRatios */
			thisElts = (size_t)roundf((float)totalElts * engine->mOpRatios[op]);
			
			/* This covers a corner case for very small ops... */
			if(thisElts == 0) {
				thisElts = 1;
			}
			/* align per engine requirements */
			opAlign = engine->alignmentForOp(op);
			size_t numStripes = (thisElts + opAlign - 1) / opAlign;
			thisElts = numStripes * opAlign;
			if(thisElts > eltsToGo) {
				thisElts = eltsToGo;
			}
		}
		RFASSERT(thisElts != 0);
		
		/* prepare a TP_Task for this engine */
		TP_Task *task		= &perThread->task;
		task->op			= TPO_App;
		task->threadFcn		= fftEngineThreadCall;
		TP_Engine *ce		= &task->u->engine;
		ce->engine			= engine;
		ce->mfftPlan		= mfftPlan;
		ce->op				= op;
		ce->forward			= forward;
		ce->inBuf			= inBuf;
		ce->outBuf			= outBuf;
		ce->startElt		= currStart;
		ce->eltsToProcess	= thisElts;
		ce->normFact        = normFact;
		
		currStart  += thisElts;
		eltsToGo   -= thisElts;

		numThreads++;
		perThread++;
		
		if(eltsToGo == 0) {
			/* could happen due to rounding to stripe boundary */
			break;
		}
	}

	/* GO */
	MFFTReturn mrtn = tpThreadDispatch(mfftPlan->threadPoolP, numThreads);
	if(mrtn) {
		return mrtn;
	}
	mrtn = tpThreadFinish(mfftPlan->threadPoolP, numThreads);
	
	if(mrtn == MR_Success) {
		fftEngineLogTimes(mfftPlan, op, numThreads);
		
		/* rebalance the load for next execution of this op */
		fftEngineBalance(mfftPlan, op, numThreads);
	}
	#if CE_LOG_SUMMARY
	endTime = CFAbsoluteTimeGetCurrent();
	csprintf("%s: %.2fs\n", fftEngineOps[op], endTime - startTime);
	#endif
	return mrtn;
}

#pragma mark --- Public SPI --- 

/* basic "FFT all rows" */
MFFTReturn feFftAllRows(
	MatrixFFTPlan		mfftPlan,	
	bool				forward,		
	size_t				startRow,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact)
{
	return fftEngineDispatch(mfftPlan, FEO_FftRows, forward, inBuf, outBuf, normFact, startRow);
}

/* 1-D forward twist, then FFT rows */
MFFTReturn feTwistPlusFftRow(
	MatrixFFTPlan		mfftPlan,	
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	return fftEngineDispatch(mfftPlan, FEO_TwistFft, true, inBuf, outBuf, 0.0);
}

/* (inverse) FFT, then inverse twist rows */
MFFTReturn feFftRowPlusInvTwist(
	MatrixFFTPlan		mfftPlan,	
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	return fftEngineDispatch(mfftPlan, FEO_FftInvTwist, false, inBuf, outBuf, 0.0);
}
	
/* FFT columns with optional startElt; input and output in row order */
MFFTReturn feFftColumns(
	MatrixFFTPlan		mfftPlan,	
	bool				forward,	
	size_t				startCol,
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact)
{
	return fftEngineDispatch(mfftPlan, FEO_FftCols, forward, inBuf, outBuf, normFact, startCol);
}

/* real-to-complex FFT all rows */
MFFTReturn feFftAllRowsReal(
	MatrixFFTPlan		mfftPlan,	
	bool				forward,	
	FFTComplex			*inBuf,
	FFTComplex			*outBuf,
	FFTFloat			normFact)
{
	return fftEngineDispatch(mfftPlan, FEO_FftRowsReal, forward, inBuf, outBuf, normFact);
}

/* 1-D real twist */
MFFTReturn fe1DRealTwist(
	MatrixFFTPlan		mfftPlan,	
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf)
{
	return fftEngineDispatch(mfftPlan, FEO_RealTwist, forward, inBuf, outBuf, 0.0);
}

