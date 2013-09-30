/*	File: FFTEngine.cpp  
	
	Description:
		FFTEngine base class.
	
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
 * FFTEngine.cpp - FFTEngine base class
 *			  
 * Created 10/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include "FFTEngine.h"
#include "fftPriv.h"
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <stdlib.h>
#include "fftDebug.h"

const char *fftEngineOps[] = 
{
	"FftRows",
	"TwistFft",
	"FftInvTwist",
	"FftCols",
	"FftRowsReal",
	"RealTwist"
};


FFTEngine::FFTEngine(
	MatrixFFTPlan		_mfftPlan,
	size_t				_numRows,
	size_t				_numCols,
	bool				_real) :
		mOpRatios(), mEngineName(),
		mNumRows(_numRows), mNumCols(_numCols),
		mNumColsComplex(_real ? _numCols >> 1 : _numCols),
		mReal(_real), mMfftPlan(_mfftPlan),
		mNumStripeBufs(0), mStripeBufSize(0), mStripeBuf(NULL), mStripeBufFree(NULL),
		mStartTime(0.0), mEndTime(0.0)
{
	for(unsigned dex=0; dex<FEO_NumOps; dex++) {
		mOpRatios[dex] = 0.0;
	}
}

FFTEngine::~FFTEngine()
{
	for(unsigned dex=0; dex<mNumStripeBufs; dex++) {
		if(mStripeBufFree[dex] != NULL) {
			fftFreeComplexArrayAlign(mStripeBuf[dex], mStripeBufFree[dex]);
			mStripeBufFree[dex] = NULL;
		}
	}
	if(mStripeBuf != NULL) {
		free(mStripeBuf);
		mStripeBuf = NULL;
	}
	if(mStripeBufFree != NULL) {
		free(mStripeBufFree);
		mStripeBufFree = NULL;
	}
}

MFFTReturn FFTEngine::allocStripeBufs(
	unsigned	numStripeBufs, 
	size_t		_stripeBufSize,
	size_t		alignSize)
{
	RFASSERT(mStripeBuf == NULL);
	mStripeBuf = (FFTComplex **)malloc(numStripeBufs * sizeof(FFTComplex *));
	if(mStripeBuf == NULL) {
		return MR_Memory;
	}
	mStripeBufFree = (FFTComplex **)malloc(numStripeBufs * sizeof(FFTComplex *));
	if(mStripeBufFree == NULL) {
		return MR_Memory;
	}
	memset(mStripeBuf, 0, numStripeBufs * sizeof(FFTComplex *));
	memset(mStripeBufFree, 0, numStripeBufs * sizeof(FFTComplex *));
	
	for(unsigned dex=0; dex<numStripeBufs; dex++) {
		mStripeBuf[dex] = fftAllocComplexArrayAlign(_stripeBufSize, alignSize, &mStripeBufFree[dex]);
		if(mStripeBuf[dex] == NULL) {
			return MR_Memory;
		}
	}
	
	mNumStripeBufs = numStripeBufs;
	mStripeBufSize = _stripeBufSize;
	
	return MR_Success;
}

double FFTEngine::elapsed()
{
	RFASSERT((mStartTime != 0.0) && (mEndTime != 0.0));
	RFASSERT(mEndTime >= mStartTime);
	return mEndTime - mStartTime;
}

void FFTError::throwMe(MFFTReturn mrtn)	
{ 
	throw FFTError(mrtn); 
}

