/*	File: FFTEngine.h  
	
	Description:
		Abstraction of low-level FFT implementations.
	
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
 * FFTEngine.h - abstraction of low-level FFT implementations.
 *			  
 * Created 10/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_FFT_ENGINE_H_
#define _FFT_ENGINE_H_

#include <libMatrixFFT/MatrixFFT.h>
#include "fftThreadOps.h"
#include <CoreFoundation/CoreFoundation.h>

/*
 * Abstract base class which defines the functionality (as an API) of 
 * various concrete subclasses.
 * In all of the functions to be implemented by the subclass, inBuf and 
 * outBuf point to [0, 0] in the entire array (i.e. NOT the first element 
 * to be processed).
 */
class FFTEngine
{
public:
	FFTEngine(
		MatrixFFTPlan		mfftPlan,
		size_t				numRows,
		size_t				numCols,
		bool				real);
	virtual ~FFTEngine();
	
	/***
	 *** These are the ops a subclass must implement
	 ***/
	 
	/* basic "complex FFT a contiguous set of rows" with optional normalization */
	virtual MFFTReturn fftRows(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		FFTFloat			normFact,
		size_t				startRow, 
		size_t				rowsToProcess) = 0;
		
	/* 1-D forward complex twist plus FFT a contiguous set of rows */
	virtual MFFTReturn twistPlusFftRow(
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		size_t				startRow, 
		size_t				rowsToProcess) = 0;
		
	/* inverse FFT then inverse complex twist a contiguous set of rows */
	virtual MFFTReturn fftRowPlusInvTwist(
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		size_t				startRow, 
		size_t				rowsToProcess) = 0;
		
	/* complex FFT a contiguous set of columns with optional normalization */
	virtual MFFTReturn fftColumns(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		FFTFloat			normFact,
		size_t				startCol, 
		size_t				colsToProcess) = 0;
	
	/* real-to-complex FFT a contiguous set of rows with optional normalization */
	virtual MFFTReturn fftRowsReal(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		FFTFloat			normFact,
		size_t				startRow, 
		size_t				rowsToProcess) = 0;
	
	/* 1-D real twist */
	virtual MFFTReturn realTwist(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		size_t				startRow, 
		size_t				rowsToProcess) = 0;
		
	/* determine number of rows or columns needed to process given op */
	virtual size_t totalEltsForOp(
		FFTEngineOp			op) = 0;
		
	/*
	 * Determine optimal alignment for specified op, in rows or columns
	 * as appropriate. When the engine dispatch code splits up a task 
	 * among multiple FFTEngines, the resulting sizes for all engines except 
	 * possibly the last one will be multiples of this alignment. 
	 */
	virtual size_t alignmentForOp(
		FFTEngineOp			op) = 0;
		
	/*** end of subclass-implemented methods ***/
	
	/* getters for read-only params */
	size_t	numRows()			{ return mNumRows;	}
	size_t	numCols()			{ return mNumCols;	}
	size_t	numColsComplex()	{ return mNumColsComplex; }
	bool	real()				{ return mReal;		}
	MatrixFFTPlan mfftPlan()	{ return mMfftPlan; }

	/*
	 * Array indicating the percentage, relative to a total task,
	 * that this object takes on. One value per op we can perform.
	 * E.g.: if a value for a given op is 0.5, then this object 
	 * performs 50% of the total work associated with that op. 
	 *
	 * This whole array is is public since we don't maintain it. It's 
	 * filled with zeroes during our constructor. 
	 */
	float					mOpRatios[FEO_NumOps];

	/* timing support */
	void					startTimer()	{ mStartTime = CFAbsoluteTimeGetCurrent(); }
	void					stopTimer()		{ mEndTime = CFAbsoluteTimeGetCurrent(); }
	void					resetTimers()	{ mStartTime = 0.0; mEndTime = 0.0; }
	double					elapsed();
	double					startTime()		{ return mStartTime; }
	
	/* for debugging */
	char					mEngineName[20];
	size_t					stripeBufSize()	{ return mStripeBufSize; }

private:
	/* not supported */
	FFTEngine(FFTEngine &src);
	void operator = (const FFTEngine &);

	/* read-only parameters */
	size_t					mNumRows;
	size_t					mNumCols;
	size_t					mNumColsComplex;
	bool					mReal;
	MatrixFFTPlan			mMfftPlan;

protected:
	
	/* 
	 * an algorithm-dependent number of well-aligned stripe buffers; the size of 
	 * each, in FFTComplex elements, is mStripeBufSize. 
	 * Stripe buffers created (in allocStripeBufs()) and freed
	 * by base class. 
	 */
	unsigned				mNumStripeBufs;
	size_t					mStripeBufSize;
	FFTComplex				**mStripeBuf;
	FFTComplex				**mStripeBufFree;	// unaligned, to be freed
	
	/* timing support */
	CFAbsoluteTime			mStartTime;
	CFAbsoluteTime			mEndTime;
	
	MFFTReturn allocStripeBufs(
		unsigned	numStripeBufs, 
		size_t		stripeBufSize,
		size_t		alighnSize);		// alignmnet in bytes
};

/* 
 * Trivial error class, the only exception thrown by any of the above classes, and 
 * then only from their constructors. All other errors are reported by return 
 * values.
 */
class FFTError
{
public:
	static void throwMe(MFFTReturn mrtn);
	MFFTReturn err()						{ return mMrtn; }
	~FFTError()							{}
private:
	FFTError(MFFTReturn mrtn)			: mMrtn(mrtn) {}	
	MFFTReturn mMrtn;
	
	/* not supported */
	FFTError();
};

/*** debugging ***/

extern const char *fftEngineOps[];

#define CE_LOG_OPS		0
#if		CE_LOG_OPS
#define ceprintf(args...)	printf(args)
#else
#define ceprintf(args...)
#endif


#endif	/* _FFT_ENGINE_H_ */
