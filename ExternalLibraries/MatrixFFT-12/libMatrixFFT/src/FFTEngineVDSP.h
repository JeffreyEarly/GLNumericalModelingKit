/*	File: FFTEngineVDSP.h 
	
	Description:
		Accelerate.framework-based FFTEngine.
	
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
 * FFTEngineVDSP.h - Accelerate.framework-based FFTEngine.
 *			  
 * Created 10/20/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_FFT_ENGINE_VDSP_H_
#define _FFT_ENGINE_VDSP_H_

#include "FFTEngine.h"
#include <libMatrixFFT/vdspUtils.h>

/* 
 * Accelerate-based FFTEngine 
 */
class FFTEngineVDSP : public FFTEngine
{
public:
	FFTEngineVDSP(
		MatrixFFTPlan		mfftPlan,
		unsigned			dims,
		size_t				numRows,
		size_t				numCols,
		bool				real);

	virtual ~FFTEngineVDSP();
	
	/*** basic public FFTEngine methods ***/
	
	/* basic "complex FFT a contiguous set of rows" */
	virtual MFFTReturn fftRows(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		FFTFloat			normFact,
		size_t				startRow, 
		size_t				rowsToProcess);
		
	/* complex 1-D forward twist, then FFT a contiguous set of rows */
	virtual MFFTReturn twistPlusFftRow(
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		size_t				startRow, 
		size_t				rowsToProcess);
		
	/* (inverse) complex FFT, then inverse twist a contiguous set of rows */
	virtual MFFTReturn fftRowPlusInvTwist(
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		size_t				startRow, 
		size_t				rowsToProcess);
		
	/* complex FFT a contiguous set of columns */
	virtual MFFTReturn fftColumns(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		FFTFloat			normFact,
		size_t				startCol, 
		size_t				colsToProcess);

	/* real-to-complex FFT a contiguous set of rows with optional normalization */
	virtual MFFTReturn fftRowsReal(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		FFTFloat			normFact,
		size_t				startRow, 
		size_t				rowsToProcess);
		
	/* 1-D real twist */
	virtual MFFTReturn realTwist(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		size_t				startRow, 
		size_t				rowsToProcess);

	/* determine number of rows or columns need to process given op */
	virtual size_t totalEltsForOp(
		FFTEngineOp			op);

	/*
	 * Determine optimal alignment for specified op. 
	 */
	virtual size_t alignmentForOp(
		FFTEngineOp			op);
	
	/*** end of FFTEngine methods ***/

	/* low-level atom called by thread callout (that's why it's public) */
	MFFTReturn fftRowsAtom(
		bool				forward,		
		FFTComplex			*inBuf,
		FFTComplex			*outBuf,
		FFTFloat			normFact,
		size_t				startRow, 
		size_t				rowsToProcess,
		unsigned			log2Size);		/* can be numRows() when called by fftColumnsAtom() */
		
	unsigned log2NumCols()			{ return mLog2NumCols; }
	unsigned log2NumRows()			{ return mLog2NumRows; }
private:
	/* not supported */
	FFTEngineVDSP(FFTEngineVDSP &src);
	void operator = (const FFTEngineVDSP &);

	FFT_Setup				mFftSetup;		/* same used for rows & columns */
	
	/* Since vDSP works in log2 numbers.... */
	unsigned				mLog2NumRows;
	unsigned				mLog2NumCols;
	unsigned				mLog2NumColsComplex;
	
	/* column stripe size in columns */
	size_t					mColStripeSize;
	
	#if		!FFT_SPLIT_COMPLEX
	/* 
	 * Buffers for conversion to/from vDSP format, when
	 * the library's native format is interleaved.
	 * Size is the larger of (numRows, numCols).
	 * Contents of mVdspBuf are vector-aligned; we have to free
	 * the pointers in mVdspBufFree.
	 */
	vDSPComplex				mVdspBuf;
	vDSPComplex				mVdspBufFree;
	#endif	/* FFT_SPLIT_COMPLEX */

};

#endif	/* _FFT_ENGINE_VDSP_H_ */
