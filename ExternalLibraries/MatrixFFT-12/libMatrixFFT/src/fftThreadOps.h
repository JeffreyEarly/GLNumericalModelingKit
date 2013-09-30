/*	File: fftThreadOps.h  
	
	Description:
		MatrixFFT-specific definitions for ThreadPool module.
	
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
 * fftThreadOps.h - MatrixFFT-specific definitions for ThreadPool module.
 *			  
 * Created Oct. 13 2008.
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_FFT_THREAD_OPS_H_
#define _FFT_THREAD_OPS_H_

#include <libMatrixFFT/fftPrecision.h>
#include <libMatrixFFT/NUFFT.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Split complex transpose */
typedef struct {
	const FFTFloat		*src;			/* not used for FRTO_TransposeIP */
	FFTFloat			*dst;		
	/* total size of src */
	size_t				numRows;
	size_t				numCols;		/* not used for FRTO_TransposeIP */
	/* portion to move */
	size_t				startRow;
	size_t				rowsToMove;
} TP_TransposeSplit;

/* Interleaved complex transpose */
typedef struct {
	const FFTComplex	*src;			/* not used for FRTO_TransposeIP */
	FFTComplex			*dst;		
	/* total size of src */
	size_t				numRows;
	size_t				numCols;		/* not used for FRTO_TransposeIP */
	/* portion to move */
	size_t				startRow;
	size_t				rowsToMove;
} TP_TransposeInt;

/* not currently used */
typedef struct {
	FFTComplex			*buf;		
	FFTFloat			scaleFact;
	
	/* total size of buf */
	size_t				totalElts;
	/* portion to scale */
	size_t				startElt;
	size_t				eltsToScale;
} TP_ComplexScale;

/* Nonuniform brute force DFT */
typedef struct {
    /* entire signal */
    FFTComplex          *inBuf;     
    FFTFloat            *tau;
    size_t              D;      
    FFTComplex          *outBuf;
    
    /* area of current thread's work */
    size_t              startElt;
    size_t              eltsToProcess;
} TP_NUDFT;

/* Inner NUFFT loop */
typedef struct {
    NUFFTPlan           nuFftPlan;
    unsigned            m;
    FFTFloat            oneOverMFact;        /* optimization, caller calculates this */
    FFTFloat            negTwoPiOverN;       /* ditto */
    FFTComplex          *outBuf;
    size_t              startK;
    size_t              kToProcess;
} TP_NUFFT_Loop;

/* NUFFT theta scale */
typedef struct {
    FFTComplex          *x;                 /* starting at inBuf[0] */
    FFTFloat            *theta;             /* starting at theta[0], not *our* startJ */
    size_t              startJ;
    size_t              jToProcess;
} TP_NUFFT_Theta;
    
/*
 * Indices associated with FFTEngine operations. Used in ThreadPool callouts,
 * workRatio arrays, maybe others. 
 */
typedef enum {
	FEO_FftRows			= 0,
	FEO_TwistFft,		/* 1-D complex twist plus FFT rows */
	FEO_FftInvTwist,	/* 1-D inverse complex FFT plus inverse twist */
	FEO_FftCols,
	FEO_FftRowsReal,
	FEO_RealTwist,		/* 1-D real twist */
	FEO_NumOps			/* for size - not an op */
} FFTEngineOp;

/* FFTEngine ops - 'startElt' is startRow or startCol as appropriate. etc. */
class FFTEngine;

typedef struct {
	FFTEngine			*engine;
	MatrixFFTPlan		mfftPlan;
	FFTEngineOp			op;
	bool				forward;
	FFTComplex			*inBuf;			/* [0, 0] in input */
	FFTComplex			*outBuf;		/* [0, 0] in output */
	size_t				startElt;		/* startRow, startCol */
	size_t				eltsToProcess;	/* rows or columns to process */
	FFTFloat			normFact;	
} TP_Engine;

/* Other task info structs here */

union TP_TaskUnion_U {
	TP_TransposeSplit	transS;
	TP_TransposeInt		transI;
	TP_Engine			engine;
	TP_ComplexScale		scale;
    TP_NUDFT            nuDft;
    TP_NUFFT_Loop       nuFftLoop;
    TP_NUFFT_Theta      nuFftTheta;
};


#ifdef __cplusplus
}
#endif

#endif	/* _FFT_THREAD_OPS_H_ */
