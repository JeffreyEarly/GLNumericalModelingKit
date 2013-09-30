/*	File: MatrixFFTPlan.h 
	
	Description:
		Private MatrixFFTPlan definitions.
	
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
 * MatrixFFTPlan.h - Private MatrixFFTPlan definitions.
 *			  
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_CL_MATRIX_PLAN_H_
#define _CL_MATRIX_PLAN_H_

#include <libMatrixFFT/MatrixFFT.h>
#include "fftThreadOps.h"
#include "ThreadPool.h"
#include <libMatrixFFT/vdspUtils.h>
#include <stdint.h>

#ifdef	__cplusplus
extern "C" {
#endif

/*
 * Algorithm-specific callouts. 
 */

/* execute */
typedef MFFTReturn (*mfftExecFcn)(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,
	FFTComplex			*outBuf);

/* transpose */
typedef MFFTReturn (*mfftTransFcn)(
	MatrixFFTPlan		mfftPlan,
	bool				forward,		
	bool				input, 
	const FFTComplex	*inBuf,
	FFTComplex			*outBuf);

/* 
 * Description of sin/cos lookup tables. 
 */
typedef enum {
	STT_None,			// no sin tables, e.g. 2-D
	STT_Standard,		// fully populated sin tables
	STT_Optimized,		// "thin" sin and cos
	STT_External		// managed by owner of plan, e.g. 1-D real managing its
						//   1-D complex subplan
} FFTSineTableType;

class FFTEngine;

struct MatrixFFTPlanStruct
{	
	/* algorithm-specific callouts */
	mfftExecFcn				execFcn;
	mfftTransFcn			transFcn;
	
	/* These are configured during the call to initFcn. */
	MFFTFormat				inputFormat;
	MFFTFormat				outputFormat;
	
	/* 
	 * dimensions of app-level data, viewed as rectangle appropriate
	 * to algorithm in use
	 */
	size_t					numRows;
	size_t					numCols;
	size_t					numColsComplex;	/* numCols/2 if realFft */
	size_t					N;				/* size of tables, numRows*numCols, a.k.a. 2*pi */
	
	/* and the log2 versions */
	size_t					nRow;		
	size_t					nCol;
	size_t					nColComplex;
	
	bool					realFft;
	
	/* 
	 * Aux buffer for copies (e.g. asked to do an in-place transpose for 
	 * an algorithm that must be OOP). Aligned to submatrix size. 
	 */
	FFTComplex				*auxBuf;
	size_t					auxBufSize;		/* size in FFTComplex elements */
	FFTComplex				*auxBufFree;
	
	/* variable number of Host engines, normally = number of host cores */
	unsigned				numHostEngines;
	FFTEngine				**fftEngines;
	
	/* Allocated L2 cache per engine */
	uint32_t				l2CachePerEngine;
	
	/* 
	 * Parameters for the sine lookup tables.
	 * 2-D FFTs do not require sine tables. 
	 */
	FFTSineTableType		sinTableType;
	
	size_t					sinCosNMask;	/* N-1, for quick (ab) mod N */
	size_t					piOver2;		/* pi/2,   N/4 */
	size_t					pi;				/* pi,     N/2 */
	size_t					threePiOver2;	/* 3*pi/2, 3N/4 */

	/* 
	 * For STT_External only: shift for sine lookup when signal size is smaller
	 * than the sine table.
	 */
	unsigned				sineShift;
	
	/* STT_External: location of actual sin table we're using */
	MatrixFFTPlan			externSineTable;
	
	/* for optimized sin lookup; sinPeriod != 0 indicates optimized lookup */
	size_t					sinPeriod;
	size_t					valuesPerRow;
	unsigned				sinPeriodShift;	/* sinPeriod = 1 << sinPeriodShift */
	
	/* lookup table */
	FFTFloat				*sine;

	/* only for optimized - cos not derivable from sin */
	FFTFloat				*cosine;

	/* 
	 * One pthread per FFTEngine.
	 * We actually *use* threadPoolP, which normally points to our
	 * threadPool, except in the case where this is a subplan (e.g. 
	 * a 1-D complex plan used by 1-D real FFTs) in which case
	 * threadPoolP points to the superplan's threadPool.
	 */
	TP_ThreadPool			threadPool;
	TP_ThreadPool			*threadPoolP;
	
	/* one vDSP setup, shared between all engines */
	FFT_Setup				vdspSetup;
	
	/*
	 * Additional normalization factor, used when normalization is performed by 
	 * a subplan and the subplan's normalization factor differs from the 
	 * superplan's. 
	 * 1-D real needs a normalization factor of 1/2N where N is the number
	 * of real samples. It uses a 1-D complex subplan which has a 
	 * normalization factor of 1/N where N is the number of *complex* samples.
	 * Thus in this case the 1-D complex subplan would have an addNormFact
	 * of 0.25. 
	 * If addNormFact is 0.0 no additional normalization is required. 
	 */
	FFTFloat				addNormFact;
	
	#if		!FFT_SPLIT_COMPLEX
	/* 
	 * Buffer for conversion to/from vDSP format, used by 2-D real,
     * NUFFT, ops small enough to use raw vDSP. 
     * and Vector-aligned, freed via vdspBufFree.
	 */
	vDSPComplex				vdspBuf;
	vDSPComplex				vdspBufFree;
	#endif	/* FFT_SPLIT_COMPLEX */
	
	/* 
	 * Subplan, used for 1-D real FFT. 
	 */
	MatrixFFTPlan			subPlan;
    
};

/* 
 * Flags passed to mfftPlanCreateCommon to inhibit creation of 
 * various MatrixFFTPlan components.
 */
#define MFFT_PLAN_OPT_NONE          0x0000
#define MFFT_PLAN_OPT_NO_VDSP       0x0001  /* No vdspSetup */
#define MFFT_PLAN_OPT_NO_THREADS    0x0002  /* No ThreadPool */
#define MFFT_PLAN_OPT_NO_ENGINES    0x0004  /* No engines (vDSP only) */

/* 
 * Common MatrixFFTPlan creation, in in MatrixFFT.cpp.
 */
extern MFFTReturn mfftPlanCreateCommon(
	unsigned			dims,			// 1 or 2
	FFTSineTableType	sinTableType,
	bool				real,
	uint32_t			optFlags,
	/* rectangle dimensions, even for 1-D */
	unsigned			nRow,			// two dimension : N = 2 ^ (nRow + nCol)
	unsigned			nCol,
	unsigned			numThreads,
	size_t				sinPeriod,		// only for STT_Optimized
    uint32_t            options,        // MFFT_PLAN_OPT_*, above
	MatrixFFTPlan		*mfftPlanRtn);	// RETURNED

/* 
 * Algorithm-specific MatrixFFTPlan creation declarations.
 * On successful return, the following MatrixFFTPlan fields
 * are valid:
 * 
 * -- inputFormat
 * -- outputFormat
 * -- execFcn
 * -- transFcn
 */
 
/* 2-D Complex, in fft2DComplex.cpp */
extern MFFTReturn mfftCreatePlan2DComplex(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlan);		/* RETURNED */

/* 1-D complex, in fft1DComplex.cpp */
extern MFFTReturn mfftCreatePlan1DComplex(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	bool				externalSineTable,	/* don't create sine tables */
	bool				disableThreads,		/* skip ThreadPool init */
    bool                expectColOut,       /* Don't optimize by using vDSP */
	MatrixFFTPlan		*mfftPlan);			/* RETURNED */

/* 1-D Real, in fft1DReal.cpp */
extern MFFTReturn mfftCreatePlan1DReal(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlanRtn);

/* 2-D Real, in fft2DReal.cpp */
extern MFFTReturn mfftCreatePlan2DReal(
	unsigned			*n,			
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlanRtn);

#ifdef __cplusplus
}
#endif

#endif	/* _CL_MATRIX_PLAN_H_ */
