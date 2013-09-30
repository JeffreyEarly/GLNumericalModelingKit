/*	File: NUFFT.h 
	
	Description:
		Public API for Nonuniform FFT. 
	
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
 * NUFFT.h - Public API for Nonuniform FFT.
 *			  
 * Created 04/09/2009. 
 * Copyright 2009 by Apple, Inc. 
 */
 
#ifndef	_NU_FFT_H_
#define _NU_FFT_H_

#include <libMatrixFFT/MatrixFFT.h>
#include <stdint.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * Analagous to vDSP's FFTSetup, this is our state for a given 
 * operation. It can only be used for the operation for which it
 * was explicitly created; any change to sizes or type require a 
 * new NUFFTPlan. 
 * At this interface we have an opaque pointer. 
 */
struct NUFFTPlanStruct;
typedef struct NUFFTPlanStruct *NUFFTPlan;

/*
 * Flags passed to nufftCreatePlan() and nuFftExecute(). 
 *
 * There are three different nonuniform Fourier transforms implemented
 * herein. A given NUFFTPlan can be used for any or all of them depending
 * on the flags passed to nufftCreatePlan(). Pass exactly one of these
 * three flags to nuFftExecute to select one of the three implementations. 
 *
 * NF_Discrete    -- Brute force nonuniform Discrete Fourier Transform. Very slow;
 *                   reference for measuring the accuracy of other implementations. 
 * NF_Reference   -- Basic NUFFT; not optimized in any way. 
 * NF_Optimized   -- Current best NUFFT performer.
 */
#define NF_Discrete         0x0001
#define NF_Reference        0x0002
#define NF_Optimized        0x0004

/* 
 * Create an NUFFTPlan for nonuniform FFT. 
 *
 * Returned NUFFTPlan must be freed via nufftFreePlan(). 
 * 
 * Inputs:
 * ------
 * dims       : Dimensions. Must be 1 (one) for now. 
 * n          : Array of log2 of dimensions. For 2D, n[0] is (will be) 
 *              log2(numRows), n[1] is logn2(numColumns).
 * bits       : Number of bits of rigorous accuracy needed in output.
 * optFlags   : Option flags, see comments above. Specify at least one, and up
 *              to 3, implementations. 
 * maxThreads : Max number of pthreads to use. Specify 0 to use the default 
 *              (which is currently the number of CPU cores).
 */ 
extern MFFTReturn nufftCreatePlan(
	unsigned			dims,		
	unsigned			*n,			
    unsigned            bits,    
	uint32_t			optFlags,           /* NF_Optimized, etc. */
	unsigned			numThreads,
	NUFFTPlan           *nufftPlan);		/* RETURNED */

/*
 * Free an NUFFTPlan.
 */
extern void nufftFreePlan(
	NUFFTPlan           nuFftPlan);


/* 
 * 1-D Nonuniform FFT.
 *
 * Inputs:
 * ------
 * nuFftPlan  : created in nufftCreatePlan().
 * optFlags   : Option flags, see comments above. Specify exactly one 
 *              implementation. 
 * inBuf      : complex input samples, length 2^n (as specified in 
 *              mfftCreatePlanNUFFT()). NOTE: the NF_Optimized implementation
 *              modifies this input buffer. 
 * tau        : Times associated with inBuf, length 2^n. 
 * outBuf     : Caller-allocated output buffer, length 2^n complex samples.
 *              Note: implementation NF_Reference can perform NUFFT 
 *              in-place; other implementations must be OOP. 
 */
extern MFFTReturn nuFftExecute(
	NUFFTPlan       nuFftPlan,
	uint32_t        optFlags,           /* NF_Optimized, etc. */
	FFTComplex		*inBuf,	
    FFTFloat        *tau,   
    FFTComplex      *outBuf); 

#ifdef __cplusplus
}
#endif

#endif	/* _NU_FFT_H_ */
