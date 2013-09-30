/*	File: fft1DComplexTwist.cpp 
	
	Description:
		1-D Complex FFT twist.
	
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
 * fft1DComplexTwist.cpp - 1-D Complex FFT twist.
 *			  
 * Created 11/21/2008. 
 * Copyright 2008 by Apple, Inc. 
 */

#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "MatrixFFTPlan.h"
#include "fftPriv.h"
#include "fftSinCos.h"
#include "fftDebug.h"
#include "fftIntel.h"
#include "PolyComplex.h"

#pragma mark --- Configuration and debugging ---

#define CFT_DEBUG						0

#if		CFT_DEBUG
#define cftDebug(s...)		printf(s)
#else
#define cftDebug(s...)
#endif

/*
 * This is an implementation of Algortihm 1 in "Large-scale FFTs and convolutions on 
 * Apple hardware", available here:
 * http://images.apple.com/acg/pdf/FFTapps.pdf
 *
 *
 * The twist for this module is defined as:
 *
 *    buf[a,b] *= exp (-/+ 2 * pi * i * a * b / N) 
 *             *= cos ((-/+ 2 * pi * a * b / N) + i sin (-/+ 2 * pi * a * b / N))
 *             *= cos ((2 * pi * a * b / N)     + i sin (-/+ 2 * pi * a * b / N))
 *
 * numRows, numCols        : total size of array 
 * startRow, rowsToProcess : portion to process 
 */

#if		FFT_INTEL

/*
 * Intel, precision-independent. 
 */
static MFFTReturn fft1DTwistOpt(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	bool				forward,	
	size_t				numRows,	
	size_t				numCols,
	size_t				startRow,	
	size_t				rowsToProcess)
{
	FFTFloat 	imagSign = forward ? -1.0 : 1.0;
	FFTVector	vImagSign = FFTVectSet1(imagSign);
	size_t		lastRow = startRow + rowsToProcess;
	
	#if	DUMP_MATRIX
	#if FFT_SPLIT_COMPLEX
	FFTComplex start;
	fftComplexOffset(buf, startRow * numCols, &start);
	dumpMatrixRect("fft1DTwistOpt input", &start, rowsToProcess, numCols);
	#else
	FFTComplex *start = buf + (startRow * numCols);
	dumpMatrixRect("fft1DTwistOpt input", start, rowsToProcess, numCols);
	#endif	/* FFT_SPLIT_COMPLEX */
	#endif	/* DUMP_MATRIX */

	for(size_t row=startRow; row<lastRow; row++) {
		size_t rowOff = numCols * row;
		PolyComplex pc(buf, rowOff);
		
		FFTVector		vTempCos;
		FFTVector		vCurCos;
		FFTVector		vCurSin;
		FFTVector		vIncA;
		FFTVector		vIncB;
		FFTVectUnion	transferCos;
		FFTVectUnion	transferSin;
		unsigned		angleIndex;
		
		// set up initial sin & cos vectors
		if(mfftPlan->sinPeriod) {
			for(angleIndex = 0; angleIndex < FFT_FLOATS_PER_VECTOR; angleIndex++) {
				fftSinCosOpt(mfftPlan, row, angleIndex, 
					&transferCos.f[angleIndex], &transferSin.f[angleIndex]);
			}
		}
		else {
			for(angleIndex = 0; angleIndex < FFT_FLOATS_PER_VECTOR; angleIndex++) {
				fftSinCos(mfftPlan, row*angleIndex, 
					&transferCos.f[angleIndex], &transferSin.f[angleIndex]);
			}
		}
		
		vCurCos = transferCos.v;
		vCurSin = transferSin.v;

		// angle of increment between steps, FFT_FLOATS_PER_VECTOR steps, since
		// each vector has FFT_FLOATS_PER_VECTOR elements
		
		FFTFloat incA, incB;
		if(mfftPlan->sinPeriod) {
			fftSinCosOpt(mfftPlan, row, FFT_FLOATS_PER_VECTOR / 2, NULL, &incA);
			incA = incA*incA*2;
			fftSinCosOpt(mfftPlan, row, FFT_FLOATS_PER_VECTOR, NULL, &incB);
		}
		else {
			size_t incAngle = row * FFT_FLOATS_PER_VECTOR;
			fftSinCos(mfftPlan, incAngle / 2, NULL, &incA);
			incA = incA*incA*2;
			fftSinCos(mfftPlan, incAngle, NULL, &incB);
		}
		vIncA = FFTVectSet1(incA);
		vIncB = FFTVectSet1(incB);
		
		for (size_t col=0; col<numCols; col+=FFT_FLOATS_PER_VECTOR) {
			FFTVector vRTop;
			FFTVector vITop;

			// prefetch these
			pc.loadVect(vRTop, vITop);
									
			FFTVector vcosv = vCurCos;
			FFTVector vsinv = vCurSin;

			if(col < (numCols - FFT_FLOATS_PER_VECTOR - 1)) {
				/* Update vCurSin and vCurCos unless we're at end of row. */
				if((FFT_SIN_RECALC_COMPLEX > 0) && 
				   ((col % FFT_SIN_RECALC_COMPLEX) == 
				       (unsigned)(FFT_SIN_RECALC_COMPLEX - FFT_FLOATS_PER_VECTOR))) {
					size_t newCol = col + FFT_FLOATS_PER_VECTOR;
					for (angleIndex = 0; angleIndex < FFT_FLOATS_PER_VECTOR; angleIndex++) { 
						/*
						 * Note that we might be using a fully populated sine table even if 
						 * we're configured with FFT_SIN_RECALC_COMPLEX > 0. This happens
						 * when we're running as a subplan of a 1-D real FFT. In that
						 * case we're using the 1-D real's sine table, which is always fully
						 * populated.
						 */
						if(mfftPlan->sinPeriod) {
							fftSinCosOpt(mfftPlan, row, newCol+angleIndex, 
								&transferCos.f[angleIndex], &transferSin.f[angleIndex]);
						}
						else {
							fftSinCos(mfftPlan, row*(newCol+angleIndex), 
								&transferCos.f[angleIndex], &transferSin.f[angleIndex]);
						}
					}
					vCurCos = transferCos.v;
					vCurSin = transferSin.v;
				}
				else {
					// vTempCos = vCurCos - incA*curCos - incB*curSin;					
					vTempCos = FFTVectSub(vCurCos, FFTVectAdd(FFTVectMul(vIncA, vCurCos),
															  FFTVectMul(vIncB, vCurSin)));
					// curSin = curSin - incA*curSin + incB*curCos;
					vCurSin = FFTVectSub(vCurSin, FFTVectSub(FFTVectMul(vIncA, vCurSin),
															 FFTVectMul(vIncB, vCurCos)));
					vCurCos = vTempCos;
				}
			}
			
			// sinv *= imagSign;
			vsinv = FFTVectMul(vsinv, vImagSign);
			
			// real = (cosv * rTop) - (sinv * iTop);
			FFTVector vr = FFTVectSub(FFTVectMul(vcosv, vRTop), FFTVectMul(vsinv, vITop));
						
			// imag = (cosv * iTop) + (sinv * rTop);
			FFTVector vi = FFTVectAdd(FFTVectMul(vcosv, vITop), FFTVectMul(vsinv, vRTop));

			pc.storeVect(vr, vi);
			pc.offset(FFT_FLOATS_PER_VECTOR);
		}
	}
	#if FFT_SPLIT_COMPLEX
	dumpMatrixRect("fft1DTwistOpt output", &start, rowsToProcess, numCols);
	#else
	dumpMatrixRect("fft1DTwistOpt output", start, rowsToProcess, numCols);
	#endif
	return MR_Success;
}

#else	/* FFT_INTEL */

/* This code never runs, it's just here to compile */
static MFFTReturn fft1DTwistOpt(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	bool				forward,	
	size_t				numRows,	
	size_t				numCols,
	size_t				startRow,	
	size_t				rowsToProcess)
{
	return MR_Unsupported;
}

#endif	/* FFT_INTEL */

/* Plain vanilla, unoptimized, platform-independent twist */
static MFFTReturn fft1DTwistSmall(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	bool				forward,	
	size_t				numRows,
	size_t				numCols,
	size_t				startRow,
	size_t				rowsToProcess)
{
	RFASSERT((mfftPlan->sinTableType == STT_Standard) ||
	         (mfftPlan->sinTableType == STT_External));

	FFTFloat imagSign = forward ? -1.0 : 1.0;
	
	#if	DUMP_MATRIX
	#if FFT_SPLIT_COMPLEX
	FFTComplex start;
	fftComplexOffset(buf, startRow * numCols, &start);
	dumpMatrixRect("fft1DTwistSmall input", &start, rowsToProcess, numCols);
	#else
	FFTComplex *start = buf + (startRow * numCols);
	dumpMatrixRect("fft1DTwistSmall input", start, rowsToProcess, numCols);
	#endif	/* FFT_SPLIT_COMPLEX */
	#endif	/* DUMP_MATRIX */
	
	size_t row = startRow;
	for(size_t rowDex=0; rowDex<rowsToProcess; rowDex++, row++) {
		PolyComplex pc(buf, row * numCols);
		for(size_t col=0; col<numCols; col++) {
			FFTFloat cosv, sinv;
			fftSinCos(mfftPlan, row*col, &cosv, &sinv);
			
			sinv *= imagSign;
			FFTFloat r = (cosv * pc.real()) - (sinv * pc.imag());
			FFTFloat i = (cosv * pc.imag()) + (sinv * pc.real());
			pc.real(r);
			pc.imag(i);
			++pc;
			
		}
	}

	#if FFT_SPLIT_COMPLEX
	dumpMatrixRect("fft1DTwistSmall output", &start, rowsToProcess, numCols);
	#else
	dumpMatrixRect("fft1DTwistSmall output", start, rowsToProcess, numCols);
	#endif
	
	return MR_Success;
}

/* 
 * General purpose, public 1-D twist - use optimized version if we can, else
 * use generic version. 
 */
MFFTReturn fft1DTwist(
	MatrixFFTPlan		mfftPlan,
	FFTComplex			*buf,
	bool				forward,
	size_t				numRows,
	size_t				numCols,
	size_t				startRow,
	size_t				rowsToProcess)
{
	cftDebug("fft1DTwist startRow %llu\n", (unsigned long long)startRow);
	
	if((numRows > FFT_FLOATS_PER_VECTOR) && 
	   FFT_INTEL &&
	   !FFT_FWD_TWIST_COMPLEX_SMALL) {
		return fft1DTwistOpt(mfftPlan, buf, forward, numRows, numCols, startRow, rowsToProcess);
	}
	else {
		return fft1DTwistSmall(mfftPlan, buf, forward, numRows, numCols, startRow, rowsToProcess);
	}
}
