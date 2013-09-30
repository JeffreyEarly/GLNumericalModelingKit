/*	File: fft1DRealTwist.cpp 
	
	Description:
		1-D Real-signal FFT twist
	
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
 * fft1DRealTwist.cpp - 1-D real-signal FFT twist
 *			  
 * Created 1/5/2009. 
 * Copyright 2009 by Apple, Inc. 
 */

#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/complexBufUtils.h>
#include "MatrixFFTPlan.h"
#include "fftPriv.h"
#include "fftSinCos.h"
#include "PolyComplex.h"
#include "fftIntel.h"

/* 
 * This is an implementation of Algortihm 2 in "Large-scale FFTs and convolutions on 
 * Apple hardware", available here:
 * http://images.apple.com/acg/pdf/FFTapps.pdf
 *
 *
 * The twist operation for 1-D real signal FFT is performed on the output of 
 * a 1-D complex FFT in column order. The twist is defined as follows. 
 * 
 * Note "-/+" means "-" for forward, "+" for reverse. "Q@[k]" is the complex
 * conjugate of Q[k].
 *
 * 2P[k] = Q[k]+Q@[N/2 - k] -/+ (i * exp(-/+2 * pi * i * k / N) * (Q[k]-Q@[N/2 - k]))
 *       = Q[k]+Q@[N/2 - k] -/+ (i * (cos(2*pi*j/N) -/+ (i * sin(2*pi*j/N))) * (Q[k]-Q@[N/2 - k]))
 *	     = Q[k]+Q@[N/2 - k] -/+ ((+/-sin(2*pi*j/N) + (i * cos(2*pi*j/N))) * (Q[k]-Q@[N/2 - k]))
 *
 *
 * The twist implementation operates on consecutive memory locations, which 
 * do not represent consecutive values of k (or of N/2 - k). To illustrate
 * what's happening, take this array, 4 rows, 8 columns, shown in column order:
 * 
 * Q0 Q4 Q8  Q12 Q16 Q20 Q24 Q28
 * Q1 Q5 Q9  Q13 Q17 Q21 Q25 Q29
 * Q2 Q6 Q10 Q14 Q18 Q22 Q26 Q30
 * Q3 Q7 Q11 Q15 Q19 Q23 Q27 Q31
 * 
 * When you traverse row 0 in linear order (starting with element 1, since
 * element 0 is a special case that's handled out-of-line) you're dealing with
 * 
 * Q[k]     = Q4  Q8  Q12 Q16 Q20 Q24 Q28
 * Q[N/2-k] = Q28 Q24 Q20 Q16 Q12 Q8  Q4
 * 
 * The the Q[N/2-k] values are in that same row, backwards. And Q16 is its 
 * own Q[N/2-k].
 * 
 * Now take row 1:
 * 
 * Q[k]     = Q1  Q5  Q9  Q13 Q17 Q21 Q25 Q29   -- row 1
 * Q[N/2-k] = Q31 Q27 Q23 Q19 Q15 Q11 Q07 Q03   -- row 3 backwards
 * 
 * Row 2:
 * 
 * Q[k]     = Q2  Q6  Q10 Q14 Q18 Q22 Q26 Q30   -- row 2
 * Q[N/2-k] = Q30 Q26 Q22 Q18 Q14 Q10 Q06 Q02   -- row 2 backwards
 * 
 * Row 3:
 * 
 * Q[k]     = Q3  Q7  Q11 Q15 Q19 Q23 Q27 Q31   -- row 3
 * Q[N/2-k] = Q29 Q25 Q21 Q17 Q13 Q09 Q05 Q01   -- row 1 backwards
 *
 * So, row 0 is a special case. After that, processing Q[k] and Q[N/2-k] 
 * involves contiguous memory locations, incrementing for Q[k] and 
 * decrementing for Q[N/2-k] without regard to row or column boundaries.
 * The values of k and N/2-k do, however, have discontinuities at
 * row boundaries; within a row they increment and decrement by 
 * numRows.
 *
 * Subsequent to processing row 0, when processing Q[k] values in 
 * row n, the corresponding Q[N/2-k] values are in row numRows-n. 
 *
 * Also note the special case for row numRows/2, in which the 
 * Q[k] and Q[N/2-k] values are in the same row. 
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * When calculating sin and cos for successive values of k and N/2-k, we
 * use the following:
 *
 *		cos(x+y) = cos(x)*cos(y) - sin(x)*sin(y)
 *		sin(x+y) = sin(x)*cos(y) + sin(y)*cos(x)
 *
 * While traversing one row sequentially in memory, the difference
 * between successive k values is numRows ('h'); that's the y in 
 * the above formulae, and it's constant, thus so are these:
 * 
 *		d = cos(2 * pi * h / N) 
 *		t = sin(2 * pi * h / N)
 *
 * So given sin(Q[k]) and cos(Q[k]):
 *
 *		cos(Q[k+numRows]) = cos(Q[k])*d - sin(Q[k])*t
 *		sin(Q[k+numRows]) = sin(Q[k])*d + cos(Q[k])*t
 *
 * We process Q[N/2-k] in reverse order; for that we use
 *
 *		cos(x-y) = cos(x)*cos(-y) - sin(x)*sin(-y)
 *				 = cos(x)*cos(y)  + sin(x)*sin(y)
 *		sin(x-y) = sin(x)*cos(-y) + sin(-y)*cos(x)
 *				 = sin(x)*cos(y)  - sin(y)*cos(x)
 *
 * Thus given sin(Q[h]) and cos(Q[k]):
 *
 *		cos(Q[k-numRows]) = cos(Q[k])*d + sin(Q[k])*t
 *		sin(Q[k-numRows]) = sin(Q[k])*d - cos(Q[k])*t
 *
 *************************************************************************
 *
 * This 1-D real twist implementation is precision-independent as well
 * as independent of split/interleaved complex format. 
 */

/* Force fft1DRealTwistSmall() to debug the rest */
#define FFT_FWD_TWIST_REAL_SMALL		0

/* 
 * Although we use a fully populated sine lookup table, because our caller
 * shares it with its 1-D complex subplan, we still use the sin recalculation
 * technique to minimize sine lookups (i.e. we still save memory bandwidth,
 * just not overall memory usage). 
 * Sine reclaculation period - <0 to disable
 */
#if		FFT_DOUBLE_PREC
#define FFT_SIN_RECALC_REAL			(128)
//#define FFT_SIN_RECALC_REAL		(-1)
#else
#define FFT_SIN_RECALC_REAL			(32)
//#define FFT_SIN_RECALC_REAL		(-1)
#endif


/* log every twist op - very verbose */
#define DEBUG_TWIST					0
#if		DEBUG_TWIST
#define tprintf(args...)			printf(args)
#else
#define tprintf(args...)
#endif	/* DEBUG_TWIST */

#pragma mark --- Scalar twist ---

/*
 * Platform independent 1-D real FFT twist operation.
 */
static MFFTReturn fftReal1DTwistScalar(
	MatrixFFTPlan			mfftPlan,
	const FFTComplex		*src,
	FFTComplex				*dst,
	size_t					numRows,
	size_t					numCols,		// in complex elements
	bool					forward,
	
	/* The rows we actually process. */
	size_t					startRow,
	size_t					rowsToProcess)
{
	size_t nOver2 = numRows * numCols;	// N/2, N is the size of the original real array 
	
	/* working pointers for source */
	PolyComplex kIn((FFTComplex *)src);		// Q[k]
	PolyComplex nmkIn((FFTComplex *)src);	// Q[N/2-k]
	
	/* working pointers for destination */
	PolyComplex kOut(dst);					// P[k]
	PolyComplex nmkOut(dst);				// P[N/2-k]

	/* Current values of Q[k] and Q[N/2-k] */
	FFTFloat qkR;
	FFTFloat qkI;
	FFTFloat qnmkR;
	FFTFloat qnmkI;

	FFTFloat d;
	FFTFloat t;
	fftSinCos(mfftPlan, numRows, &d, &t);
	
	/* inputs currently at src; outputs currently at dst */
	if(startRow == 0) {
		/* k=0, Q[k] at start of data */
		/* N/2-k, end of row 0 */
		size_t offset = numCols - 1;
		nmkIn.offset(offset);
		nmkOut.offset(offset);
		
		/* 
		 * DC      = Q[0].real + Q[0].imag
		 * Nyquist = Q[0].real - Q[0].imag
		 *
		 * And we multiply these by 2 to match Accelerate FFT...
		 */
		qkR = kIn.real();
		qkI = kIn.imag();
		++kIn;
		
		if(forward) {
			kOut.real((qkR + qkI) * 2.0);	/* DC */
			kOut.imag((qkR - qkI) * 2.0);	/* Nyquist */
		}
		else {
			kOut.real(qkR + qkI);			/* DC */
			kOut.imag(qkR - qkI);			/* Nyquist */
		}
		++kOut;	
	}
	else {
		/* start at arbitrary row */
		size_t topOffset = startRow * numCols; 
		kIn.offset(topOffset);
		kOut.offset(topOffset);
		
		/* 
		 * k at Row 1, col 0 corresponds to N-k at end of array.
		 * row R reflects on row numRows - r.
		 */
		size_t endRow    = numRows - startRow;
		size_t endOffset = (endRow * numCols) + numCols - 1;
		nmkIn.offset(endOffset);
		nmkOut.offset(endOffset);
	}
	
	/* 
	 * Remainder:
	 *
	 * Note "-/+" means "-" for forward, "+" for reverse.
	 *
	 * 2P[k] = Q[k]+Q@[N/2 - k] -/+ (i * exp(-/+2 * pi * i * k / N) * (Q[k]-Q@[N/2 - k]))
     *       = Q[k]+Q@[N/2 - k] -/+ (i * (cos(2*pi*j/N) -/+ (i * sin(2*pi*j/N))) * (Q[k]-Q@[N/2 - k]))
	 *	     = Q[k]+Q@[N/2 - k] -/+ ((+/-sin(2*pi*j/N) + (i * cos(2*pi*j/N))) * (Q[k]-Q@[N/2 - k]))
	 * 
	 * let s = sin(2*pi*k/N)
     *     c = cos(2*pi*k/N)
	 *     Z = Q[k] + Q@[N/2 - k]
	 *     Y = Q]k] - Q@[N/2 - k]
	 *  
     * Then 
     * 
     * 2P[k] = Z -/+ ((+/-s + c*i) * Y)
	 */
	size_t endRow = startRow + rowsToProcess;
	for(size_t row=startRow; row<endRow; row++) {
		size_t startCol;
		size_t rowNmk;		/* row for Q[N/2 - k] */

		/* k, N/2-k */
		size_t k;
		size_t nmk;
	
		if(row == 0) {
			/* 
			 * Special case for first row - smaller by 1, and it reflects on itself.
			 */
			startCol = 1;
			rowNmk = 0;
			k = numRows;
		}
		else {
			/* subsequent rows */
			startCol = 0;
			rowNmk = numRows - row;
			k = row;
		}
		nmk = nOver2 - k;
	

		/* 
		 * Starting cos, sin for each row
		 * rowSinK is the running, actual value of sin[k]; sinK gets
		 * multiplied by -1 for reverse FFT. 
		 */
		FFTFloat cosK;
		FFTFloat rowSinK;
		fftSinCos(mfftPlan, row, &cosK, &rowSinK);
		FFTFloat sinK = rowSinK;
		
		/*
		 * For N/2-k we start at the end of its row.
		 */
		FFTFloat cosNmk;
		FFTFloat rowSinNmk;
		fftSinCos(mfftPlan, nmk, &cosNmk, &rowSinNmk);
		FFTFloat sinNmk = rowSinNmk;
		
		/*
		 * Two rows reflect on themselves and require less than numCols
		 * inner loops. 
		 */
		size_t lastCol = numCols - 1;
		if(row == 0) {
			/* Last element is numCols/2 + 1, which reflects upon itself */
			lastCol = numCols >> 1; 
		}
		else if(row == (numRows >> 1)) {
			/* middle row, even number of elements, process half of it. */
			lastCol = (numCols >> 1) - 1;
		}
		tprintf("row %u startCol %u lastCol %u\n",
			(unsigned)row, (unsigned)startCol, (unsigned)lastCol);
			
		for(size_t col=startCol; col<=lastCol; col++) {
		
			/*
			 * Snag current Q[k] and Q[N/2-k]
			 */
			qkR = kIn.real();
			qkI = kIn.imag();
			++kIn;	
			qnmkR = nmkIn.real();
			qnmkI = nmkIn.imag();
			--nmkIn;
			
			/* 
			 * sinK = +/- sin(2*pi*k/N)
			 * cosK =     cos(2*pi*k/N)
			 */
			if(col > 0) {
				/* for col=0 we already know sin/cos, exactly */
				if((FFT_SIN_RECALC_REAL > 0) && 
					/* re-sync to actual sin/cos periodically */
				   ((col % FFT_SIN_RECALC_REAL) == 0)) {
					fftSinCos(mfftPlan, k, &cosK, &rowSinK);
					sinK = rowSinK;
					fftSinCos(mfftPlan, nmk, &cosNmk, &rowSinNmk);
					sinNmk = rowSinNmk;
				}
				else {
					/* Infer new sin/cos */
					FFTFloat f = (cosK * d) - (rowSinK * t);
					rowSinK = sinK = (rowSinK * d) + (cosK * t);
					cosK = f;
					
					if(!((row == 0) && (col == 1))) {
						/* 
						 * Skip this for special case of first element we're
						 * processing; we already know sin and cos exactly 
						 */
						f = (cosNmk * d) + (rowSinNmk * t);
						rowSinNmk = sinNmk = (rowSinNmk * d) - (cosNmk * t);
						cosNmk = f;
					}
				}
			}
			
			if(!forward) {
				sinK = -sinK;
				sinNmk = -sinNmk;
			}
			
			tprintf("k %u row %u col %u sin %f cos %f\n",
				(unsigned)k, (unsigned)row, (unsigned)col,
				sinK, cosK);

			/* 
			 * Z = Q[k] + Q@[N/2 - k]
			 * Y = Q]k] - Q@[N/2 - k]
			 */
			FFTFloat Zr, Zi, Yr, Yi;
			Zr = qkR + qnmkR;
			Zi = qkI - qnmkI;
			Yr = qkR - qnmkR;
			Yi = qkI + qnmkI;
			tprintf("   Q[k] = [%f, %f]  Q[N-k] = [%f, %f]\n", qkR, qkI, qnmkR, qnmkI);
			tprintf("   Z = [%f, %f];  Y = [%f, %f]\n", Zr, Zi, Yr, Yi);
			
			/* 
			 * prod = (sinK + cosK*i) * Y
			 */
			FFTFloat prodR, prodI;
			prodR = (sinK * Yr) - (cosK * Yi);
			prodI = (sinK * Yi) + (cosK * Yr);
			tprintf("   (sinK + cosK*i) * Y = [%f, %f]\n", prodR, prodI);
			
			
			/* 
			 * 2P[k] = Z -/+ ((s + c*i) * Y)
			 *       = Z -/+ prod
			 */
			if(forward) {
				kOut.real(Zr - prodR);
				kOut.imag(Zi - prodI);
			}
			else {
				kOut.real(Zr + prodR);
				kOut.imag(Zi + prodI);
			}
			tprintf("   2P[k] = [%f, %f]\n", kOut.real(), kOut.imag());
			++kOut;	
			
			if((row == 0) && (col == (numCols >> 1))) {
				/*
				 * Special case: "middle" element in row 0, reflects on 
				 * itself, we're done with the row. 
				 */
				break;
			}
			
			/* 
			 * Repeat for P[N/2 - k]; 'k' and 'nmk' are swapped
			 * 
			 * Z = Q[k] + Q@[N/2 - k] 
			 * Y = Q]k] - Q@[N/2 - k]
			 */
			tprintf("k %u row %u col %u sin %f cos %f\n",
				(unsigned)nmk, (unsigned)rowNmk, (unsigned)(numCols - col - 1),
				sinNmk, cosNmk);
			// Zr = qkR + qnmkR;		(we already have this)
			Zi = qnmkI - qkI;
			Yr = qnmkR - qkR;
			// Yi = qkI + qnmkI;		(we already have this)
			tprintf("   Q[k] = [%f, %f]  Q[N-k] = [%f, %f]\n", qnmkR, qnmkI, qkR, qkI);
			tprintf("   Z = [%f, %f];  Y = [%f, %f]\n", Zr, Zi, Yr, Yi);
			
			/* 
			 * prod = (sinK + cosK*i) * Y
			 */
			prodR = (sinNmk * Yr) - (cosNmk * Yi);
			prodI = (sinNmk * Yi) + (cosNmk * Yr);
			tprintf("   (sinK + cosK*i) * Y = [%f, %f]\n", prodR, prodI);
			
			
			/* 
			 * 2P[k] = Z -/+ ((s + c*i) * Y)
			 *       = Z -/+ prod
			 */
			if(forward) {
				nmkOut.real(Zr - prodR);
				nmkOut.imag(Zi - prodI);
			}
			else {
				nmkOut.real(Zr + prodR);
				nmkOut.imag(Zi + prodI);
			}
			tprintf("   2P[k] = [%f, %f]\n", nmkOut.real(), nmkOut.imag());
			--nmkOut;
			
			/* these increments get discarded at the end of a row */
			k   += numRows;
			nmk -= numRows;
		}	/* for col */
	
		if((row == 0) && (rowsToProcess > 1)) {
			/* 
			 * now set pointers for processing Q[1] and P[1], the elements at the 
			 * start of row 1. End pointers point to the end of array, which is  
			 * Q[N/2-k]. 
			 */
			kIn.set((FFTComplex *)src, numCols);
			kOut.set(dst, numCols);
			size_t offset = nOver2 - 1;
			nmkIn.set((FFTComplex *)src, offset);
			nmkOut.set(dst, offset);
		}
	}	/* for row */
	return MR_Success;
}

#pragma mark --- Vector twist ---

#if		FFT_INTEL

/*
 * Intel optimized 1-D real FFT twist operation.
 */
static MFFTReturn fftReal1DTwistOpt(
	MatrixFFTPlan			mfftPlan,
	const FFTComplex		*src,
	FFTComplex				*dst,
	size_t					numRows,
	size_t					numCols,		// in complex elements
	bool					forward,
	
	/* The rows we actually process. */
	size_t					startRow,
	size_t					rowsToProcess)
{
	size_t nOver2 = numRows * numCols;	// N/2, N is the size of the original real array 
	
	/* Current values of Q[k] and Q[N/2-k] */
	FFTVector qkRV;
	FFTVector qkIV;
	FFTVector qnmkRV;
	FFTVector qnmkIV;

	FFTVector minus1V = FFTVectSet1(-1.0);
	MFFTReturn mrtn;
	
	/*
	 * First off, if we've been asked to process row 0, farm that out to the scalar
	 * implementation. Row 0 is not amenable to vector processing: its elements start
	 * on an odd boundary and have an odd length.
	 */
	if(startRow == 0) {
		mrtn = fftReal1DTwistScalar(mfftPlan, src, dst, numRows, numCols, forward, 0, 1);
		if(mrtn) {
			return mrtn;
		}
		startRow++;
		if(--rowsToProcess == 0) {
			return MR_Success;
		}
	}
	
	/* 
	 * Set up constant 'y' increment angle
	 */
	FFTFloat d;
	FFTFloat t;
	fftSinCos(mfftPlan, numRows * FFT_FLOATS_PER_VECTOR, &d, &t);
	FFTVector dV = FFTVectSet1(d);
	FFTVector tV = FFTVectSet1(t);
	
	/* working pointers */
	size_t topOffset = startRow * numCols; 
	PolyComplex kIn((FFTComplex *)src, topOffset);	// Q[k]
	PolyComplex kOut(dst, topOffset);				// P[k]
	
	/* 
	 * k at Row 1, col 0 corresponds to N-k at end of array.
	 * row R reflects on row numRows - r.
	 */
	size_t endRow    = numRows - startRow;
	size_t endOffset = (endRow * numCols) + numCols - FFT_FLOATS_PER_VECTOR;
	PolyComplex nmkIn((FFTComplex *)src, endOffset);// Q[N/2-k]
	PolyComplex nmkOut(dst, endOffset);				// P[N/2-k]
	
	/* 
	 * Note "-/+" means "-" for forward, "+" for reverse.
	 *
	 * 2P[k] = Q[k]+Q@[N/2 - k] -/+ (i * exp(-/+2 * pi * i * k / N) * (Q[k]-Q@[N/2 - k]))
     *       = Q[k]+Q@[N/2 - k] -/+ (i * (cos(2*pi*j/N) -/+ (i * sin(2*pi*j/N))) * (Q[k]-Q@[N/2 - k]))
	 *	     = Q[k]+Q@[N/2 - k] -/+ ((+/-sin(2*pi*j/N) + (i * cos(2*pi*j/N))) * (Q[k]-Q@[N/2 - k]))
	 * 
	 * let s = sin(2*pi*k/N)
     *     c = cos(2*pi*k/N)
	 *     Z = Q[k] + Q@[N/2 - k]
	 *     Y = Q]k] - Q@[N/2 - k]
	 *  
     * Then 
     * 
     * 2P[k] = Z -/+ ((+/-s + c*i) * Y)
	 */
	endRow = startRow + rowsToProcess;
	for(size_t row=startRow; row<endRow; row++) {
		size_t startCol = 0;
		#if DEBUG_TWIST
		//size_t rowNmk = numRows - row;		/* row for Q[N/2 - k] */
		#endif
		
		/* k, N/2-k */
		size_t k = row;
		size_t nmk = nOver2 - k;
	
		/* 
		 * Starting cos, sin for each row
		 * rowSinK is the running, actual value of sin[k]; sinK gets
		 * multiplied by -1 for reverse FFT. 
		 *
		 * We have to init these to *something* to avoid a compiler error.
		 */
		FFTVector cosKV = {0};
		FFTVector rowSinKV = {0};
		FFTVector sinKV;
		
		/*
		 * For N/2-k, which is at the end of its row.
		 */
		FFTVector cosNmkV = {0};
		FFTVector rowSinNmkV = {0};
		FFTVector sinNmkV;
		
		size_t lastCol;
		if(row == (numRows >> 1)) {
			/* middle row, even number of elements, process half of it. */
			lastCol = (numCols >> 1) - 1;
		}
		else {
			/* normal case, all columns */
			lastCol = numCols - 1;
		}
		tprintf("row %u startCol %u lastCol %u\n",
			(unsigned)row, (unsigned)startCol, (unsigned)lastCol);
			
		for(size_t col=startCol; col<=lastCol; col+=FFT_FLOATS_PER_VECTOR) {
		
			/*
			 * Prefetch current Q[k] and Q[N/2-k]
			 */
			kIn.prefetch();
			nmkIn.prefetch();
			
			/* 
			 * sinK = +/- sin(2*pi*k/N)
			 * cosK =     cos(2*pi*k/N)
			 */
			 
			/* 
			 * Calculate actual sin/cos for column zero and then periodically 
			 * after that; infer new values in between "actual" updates. 
			 * FIXME this needs work to use optimized sin/cos lookup tables. 
			 */
			if( (col == 0) ||							/* always on column 0 */
			    ( (FFT_SIN_RECALC_REAL > 0) &&			/* re-sync enabled, and... */
				  ((col % FFT_SIN_RECALC_REAL) == 0)	/* it's time to resync */
				)
			   ) {
				FFTVectUnion cosKU;
				FFTVectUnion sinKU;
				FFTVectUnion cosNmkU;
				FFTVectUnion sinNmkU;
				for(unsigned dex=0; dex<FFT_FLOATS_PER_VECTOR; dex++) {
					size_t incr = dex * numRows;
					fftSinCos(mfftPlan, k+incr,   &cosKU.f[dex],   &sinKU.f[dex]);
					/* note f[0] in cosNmkU is the cos[nmk], and we decrement from there */
					fftSinCos(mfftPlan, nmk-incr, &cosNmkU.f[dex], &sinNmkU.f[dex]);
				}
				cosKV   = cosKU.v;
				cosNmkV = cosNmkU.v;
				sinKV   = rowSinKV   = sinKU.v;
				sinNmkV = rowSinNmkV = sinNmkU.v;
			}
			else {
				/* 
				 * Infer new sin/cos 
				 * FFTFloat f = (cosK * d) - (rowSinK * t);
				 * rowSinK = sinK = (rowSinK * d) + (cosK * t);
				 * cosK = f;
				 */
				FFTVector fV     = FFTVectSub(FFTVectMul(cosKV, dV),    FFTVectMul(rowSinKV, tV));
				rowSinKV = sinKV = FFTVectAdd(FFTVectMul(rowSinKV, dV), FFTVectMul(cosKV, tV));
				cosKV = fV;
				
				/* 
				 * f = (cosNmk * d) + (rowSinNmk * t);
				 * rowSinNmk = sinNmk = (rowSinNmk * d) - (cosNmk * t);
				 * cosNmk = f;
				 */
				fV                   = FFTVectAdd(FFTVectMul(cosNmkV, dV), FFTVectMul(rowSinNmkV, tV));
				rowSinNmkV = sinNmkV = FFTVectSub(FFTVectMul(rowSinNmkV, dV), FFTVectMul(cosNmkV, tV));
				cosNmkV = fV;
			}
			
			if(!forward) {
				sinKV   = FFTVectMul(sinKV, minus1V);
				sinNmkV = FFTVectMul(sinNmkV, minus1V);
			}
			
			#if 0
			tprintf("k %u row %u col %u sin[k] %f sin [k+1] %f cos[k] %f cos[k+1] %f\n",
				(unsigned)k, (unsigned)row, (unsigned)col, 
				VEC0(sinKV), VEC1(sinKV), VEC0(cosKV), VEC1(cosKV));
			#endif
			
			/*
			 * Snag current Q[k] and Q[N/2-k]
			 * Q[N/2-k] is in reverse order
			 */
			kIn.loadVect(qkRV, qkIV);
			nmkIn.loadVectRev(qnmkRV, qnmkIV);
			kIn.offset(FFT_FLOATS_PER_VECTOR);
			nmkIn.offset(-FFT_FLOATS_PER_VECTOR);
			
			/* 
			 * Z = Q[k] + Q@[N/2 - k]
			 * Y = Q]k] - Q@[N/2 - k]
			 */
			FFTVector ZrV, ZiV, YrV, YiV;
			ZrV = FFTVectAdd(qkRV, qnmkRV);
			ZiV = FFTVectSub(qkIV, qnmkIV);
			YrV = FFTVectSub(qkRV, qnmkRV);
			YiV = FFTVectAdd(qkIV, qnmkIV);
			#if 0
			tprintf("   Q[k] = [%f, %f]  Q[N-k] = [%f, %f]\n", 
				VEC0(qkRV), VEC0(qkIV), VEC0(qnmkRV), VEC0(qnmkIV));
			tprintf("   Z = [%f, %f];  Y = [%f, %f]\n", 
				VEC0(ZrV), VEC0(ZiV), VEC0(YrV), VEC0(YiV));
			#endif
			
			/* 
			 * prod = (sinK + cosK*i) * Y
			 *
			 * FFTFloat prodR, prodI;
			 * prodR = (sinK * Yr) - (cosK * Yi);
			 * prodI = (sinK * Yi) + (cosK * Yr);
			 */
			FFTVector prodRV, prodIV;
			prodRV = FFTVectSub(FFTVectMul(sinKV, YrV), FFTVectMul(cosKV, YiV));
			prodIV = FFTVectAdd(FFTVectMul(sinKV, YiV), FFTVectMul(cosKV, YrV));
			//tprintf("   (sinK + cosK*i) * Y = [%f, %f]\n", VEC0(prodRV), VEC0(prodIV));
			
			
			/* 
			 * 2P[k] = Z -/+ ((s + c*i) * Y)
			 *       = Z -/+ prod
			 */
			FFTVector outRV, outIV;
			if(forward) {
				/* 
				 * *kOutR = Zr - prodR;
				 * *kOutI = Zi - prodI;
				 */
				outRV = FFTVectSub(ZrV, prodRV);
				outIV = FFTVectSub(ZiV, prodIV);
			}
			else {
				/* 
				 * *kOutR = Zr + prodR;
				 * *kOutI = Zi + prodI;
				 */
				outRV = FFTVectAdd(ZrV, prodRV);
				outIV = FFTVectAdd(ZiV, prodIV);
			}
			//tprintf("   2P[k] = [%f, %f]\n", VEC0(outRV), VEC0(outIV));
			
			kOut.storeVect(outRV, outIV);
			kOut.offset(FFT_FLOATS_PER_VECTOR);
			
			/* 
			 * Repeat for P[N/2 - k]; 'k' and 'nmk' are swapped
			 * 
			 * Z = Q[k] + Q@[N/2 - k] 
			 * Y = Q]k] - Q@[N/2 - k]
			 */
			#if 0
			tprintf("k %u row %u col %u sin[k] %f sin [k+1] %f cos[k] %f cok[k+1] %f\n",
				(unsigned)nmk, (unsigned)rowNmk, (unsigned)(numCols - col - 1),
				VEC0(sinNmkV), VEC0(cosNmkV), VEC1(sinNmkV), VEC1(cosNmkV));
			#endif
			
			// Zr = qkR + qnmkR;		(we already have this)
			ZiV = FFTVectSub(qnmkIV, qkIV);
			YrV = FFTVectSub(qnmkRV, qkRV);
			#if 0
			// Yi = qkI + qnmkI;		(we already have this)
			tprintf("   Q[k] = [%f, %f]  Q[N-k] = [%f, %f]\n", 
				VEC0(qnmkRV), VEC0(qnmkIV), VEC0(qkRV), VEC0(qkIV));
			tprintf("   Z = [%f, %f];  Y = [%f, %f]\n", 
				VEC0(ZrV), VEC0(ZiV), VEC0(YrV), VEC0(YiV));
			#endif
			
			/* 
			 * prod = (sinK + cosK*i) * Y
			 */
			prodRV = FFTVectSub(FFTVectMul(sinNmkV, YrV), FFTVectMul(cosNmkV, YiV));
			prodIV = FFTVectAdd(FFTVectMul(sinNmkV, YiV), FFTVectMul(cosNmkV, YrV));
			//tprintf("   (sinK + cosK*i) * Y = [%f, %f]\n", VEC0(prodRV), VEC0(prodIV));			
			
			/* 
			 * 2P[k] = Z -/+ ((s + c*i) * Y)
			 *       = Z -/+ prod
			 */
			if(forward) {
				outRV = FFTVectSub(ZrV, prodRV);
				outIV = FFTVectSub(ZiV, prodIV);
			}
			else {
				outRV = FFTVectAdd(ZrV, prodRV);
				outIV = FFTVectAdd(ZiV, prodIV);
			}
			//tprintf("   2P[k] = [%f, %f]\n", VEC0(outRV), VEC0(outIV));
		
			/* 
			 * Everthing at this end has been in reverse order; swap back to 
			 * natural order here
			 */
			nmkOut.storeVectRev(outRV, outIV);
			nmkOut.offset(-FFT_FLOATS_PER_VECTOR);
			/* these increments of k get discarded at the end of a row */
			k   += (numRows * FFT_FLOATS_PER_VECTOR);
			nmk -= (numRows * FFT_FLOATS_PER_VECTOR);
		}	/* for col */
	}	/* for row */
	return MR_Success;
}

#else	/* !FFT_INTEL */

/* This is just for compiling/linking. It never runs. */
static MFFTReturn fftReal1DTwistOpt(
	MatrixFFTPlan			mfftPlan,
	const FFTComplex		*src,
	FFTComplex				*dst,
	size_t					numRows,
	size_t					numCols,		// in complex elements
	bool					forward,
	
	/* The rows we actually process. */
	size_t					startRow,
	size_t					rowsToProcess)
{
	return MR_Internal;
}

#endif	/* FFT_INTEL */

#pragma mark --- Public SPI ---

MFFTReturn fft1DRealTwist(
	MatrixFFTPlan		mfftPlan,
	const FFTComplex	*src,
	FFTComplex			*dst,
	bool				forward,
	size_t				numRows,
	size_t				numCols,		// in complex elements
	size_t				startRow,
	size_t				rowsToProcess)
{
	if((numRows > FFT_FLOATS_PER_VECTOR) && 
	   FFT_INTEL &&
	   !FFT_FWD_TWIST_REAL_SMALL) {
		return fftReal1DTwistOpt(mfftPlan, src, dst, numRows, numCols, forward, startRow, rowsToProcess);
	}
	else {
		return fftReal1DTwistScalar(mfftPlan, src, dst, numRows, numCols, forward, startRow, rowsToProcess);
	}
}


