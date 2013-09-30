/*	File: MatrixFFT.h 
	
	Description:
		Primary public API for MatrixFFT module.
	
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
 * MatrixFFT.h - Primary public API for MatrixFFT module.
 *			  
 * Created 12/16/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#ifndef	_MATRIX_FFT_H_
#define _MATRIX_FFT_H_

#include <libMatrixFFT/MatrixFFTConfig.h>
#include <libMatrixFFT/complexBuf.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/***
 *** Note: the precision of the floating point types, and the format of complex
 ***       data, used in this library are determined a compile time by the flags
 ***       in MatrixFFTConfig.h. Be sure you are aware of the configuration 
 ***       defined in that file when using this library.
 ***/
 
#pragma mark --- types and #defines ---

/* 
 * Return codes specific to this module.
 */
typedef enum {
	MR_Success = 0,
	MR_Unsupported = 1,			// unsupported configuration or function
	MR_VDSP = 2,				// vDSP error
	MR_Memory = 3,				// out of memory
	MR_IllegalArg = 4,			// illegal argument
	MR_Internal = 5				// internal error (in this library) 
} MFFTReturn;

/*
 * Describe signal format. Time domain data - the input to a forward FFT - is
 * in always in row order. Frequency domain data - the output from a
 * forward FFT - can be in one of the following three formats.
 */
typedef enum {
	MF_RowOrder = 0,
	MF_ColumnOrder = 1,
	MF_CustomOrder = 2
} MFFTFormat;

/* 
 * Flags passed to mfftCreatePlan(). 
 *
 * MCF_HintTranspose, when set, indicates that the output of the forward FFT
 *    will be transposed to normal row order. Currently used in 1-D complex
 *    and 1-D real-signal FFTs to optimize the internal matrix organization 
 *    for some configurations; when this flag is set in these configurations, 
 *    FFT performance proper decreases, but transpose performance increases 
 *    much more, so that overall throughput is greater. This flag has no effect 
 *    on actual functionality - only on performance. This flag is ignored in 
 *    2-D complex and real-signal FFTs. 
 */
#define MCF_HintTranspose		0x0001	

/* 
 * Flags passed to mfftExecute().
 */
#define MEF_TransposeOutput     0x0001		/* only valid for forward FFT */
#define MEF_TransposeInput		0x0002		/* only valid for inverse FFT */
#define MEF_NormOutput			0x0004		/* scale output of inverse FFT as 
                                             * appropriate to recover original data */

/* 
 * Analagous to vDSP's FFTSetup, this is our state for a given 
 * operation. It can only be used for the operation for which it
 * was explicitly created; any change to sizes or type require a 
 * new MatrixFFTPlan. 
 * At this interface we have an opaque pointer. 
 */
struct MatrixFFTPlanStruct;
typedef struct MatrixFFTPlanStruct *MatrixFFTPlan;

#pragma mark --- Public FFT API ---

/* 
 * Create a MatrixFFTPlan.
 * 
 * Returned MatrixFFTPlan must be freed via mfftFreePlan(). 
 * 
 * dims       : Dimesions; currently 1 and 2 are supported.
 * n          : Array of log2 of dimensions. For 2D, n[0] is log2(numRows), n[1]
 *              is logn2(numColumns).
 * real       : True if forward FFT is real-to-complex; false if complex-to-complex.
 * optFlags   : Option flags, MCF_* defined above. 
 * maxThreads : Max number of pthreads to use for host FFT processing. Specify 
 *			    0 to use the default (which is currently the number of CPU cores).
 */
extern MFFTReturn mfftCreatePlan(
	unsigned			dims,		
	unsigned			*n,			
	bool				real,
	uint32_t			optFlags,
	unsigned			numThreads,
	MatrixFFTPlan		*mfftPlan);		/* RETURNED */
	
/*
 * Free a MatrixFFTPlan.
 */
extern void mfftFreePlan(
	MatrixFFTPlan		mfftPlan);

/* 
 * Execute FFT for the specified MatrixFFTPlan.
 *
 * MEF_TRANSPOSE_xxx can be specified in the optFlags argument
 * to transpose data to/from the native format as required. 
 *
 * In the absence of MEF_TRANSPOSE_xxx flags, input is assumed to 
 * be in native format for input and will be in native
 * format on output. 
 *
 * If a transpose flag is set for an operation for which no transposition is 
 * required (e.g. the native input format to a forward FFT is in fact in 
 * row order), the transpose flag is ignored. 
 *
 * Note that in some configurations, the input buffer is overwritten
 * during the execution of this function.
 *
 * For real-signal FFTs, the real data is split into an FFTComplex in the 
 * same manner is with vDSP: the FFTComplex.real components take the even-index
 * real data, and the FFTComplex.imag components take the odd-index real data. 
 * When configured to use interleaved data (FFT_SPLIT_COMPLEX = 0 in
 * MatrifFFTConfig.h), real data will be in memory in normal sequential
 * order, and it's OK to (carefully) cast your FFTFloat* to a FFTComplex* for
 * this operation. 
 *
 * Buffers should be aligned to FFT_MEM_ALIGNMENT bytes for maximum performance. 
 *
 * The outputs of inverse FFTs have to be normalized in order to obtain the actual
 * time-domain data. Normalization involves dividing each element of the output
 * by some multiple of the total number of elements in the signal as described below. 
 * If you specify the MEF_NormOutput on the mfftExecute() call for the inverse FFT,
 * the library will perform the required normalization in a highly optimized manner. 
 * The MEF_NormOutput flag is ignored during forward FFTs. 
 *
 *
 * Notes on the data produced by mfftExecute()
 * -------------------------------------------
 *
 * When the output format of a forward FFT is in column order, that data can be transposed 
 * to normal row order during mfftExecute() via the MEF_TransposeOutput option flag, or 
 * separately using mfftTranspose(). Note that specifying e.g. the MEF_TransposeOutput
 * option flag on a forward FFT when the native output is in row order has no effect; 
 * there is no performance penalty in doing so. 
 *
 * For all real-signal FFTs, when the output of a forward FFT is described as being in "row 
 * order", the format of the data is the same as the format of the output of the corresponding
 * FFTs implemented in the vDSP library. These formats are described in the following
 * document:
 *      http://developer.apple.com/hardwaredrivers/ve/downloads/vDSP_Library.pdf
 *
 * See the sections entitled "Data Packing for Real FFTs" and "Scaling Fourier Transforms"
 * for details.
 *
 * Descriptions of the frequency-domain data (the output of the forward FFT) and the 
 * required normalization (for the output of the inverse FFT) follow. 
 * Note for 1-D FFTs, 'numRows' and 'numCols' refer to the internal representation of the 
 * data; these values can be obtained from mfftRectangle().
 *
 * 
 *   1-D complex
 *   ---------------------
 *   Forward output format : Row order, or column order with numRows and numCols 
 *                           swapped, depending on signal size and configuration.
 *   Inverse normalization : Requires normalization by 1/N where N is the total number of 
 *                           complex elements. Use the MEF_NormOutput flag to perform 
 *                           normalization during the mfftExecute() call.
 *
 *   2-D complex
 *   ---------------------
 *   Forward output format : Row order or column order, depending on signal size and 
 *                           configuration.  
 *   Inverse normalization : Requires normalization by 1/N where N is the total number of 
 *                           complex elements. Use the MEF_NormOutput flag to perform 
 *                           normalization during the mfftExecute() call.
 *
 *   1-D real signal
 *   ---------------------
 *   Forward output format : Row order, or column order with numRows and numCols 
 *                           swapped, depending on signal size and configuration.
 *   Inverse normalization : Requires normalization by 1/2N where N is the total number of 
 *                           real elements. Use the MEF_NormOutput flag to perform 
 *                           normalization during the mfftExecute() call.
 *
 *   2-D real signal
 *   ---------------------
 *   Forward output format : Row order or column order, depending on signal size and 
 *                           configuration. You can determine which format will obtain for 
 *                           a given plan via mfftNativeFormat(). Note that specifying e.g. 
 *                           the MEF_TransposeOutput option flag on a forward FFT when the 
 *                           native output is in row order has no effect; there is no 
 *                           performance penalty in doing so. 
 *   Inverse normalization : Requires normalization by 1/2N where N is the total number of 
 *                           real elements. Use the MEF_NormOutput flag to perform 
 *                           normalization during the mfftExecute() call.

 */
extern MFFTReturn mfftExecute(
	MatrixFFTPlan		mfftPlan,
	uint32_t			optFlags,
	bool				forward,		
	FFTComplex			*inBuf,		/* NOT const! */
	FFTComplex			*outBuf);

/*
 * Transpose a signal to/from normal row order.
 *
 * If 'input' is true, convert from row order to native input format, else
 * convert from native output format to row order. 
 *
 * In-place is supported - i.e. inBuf == outBuf - but may require
 * extra internal memory allocation to execute. 
 *
 * If the current algorithm does not require a transposition for 
 * the specified input/output, this routine either does nothing (if
 * inBuf == outBuf) or does a simple copy (from inBuf to outBuf).  
 */
extern MFFTReturn mfftTranspose(
	MatrixFFTPlan		mfftPlan,
	bool				forward,		
	bool				input, 
	const FFTComplex	*inBuf,
	FFTComplex			*outBuf);
	
#pragma mark --- Obtain MatrixFFTPlan attributes ---

/* 
 * Obtain optimal native format for the input to, or the output from,
 * the forward FFTs for the specified MatrixFFTPlan. 
 *
 * These are the formats required and produced by the low-level algorithms; 
 * signals can be converted between row order and column order or custom 
 * format either at mfftExecute() time, via the MEF_TRANSPOSE_xxx flags, or 
 * explicitly via mfftTranspose(). 
 */
extern MFFTReturn mfftNativeFormat(
	MatrixFFTPlan		mfftPlan,
	bool				input,
	MFFTFormat			*format);		/* RETURNED */

/*
 * Return the number of threads in use by a MatrixFFTPlan.
 * If mfftPlan is NULL, return the default number of threads used
 * (i.e. the result of passing a numThreads of 0 to mfftCreatePlan()).
 */
extern MFFTReturn mfftNumThreads(
	MatrixFFTPlan		mfftPlan,
	unsigned			*numThreads);	/* RETURNED */

/*
 * Obtain the number of rows and columns used in a MatrixFFTPlan.
 *
 * For plans configured for 2-D FFTs, numRows and numCols
 * correspond to the n[0] and n[1] values passed to mfftCreatePlan()
 * (though the values obtained here are true sizes, not log2 of the
 * sizes). 
 *
 * For plans configured for 1-D FFTs, numRows and numCols are the 
 * values used in representing the data during matrix decomposition 
 * operations. In such cases, (numRows * numCols) will always equal 
 * to 2^(n[0]) where n[0] is the log2 of the array size passed to 
 * mfftCreatePlan(). 
 *
 * The output of this function is undefined when number of dimensions
 * is larger than 2 (which is currently unimplemented).
 *
 * These values are of interest mainly when the application needs to 
 * access column-order frequency-domain data.
 */
extern MFFTReturn mfftRectangle(
	MatrixFFTPlan		mfftPlan,
	size_t				*numRows,		/* RETURNED */
	size_t				*numCols);		/* RETURNED */
	
#pragma mark --- Utility functions ---

/* 
 * Obtain C-string version of a MFFTReturn.
 */
extern const char *mfftErrStr(
	MFFTReturn			mrtn);

/* 
 * Dump error info to stdout (like perror()). 
 */
extern void mfftPrintErrInfo(
	const char			*func, 
	MFFTReturn			mrtn);

/*
 * For best performance, all memory passed to mfftExecute() and mfftTranspose
 * should be aligned to FFT_MEM_ALIGNMENT bytes.
 */
#define FFT_MEM_ALIGNMENT	64

#ifdef __cplusplus
}
#endif

#endif	/* _MATRIX_FFT_H_ */
