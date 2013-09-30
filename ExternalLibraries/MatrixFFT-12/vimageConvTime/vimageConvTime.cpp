/*	File: vimageConvTime.cpp
	
	Description:
		Measure performance of vImage convolution.
	
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
 * Copyright (c) 2007 Apple, Inc. 
 * 
 * vimageConvTime.cpp - measure performance of vImage convolution
 *
 * Created Oct 1 2007. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/devRandom.h>
#include <Accelerate/Accelerate.h>

#define VT_IMAGE_ROWS_DEF	300		
#define VT_IMAGE_COLS_DEF	400

/* 
 * Kernel is a square of dimension k x k where k is odd.
 * If min != max we iterate for(k=min; k<=max; k+=incr).
 */
#define VT_KERNEL_SIZE_MIN	3	
#define VT_KERNEL_SIZE_MAX	17
#define VT_KERNEL_SIZE_INCR	2

/* calculate rowBytes given width in pixels, rounding up so rows are 16-byte aligned */
#define VT_ROW_BYTES(w)		(((w * sizeof(float)) + 0xf) & ~0xf)

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -r imageRows      -- image rows; default %u\n", VT_IMAGE_ROWS_DEF);
	printf("  -c imageCols      -- image columns; default %u\n", VT_IMAGE_COLS_DEF);
	printf("  -k kernelSizeMin  -- min kernel size, must be odd, default %u\n", 
		VT_KERNEL_SIZE_MIN);
	printf("  -K kernelSizeMax  -- max kernel size, must be odd, default %u\n", 
		VT_KERNEL_SIZE_MAX);
	printf("  -i kernelSizeIncr -- kernel size increment, must be even, default %u\n",
		VT_KERNEL_SIZE_INCR);
	printf("  -p                -- preallocate temp buffer\n");
	printf("  -e                -- edge extend; default is black background fill\n");
	printf("  -w                -- wall time; default is user time\n");
	printf("  -C                -- output in Crandall format; default is table\n");
	printf("  -n                -- no banner\n");
	exit(1);
}

/*
 * Parameters for one test run. All memory preallocated and written with 
 * test pattern before doTest() called. 
 */
typedef struct {
	/* inputs */
	void				*srcBuf;			/* source vImage_Buffer.data */
	void				*dstBuf;			/* destination vImage_Buffer.data */
	float				*kernel;
	vImagePixelCount	kernelSize;			/* 'k', as in k x k */
	vImagePixelCount	imageRows;	
	vImagePixelCount	imageCols;	
	void				*tempBuffer;		/* optional */
	bool				edgeExtend;			/* false --> kvImageBackgroundColorFill
											 * true  --> kvImageEdgeExtend */
	bool				wallTime;
	
	/* output */
	double				runTime;			/* in seconds */
} VTParams;


/* 
 * Core test routine; performs one convolution. All data buffers allocated by caller. 
 * Returns nonzero on error. 
 */
static int doTest(VTParams *vtp)
{
	vImage_Buffer src;
	vImage_Buffer dst;
	

	src.height = dst.height = vtp->imageRows;
	src.width  = dst.width  = vtp->imageCols;
	src.rowBytes = dst.rowBytes = VT_ROW_BYTES(src.width);
	src.data = vtp->srcBuf;
	dst.data = vtp->dstBuf;
	
	vImage_Flags flags = 0;
	if(vtp->edgeExtend) {
		flags = kvImageEdgeExtend;
	}
	else {
		flags = kvImageBackgroundColorFill;
	}

	double startTime = fftGetTime(vtp->wallTime);
	vImage_Error verr = vImageConvolve_PlanarF(&src, &dst,  vtp->tempBuffer,
		0, 0,							// region of interest
		vtp->kernel, vtp->kernelSize, vtp->kernelSize, 
		/* background fill 0.0 */
		0.0, flags);
	double endTime = fftGetTime(vtp->wallTime);
	if(verr) {
		printf("***vImageConvolve_PlanarF() returned %ld\n", (long)verr);
		return -1;
	}
	vtp->runTime = endTime - startTime;
	return 0;
}

/* 
 * Generate a string for exporting results to external tools for analysis.
 */
static unsigned crandCtr = 0;

static char *appendCrandallFormat(
	vImagePixelCount imageRows,
	vImagePixelCount imageCols,
	vImagePixelCount kernelSize,
	double runTime,
	char *crandallStr)				/* optional existing entry */
{
	/* nSqkSq = R * C * m^2 */
	unsigned long long rckSq = imageRows * imageCols;
	rckSq *= (unsigned long long)kernelSize;
	rckSq *= (unsigned long long)kernelSize;
	
	char newEntry[200];
	
	snprintf(newEntry, 200, "{%f, %llu}, ", runTime, rckSq);
	if((++crandCtr % 4) == 0) {
		strcat(newEntry, "\n");
	}
	
	unsigned newLen = strlen(newEntry) + 1;
	if(crandallStr) {
		newLen += strlen(crandallStr);
	}
	char *newStr = (char *)malloc(newLen);
	if(crandallStr) {
		strcpy(newStr, crandallStr);
		free(crandallStr);
	}
	else {
		newStr[0] = '\0';
	}
	strcat(newStr, newEntry);
	return newStr;
}

#define RRAND_BUFSIZE	2000

/* generate random float signal */
static void genRandFloats(
	float *buf,
	size_t numSamples)
{
	unsigned bufSize = RRAND_BUFSIZE;
	
	if(bufSize > numSamples) {
		bufSize = numSamples;
	}
	
	/* fill up our buffer */
	float *inBuf = (float *)malloc(bufSize * sizeof(float));
	float *inp = inBuf;
	for(size_t dex=0; dex<bufSize; dex++) {
		*inp++ = getRandomDouble();
	}

	/* copy to caller's buffer */
	float *outp = buf; 
	inp = inBuf;
	float *inEnd = inBuf + bufSize;
	
	for(unsigned dex=0; dex<numSamples; dex++) {
		*outp++ = *inp++;
		if(inp == inEnd) {
			inp = inBuf;
		}
	}
	free(inBuf);
}

int main(int argc, char **argv)
{
	/* user-spec'd variables */
	vImagePixelCount imageRows = VT_IMAGE_ROWS_DEF;
	vImagePixelCount imageCols = VT_IMAGE_COLS_DEF;
	vImagePixelCount kernelSizeMin = VT_KERNEL_SIZE_MIN;
	vImagePixelCount kernelSizeMax = VT_KERNEL_SIZE_MAX;
	unsigned kernelSizeIncr = VT_KERNEL_SIZE_INCR;
	bool preallocTempBuf = false;
	bool edgeExtend = false;
	bool wallTime = false;
	bool crandallFormat = false;
	bool printBanner = true;
	
	int arg;
	while ((arg = getopt(argc, argv, "r:c:k:K:i:pewCnh")) != -1) {
		switch (arg) {
			case 'r':
				imageRows = atoi(optarg);
				break;
			case 'c':
				imageCols = atoi(optarg);
				break;
			case 'k':
				kernelSizeMin = atoi(optarg);
				break;
			case 'K':
				kernelSizeMax = atoi(optarg);
				break;
			case 'i':
				kernelSizeIncr = atoi(optarg);
				break;
			case 'p':
				preallocTempBuf = true;
				break;
			case 'e':
				edgeExtend = true;
				break;
			case 'w':
				wallTime = true;
				break;
			case 'C':
				crandallFormat = true;
				break;
			case 'n':
				printBanner = false;
				break;
			case 'h':
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
	
	/* validate inputs */
	if(!(kernelSizeMin & 0x01) || !(kernelSizeMax & 0x01)) {
		printf("***kernelSize min and max (%lu, %lu) must be odd\n",
			(unsigned long)kernelSizeMin, (unsigned long)kernelSizeMax);
		exit(1);
	}
	if(kernelSizeIncr & 0x01) {
		printf("***kernelSizeIncr (%u) must be even\n", kernelSizeIncr);
		exit(1);
	}
	if(kernelSizeMax < kernelSizeMin) {
		printf("***kernelSizeMax must be >= kernelSizeMin\n");
		exit(1);
	}
	
	if(printBanner) {
		fftPrintTestBanner("Two-dimension Real Convolve", "vImage", false, "Random", NULL, 0, 0);
		printf("\n");
		printf("Image Rows | Image Cols | Kernel size | Convolution time (s)\n");
		printf("-----------+------------+-------------+---------------------\n");
	}
	
	/* allocate some memory */
	vImagePixelCount rowBytesMax = VT_ROW_BYTES(imageCols);
	size_t bufSize = rowBytesMax * imageRows;
	size_t imageSize = imageRows * imageCols;
	void *srcBuf = malloc(bufSize);
	void *dstBuf = malloc(bufSize);
	if((srcBuf == NULL) || (dstBuf == NULL)) {
		printf("***malloc failure (bufSize %lu\n", bufSize);
		exit(1);
	}
	size_t kernSize = kernelSizeMax * kernelSizeMax * sizeof(float);
	float *kernel = (float *)malloc(kernSize);
	if(kernel == NULL) {
		printf("***malloc failure (kernSize %lu\n", kernSize);
		exit(1);
	}	
	
	/* and the optional temp buffer */
	void *tempBuf = NULL;
	if(preallocTempBuf) {
		vImage_Buffer src;
		src.height = imageRows;
		src.width = imageCols;
		src.rowBytes = rowBytesMax;
		src.data = srcBuf;
		vImage_Flags flags = kvImageGetTempBufferSize;
		if(edgeExtend) {
			flags |= kvImageEdgeExtend;
		}
		else {
			flags |= kvImageBackgroundColorFill;
		}
		vImage_Error vrtn = vImageConvolve_PlanarF(&src, &src,  
				NULL,							// tempBuffer
				0, 0,							// region of interest
				kernel, kernelSizeMax, kernelSizeMax, 
				0.0,							// background fill 
				flags);
		if(vrtn <= 0) {
			printf("***vImageConvolve_PlanarF() returned %ld on get temp bufsize op\n",
				(long)vrtn);
			exit(1);
		}
		printf("...temp buffer size %ld\n", (long)vrtn);
		tempBuf = malloc((size_t)vrtn);
		if(tempBuf == NULL) {
			printf("***malloc failure (temp bufsize %ld\n", (long)vrtn);
			exit(1);
		}
	}
	
	genRandFloats((float *)srcBuf, imageSize);
	genRandFloats(kernel, kernelSizeMax);
	
	/* collect everything into VTParams */
	VTParams vtp;
	vtp.srcBuf     = srcBuf;
	vtp.dstBuf     = dstBuf;
	vtp.kernel     = kernel;
	vtp.imageRows  = imageRows;
	vtp.imageCols  = imageCols;
	vtp.tempBuffer = tempBuf;
	vtp.edgeExtend = edgeExtend;
	vtp.wallTime   = wallTime;
	
	/* here we go */
	vImagePixelCount kernelSize;
	char *crandallStr = NULL;

	if(crandallFormat) {
		crandallStr = strdup("\n/* {time(seconds), n^2 * m^2} */\n");
	}
	
	
	for(kernelSize=kernelSizeMin; kernelSize<=kernelSizeMax; kernelSize+=kernelSizeIncr) {
		vtp.kernelSize = kernelSize;
		int irtn = doTest(&vtp);
		if(irtn) {
			printf("***test failure; aborting\n");
			exit(1);
		}
		
		printf( "%8lu   | %8lu   | %8lu    | %10.5f\n", 
			(unsigned long)imageRows, (unsigned long)imageCols, 
			(unsigned long)kernelSize,
			vtp.runTime);
			
		if(crandallFormat) {
			/* save this result by appending it to the string we'll output when finished */
			crandallStr = appendCrandallFormat(imageRows, imageCols, kernelSize, vtp.runTime, 
				crandallStr);
		}
	}
	
	if(crandallFormat) {
		printf("\n%s\n", crandallStr);
	}
	
	/* clean up */
	COND_FREE(srcBuf);
	COND_FREE(dstBuf);
	COND_FREE(kernel);
	COND_FREE(tempBuf);
	return 0;
}
