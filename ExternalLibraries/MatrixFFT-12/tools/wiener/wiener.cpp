/*
 * Copyright (c) 2008 Apple Computer, Inc. All Rights Reserved.
 * 
 * wiener.cpp - Tool to apply blur and deblur via Wiener algorithm.
 * THIS IS NOT COMPLETE - it still needs the code to generate kernels and 
 * adjust the FFT'd kernel for the Wiener case. 
 *
 * Created Jan. 1 2008. 
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include "bitmapImageIo.h"
#include <libMatrixFFT/MatrixFFT.h>
#include <libMatrixFFT/fftConvolve.h>
#include <libMatrixFFT/complexBufUtils.h>
#include <libMatrixFFT/fftUtils.h>
#include <libMatrixFFT/src/PolyComplex.h>       /* private API */

#define RADIUS_DEF  3.0
#define C_DEF       1.0     /* ?? */
#define PAD_DEF     5

static void usage(char **argv)
{
	printf("usage: %s [options]\n", argv[0]);
	printf("Options:\n");
	printf("  -i inFile            -- image file input, required\n");
	printf("  -o outFile           -- image file output, required\n");
    printf("  -b                   -- blur; default is deblur\n");
    printf("  -r radius            -- blur radius; default %.1f\n", RADIUS_DEF);
    printf("  -c cVal              -- default %.2f\n", C_DEF);
    printf("  -p padding           -- default %u\nn", PAD_DEF);

	exit(1);
}

/*
 * Copy one channel of an RGB bitmap into a float array, with zero padding.
 */
static void copyImageToFloat(
    const unsigned char *rgb,
    const BIImageInfo *imageInfo,
    unsigned offset,            /* R, G, B --> 0, 1, 2 */
    unsigned padding,
    FFTFloat *f)
{
    size_t x, y;
    size_t floatsPerRow = imageInfo->imageWidth + (2 * padding);
    
    /* top pad */
    FFTFloat *fOut = f;
    for(y=0; y<padding; y++) {
        for(x=0; x<floatsPerRow; x++) {
            *fOut++ = 0.0;
        }
    }
    
    for(y=0; y<imageInfo->imageHeight; y++) {

        /* left pad */
        fOut = f + ((y + padding) * floatsPerRow);
        for(x=0; x<padding; x++) {
            *fOut++ = 0.0;
        }
        
        /* image row */
        const unsigned char *rgbIn = rgb + (y * imageInfo->bytesPerRow) + offset;
        for(x=0; x<imageInfo->imageWidth; x++) {
            *fOut++ = (FFTFloat)(*rgbIn);
            rgbIn += 4;
        }
        
        /* right pad */
        for(x=0; x<padding; x++) {
            *fOut++ = 0.0;
        }
    }
    
    /* bottom pad */
    for(y=0; y<padding; y++) {
        for(x=0; x<floatsPerRow; x++) {
            *fOut++ = 0.0;
        }
    }
}

/*
 * Copy an image from a FFTComplex to one RGB channel.
 */
static void copyImageFromComplex(
    const FFTComplex *cIn,
    unsigned fftN,              /* cIn is 2^fftN elements on each side */
    const BIImageInfo *imageInfo,
    unsigned offset,            /* R, G, B --> 0, 1, 2 */
    unsigned padding,
    unsigned char *rgb)
{
    size_t complexPerRow = 1 << (fftN - 1);
    size_t rowIn = padding;
    for(size_t y=0; y<imageInfo->imageHeight; y++) {
        PolyComplex pcIn((FFTComplex *)cIn, (rowIn * complexPerRow) + padding);
        unsigned char *rgbOut = rgb + (y * imageInfo->bytesPerRow) + offset;
        for(size_t x=0; x<imageInfo->imageWidth; x+=2) {
            FFTFloat f = pcIn.real();
            *rgbOut = (unsigned char)(f + 0.5);
            rgbOut += 4;
            if(x == (imageInfo->imageWidth - 1)) {
                /* odd size; discard imaginary */
                break;
            }
            f = pcIn.imag();
            *rgbOut = (unsigned char)(f + 0.5);
            rgbOut += 4;
            ++pcIn;
        }
        rowIn++;
    }
}

static int processChannel(
    unsigned char *rgb,             /* IN/OUT */
    const BIImageInfo *imageInfo,
    unsigned offset,                /* R, G, B --> 0, 1, 2 */
    unsigned padding,
    unsigned fftN,                  /* FFTs are 2^fftN elements on each side */
    MatrixFFTPlan mfftPlan,
    FFTFloat *imageFloat,
	FFTComplex *imageComp,          /* FFT size, allocated by caller */
	const FFTComplex *kernelFFT,    /* FFT(kernel) */
	FFTComplex *fftResult,          /* FFT size, allocated by caller */
    FFTFloat cVal)
{
    /* RGB channel --> float array --> complex array */
    copyImageToFloat(rgb, imageInfo, offset, padding, imageFloat);
    fftConvCopyImage(imageFloat, imageComp, 
        imageInfo->imageHeight + padding, imageInfo->imageWidth + padding,
        fftN, fftN,
        true);      // doZero
        
    /* FFT image, output in native order */
    MFFTReturn mrtn;
    mrtn = mfftExecute(mfftPlan, 0, true, imageComp, imageComp);
    if(mrtn) {
        mfftPrintErrInfo("mfftExecute", mrtn);
        return -1;
    }
    
    /* dyadic mul */
    MFFTFormat format;
    mrtn = mfftNativeFormat(mfftPlan, false, &format);
    if(mrtn) {
        mfftPrintErrInfo("mfftNativeFormat", mrtn);
        return -1;
    }
    mrtn = fftConvDyadicMul(imageComp, kernelFFT, fftResult, fftN, fftN, format);
    if(mrtn) {
        mfftPrintErrInfo("fftConvDyadicMul", mrtn);
        return -1;
    }
    
    /* Inverse FFT */
    mrtn = mfftExecute(mfftPlan, MEF_NormOutput, false, fftResult, fftResult);
    if(mrtn) {
        mfftPrintErrInfo("mfftExecute", mrtn);
        return -1;
    }
    
	/* 
	 * Scale. The output of one forward 2-D real FFT contains values
	 * that are 2x the actual FFT. When executing a reverse FFT 
	 * this 2x is handled via MEF_NormOutput, but we've done 2 forward 
	 * FFTs and multiplied the results together, followed by only one 
	 * reverse FFT, so we have to divide by two once more. 
     * The size argument is in complex elements. 
	 */
    size_t fftSize = 1 << (2 * fftN);
	fftScaleComplex(fftResult, 0.5, fftSize / 2);
    
    /* complex array --> RGB */
    copyImageFromComplex(fftResult, fftN, imageInfo, offset, padding, rgb);
    
    return 0;
}
    
int main(int argc, char **argv)
{
    FFTFloat radius = RADIUS_DEF;
    FFTFloat cVal = C_DEF;
    const char *inFile = NULL;
    const char *outFile = NULL;
    bool doDeblur = true;
    unsigned padding = PAD_DEF;
    
	int arg;
	while ((arg = getopt(argc, argv, "i:o:br:c:p:h")) != -1) {
		switch (arg) {
			case 'i':
				inFile = optarg;
				break;
            case 'o':
                outFile = optarg;
                break;
            case 'b':
                doDeblur = false;
                break;
            case 'r':
                radius = strtof(optarg, NULL);
                break;
            case 'c':
                cVal = strtof(optarg, NULL);
                break;
			case 'p':
				padding = atoi(optarg);
				break;
			default:
			case 'h':
				usage(argv);
		}
	}
	if(optind != argc) {
		usage(argv);
	}
	if(inFile == NULL) {
        printf("***no infile specified.\n");
        usage(argv);
    }
	if(outFile == NULL) {
        printf("***no infile specified.\n");
        usage(argv);
    }
    
    BIImageInfo imageInfo;
    unsigned char *rgb = NULL;
    if(BIReadFile(inFile, BI_FT_Default, 
            PM_None, 0,     // no padding
            &imageInfo, &rgb)) {
        printf("***Error reading %s. Aborting.\n", inFile);
        exit(1);
    }
        
    unsigned roundRadius = (unsigned)(radius + 0.5);
    unsigned kernelSize = (6 * roundRadius) + 1;
    
    /*
     * Calculate FFT size. It's a square, next power of 2 larger than
     * (image size + kernel size + padding - 1).
     */
    size_t imageSide = imageInfo.imageWidth;
    if(imageSide < imageInfo.imageHeight) {
        imageSide = imageInfo.imageHeight;
    }   
    unsigned fftN;            /* log2(fftRowsCols) */
    size_t fftRowsCols = fftRoundNumSamples(imageSide + kernelSize + padding - 1, &fftN);
    size_t imageSamples = (imageInfo.imageWidth + padding) * (imageInfo.imageHeight + padding);
    size_t totalFftSamples = fftRowsCols * fftRowsCols;
    
	/*
	 * Alloc and init FFTComplexes and FFTFloat arrays. 
	 *  -- one FFTFloat for image
	 *  -- one FFTFloat for kernel
	 *  -- one FFTComplex for split/padded image
	 *  -- one FFTComplex for split/padded kernel
	 *  -- one FFTComplex for fftConvolve() output 
	 */
	FFTFloat *imageFloat = NULL;
	FFTFloat *kernelFloat = NULL;
	/* aligned buffers */
	FFTComplex *imageComp  = NULL;
	FFTComplex *kernelComp = NULL;
	FFTComplex *fftResult = NULL;
	/* to-be-freed buffers */
	FFTComplex *imageFree  = NULL;
	FFTComplex *kernelFree = NULL;
	FFTComplex *resultFree = NULL;
	
	/* Raw image and kernel */
	imageFloat  = (FFTFloat *)malloc(sizeof(FFTFloat) * imageSamples);
	kernelFloat = (FFTFloat *)malloc(sizeof(FFTFloat) * kernelSize * kernelSize);
	
	/* complex buffers */
	size_t numComplex = totalFftSamples / 2;
	imageComp  = fftAllocComplexArrayAlign(numComplex, FFT_MEM_ALIGNMENT, &imageFree);
	kernelComp = fftAllocComplexArrayAlign(numComplex, FFT_MEM_ALIGNMENT, &kernelFree);
	fftResult  = fftAllocComplexArrayAlign(numComplex, FFT_MEM_ALIGNMENT, &resultFree);
	
	if(!imageFloat || !kernelFloat ||  !imageComp  || !kernelComp  || !fftResult) {
		printf("***Malloc failure for numFftSamples = %lu\n", (unsigned long)totalFftSamples);
		exit(1);
	}
	
	genConstComplex(fftResult, numComplex, 0.0);

    MatrixFFTPlan mfftPlan = NULL;
    MFFTReturn mrtn;
    unsigned n[2] = {fftN, fftN};
    mrtn = mfftCreatePlan(2, n, true, 0, 0, &mfftPlan);
    if(mrtn) {
        mfftPrintErrInfo("mfftCreatePlan", mrtn);
        exit(1);
    }
    
    /* TBD: generate kernel, FFT it, process it if doDeblur */
    
    /* now grind thru the channels, using the same FFT'd kernel for each one */
    if(processChannel(rgb, &imageInfo, 
            0, padding, fftN, mfftPlan,
            imageFloat, imageComp, kernelComp, fftResult, cVal)) {
        exit(1);
    }
    if(processChannel(rgb, &imageInfo,
            1, padding, fftN, mfftPlan,
            imageFloat, imageComp, kernelComp, fftResult, cVal)) {
        exit(1);
    }
    if(processChannel(rgb, &imageInfo,
            2, padding, fftN, mfftPlan,
            imageFloat, imageComp, kernelComp, fftResult, cVal)) {
        exit(1);
    }
    
    /* write result */
    if(BIWriteFile(outFile, BI_FT_Default, &imageInfo, rgb)) {
        printf("***error writing result to %s\n", outFile);
        exit(1);
    }
    else {
        printf("...wrote result to %s\n", outFile);
    }
	return 0;
}
