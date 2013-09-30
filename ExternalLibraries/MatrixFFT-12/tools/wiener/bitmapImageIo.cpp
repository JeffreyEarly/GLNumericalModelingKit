/*
 * bitmapImageIo.h - high-level support for reading and writing image files.
 *
 * Created 7/13/2007. 
 */

#include "bitmapImageIo.h"
#include <CoreServices/CoreServices.h>
#include <CoreFoundation/CoreFoundation.h>
#include <ApplicationServices/ApplicationServices.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>

#define CFRELEASE(cf)	if(cf != NULL) { CFRelease(cf); cf = NULL; }

/* obtain a list of file type args available for use in BIFileTypeFromArg() */
extern const char *BIFileTypeArgString()
{
	return "jpeg|jpeg-l|2000|2000-l|tiff";
}

/* Map a cmd line file type arg to BIFileType */
typedef struct {
	const char *arg;
	BIFileType fileType;
} BIFileTypeArgs;

static const BIFileTypeArgs fileTypeArgs[] = 
{
	{ "jpeg",		BI_FT_JPEG				},
	{ "jpeg-l",		BI_FT_JPEG_Lossless		},
	{ "tiff",		BI_FT_TIFF				},
	{ "2000",		BI_FT_JPEG2000			},
	{ "2000-l",		BI_FT_JPEG2000_Lossless }
};
#define NUM_FILE_TYPES	(sizeof(fileTypeArgs) / sizeof(fileTypeArgs[0]))

BIFileType BIFileTypeFromArg(
	const char *arg)
{
	const BIFileTypeArgs *ftap = fileTypeArgs;
	for(unsigned dex=0; dex<NUM_FILE_TYPES; dex++) {
		if(!strcmp(arg, ftap->arg)) {
			return ftap->fileType;
		}
		ftap++;
	}
	printf("***Invalid file type specification (%s)\n", arg);
	return BI_FT_Invalid;
}

/* Map of filename extensions to UTI string. */
typedef struct {
	CFStringRef exten;		/* e.g. CFSTR("tif") */
	const CFStringRef uti;	/* e.g. kUTTypeTIFF */
} UTIMap;

static const UTIMap utiMap[] = 
{
	{ CFSTR("tif"),		kUTTypeTIFF},
	{ CFSTR("tiff"),	kUTTypeTIFF},
	{ CFSTR("TIF"),		kUTTypeTIFF},
	{ CFSTR("TIFF"),	kUTTypeTIFF},
	{ CFSTR("jpg"),		kUTTypeJPEG},
	{ CFSTR("jpeg"),	kUTTypeJPEG},
	{ CFSTR("JPG"),		kUTTypeJPEG},
	{ CFSTR("JPEG"),	kUTTypeJPEG},
	{ CFSTR("jp2"),		kUTTypeJPEG2000},
	{ CFSTR("JP2"),		kUTTypeJPEG2000},
	{ CFSTR("png"),		kUTTypePNG},
	{ CFSTR("PNG"),		kUTTypePNG}
};
#define NUM_UTI_MAPS	(sizeof(utiMap) / sizeof(utiMap[0]))

/*
 * Map a filename and/or a BIFileType to a UTI string.
 */
static CFStringRef BIGetUTI(
	CFURLRef url,
	BIFileType fileType)
{
	switch(fileType) {
		case BI_FT_JPEG:
		case BI_FT_JPEG_Lossless:
			return kUTTypeJPEG;
		case BI_FT_TIFF:
			return kUTTypeTIFF;
		case BI_FT_JPEG2000:
		case BI_FT_JPEG2000_Lossless:
			return kUTTypeJPEG2000;
		case BI_FT_Default:
			/* figure it out from file extension */
			break;
		default:
			fprintf(stderr, "***BIGetUTI internal error\n");
			return NULL;
	}
	
	CFStringRef fileExt = CFURLCopyPathExtension(url);
	if(fileExt == NULL) {
		fprintf(stderr, "***BIGetUTI(BI_FT_Default) - no extension available\n");
		return NULL;
	}
	
	for(unsigned mapDex=0; mapDex<NUM_UTI_MAPS; mapDex++) {
		const UTIMap *map = &utiMap[mapDex];
		if(CFEqual(fileExt, map->exten)) {
			return map->uti;
		}
	}

	fprintf(stderr, "***BIGetUTI(BI_FT_Default) - unknown extension\n");
	return NULL;
}

#define ROUND_TO(num, quanta)	(((num + quanta - 1) / quanta) * quanta)

/* read image file */
int BIReadFile(
	const char			*fileName,
	BIFileType			fileType,
	BIPadMode			padMode,
	unsigned			padSize,
	BIImageInfo			*imageInfo,		/* RETURNED */
	unsigned char		**bitmap)		/* mallocd and RETURNED; caller must free */
{
	CFURLRef fileURL = CFURLCreateFromFileSystemRepresentation(kCFAllocatorDefault, 
		(unsigned char*)fileName, strlen(fileName), FALSE);
	if(fileURL == NULL) {
		fprintf(stderr, "***BIReadFile: Error on CFURLCreateFromFileSystemRepresentation\n");
		return -1;
	}

	CFStringRef keys[1] = {kCGImageSourceTypeIdentifierHint};
	CFStringRef values[1] = {BIGetUTI(fileURL, fileType)};
	
	if(values[0] == NULL) {
		return -1;
	}
	
	CFDictionaryRef		optionsDict = NULL;
	CGImageSourceRef	imageSourceRef = NULL;
	CGImageRef			imageRef = NULL;
	CGColorSpaceRef		rgbColorSpaceRef = NULL;
	CGContextRef		bitmapContextRef = NULL;
	CGBitmapInfo		bitmapInfo = 0;
	CGImageAlphaInfo	alpha;
	unsigned			bytesPerPixel = 4;
	
	optionsDict = CFDictionaryCreate( kCFAllocatorDefault, 
		(const void **)keys, (const void **)values, 1,  
		&kCFTypeDictionaryKeyCallBacks,  &kCFTypeDictionaryValueCallBacks );
	/* subsequent errors to errOut: */
	
	int ourRtn = 0;
	
	/* source file --> CGImageRef */
	imageSourceRef = CGImageSourceCreateWithURL(fileURL, optionsDict);
	if(imageSourceRef == NULL) {
		fprintf(stderr, "***BIReadFile: Error on CGImageSourceCreateWithURL\n");
		ourRtn = 1;
		goto errOut;
	}
	CFRELEASE(fileURL);
	imageRef = CGImageSourceCreateImageAtIndex(imageSourceRef, 0, optionsDict );
	if(imageRef == NULL) {
		fprintf(stderr, "***BIReadFile: Error on CGImageSourceCreateImageAtIndex\n");
		ourRtn = 1;
		goto errOut;
	}
	
	imageInfo->imageWidth			= CGImageGetWidth(imageRef);
	imageInfo->imageHeight			= CGImageGetHeight(imageRef);
	imageInfo->bitsPerComponent		= CGImageGetBitsPerComponent(imageRef);
	imageInfo->bitsPerPixel			= CGImageGetBitsPerPixel(imageRef);
	
	if(imageInfo->bitsPerPixel == 8) {
		/* the image is gray */
		rgbColorSpaceRef = CGColorSpaceCreateWithName(kCGColorSpaceGenericGray);
		imageInfo->bytesPerRow = imageInfo->imageWidth; 
		alpha = kCGImageAlphaNone;
		bitmapInfo = CGImageGetBitmapInfo(imageRef);
		bytesPerPixel = 1;
	}
	else {
        rgbColorSpaceRef = CGColorSpaceCreateDeviceRGB();
		imageInfo->bytesPerRow = imageInfo->imageWidth * 4;
		alpha = kCGImageAlphaPremultipliedLast;
		bitmapInfo = kCGBitmapByteOrder32Big | alpha;
	}
	if(rgbColorSpaceRef == NULL) {
		fprintf(stderr, "***BIReadFile: Error on CGColorSpaceCreateWithName\n");
		ourRtn = 1;
		goto errOut;
	}
	
	/* optionally pad */
	imageInfo->effectHeight = imageInfo->imageHeight;
	if(padMode != PM_None) {
		if(padSize == 0) {
			fprintf(stderr, "***Pad size of 0 invalid\n");
			ourRtn = 1;
			goto errOut;
		}
		unsigned padSizeBytes = padSize;
		if(padMode == PM_Pixels) {
			padSizeBytes *= bytesPerPixel;
		}
		imageInfo->bytesPerRow = ROUND_TO(imageInfo->bytesPerRow, padSizeBytes);
		/* also round up row count */
		imageInfo->effectHeight = ROUND_TO(imageInfo->imageHeight, padSize);
	}

	*bitmap = (unsigned char *)malloc(imageInfo->bytesPerRow * imageInfo->effectHeight);
	
	bitmapContextRef = CGBitmapContextCreate(*bitmap, 
		imageInfo->imageWidth, imageInfo->imageHeight, 
		imageInfo->bitsPerComponent, imageInfo->bytesPerRow, 
		rgbColorSpaceRef, 
		bitmapInfo);
	if(bitmapContextRef == NULL) {
		fprintf(stderr, "***BIReadFile: Error creating CGBitmapContext\n");
		ourRtn = 1;
		goto errOut;
	}

	/* enable high quality interpolation */
	CGContextSetInterpolationQuality(bitmapContextRef, kCGInterpolationHigh);

	/* Draw into the context */
	CGContextDrawImage(bitmapContextRef, CGRectMake(0, 0, imageInfo->imageWidth, imageInfo->imageHeight), 
		imageRef);

errOut:
	CFRELEASE(optionsDict);
	CFRELEASE(imageSourceRef);
	CFRELEASE(imageRef);
	CFRELEASE(rgbColorSpaceRef);
	CFRELEASE(bitmapContextRef);
	if(ourRtn) {
		if(*bitmap) {
			free(*bitmap);
			*bitmap = NULL;
		}
	}
	return ourRtn;
}

/* write image file */
int BIWriteFile(
	const char			*fileName,
	BIFileType			fileType,
	const BIImageInfo	*imageInfo,
	const unsigned char *bitmap)
{
	int						ourRtn = 0;
	CGImageAlphaInfo		alpha;
	CGBitmapInfo			bitmapInfo = 0;
	CGContextRef			bitmapContextRef = NULL;
	CGImageRef				imageRef = NULL;
	CFURLRef				fileURL = NULL;
	CGImageDestinationRef	imageDestRef = NULL;
	CGColorSpaceRef			rgbColorSpaceRef = NULL;
	CFStringRef				uti = NULL;
	CFDictionaryRef			propDict = NULL;
	
	if(imageInfo->bitsPerPixel == 8) {
		/* grayscale image */
		rgbColorSpaceRef = CGColorSpaceCreateWithName(kCGColorSpaceGenericGray);
		alpha = kCGImageAlphaNone;
		bitmapInfo = 0;
	}
	else {
		/* RGB no alpha */
        rgbColorSpaceRef = CGColorSpaceCreateDeviceRGB();
		alpha = kCGImageAlphaNoneSkipLast;
		bitmapInfo = kCGBitmapByteOrder32Big | alpha;
	}
	if(rgbColorSpaceRef == NULL) {
		fprintf(stderr, "***BIWriteFile: Error on CGColorSpaceCreateWithName\n");
		ourRtn = 1;
		goto errOut;
	}
	
	/* A bitmap-oriented CGContextRef based on bitmap data */
	bitmapContextRef = CGBitmapContextCreate((void *)bitmap, 
		imageInfo->imageWidth, imageInfo->imageHeight, 
		imageInfo->bitsPerComponent, imageInfo->bytesPerRow, 
		rgbColorSpaceRef, 
		bitmapInfo);
	if(bitmapContextRef == NULL) {
		fprintf(stderr, "***BIWriteFile: Error creating CGBitmapContext\n");
		ourRtn = 1;
		goto errOut;
	}
	
	/* CGContextRef --> CGImageRef */
    imageRef = CGBitmapContextCreateImage(bitmapContextRef);

	fileURL = CFURLCreateFromFileSystemRepresentation(kCFAllocatorDefault, 
		(unsigned char*)fileName, strlen(fileName), FALSE);
	if(fileURL == NULL) {
		fprintf(stderr, "***BIWriteFile: Error on CFURLCreateFromFileSystemRepresentation\n");
		ourRtn = 1;
		goto errOut;
	}

	uti = BIGetUTI(fileURL, fileType);
	if(uti == NULL) {
		ourRtn = 1;
		goto errOut;
	}
	imageDestRef = CGImageDestinationCreateWithURL(fileURL, 
		uti, 
		1, NULL );
	if(imageDestRef == NULL) {
		fprintf(stderr, "***BIWriteFile: Error on CGImageDestinationCreateWithURL\n");
		ourRtn = 1;
		goto errOut;
	}
	
	/* Some BIFileTypes require the specification of a "Lossless" property */
	switch(fileType) {
		case BI_FT_JPEG2000_Lossless:
		case BI_FT_JPEG_Lossless:
		{
			CFStringRef key = kCGImageDestinationLossyCompressionQuality;
			float valf = 1.0;
			CFNumberRef val = CFNumberCreate(NULL, kCFNumberFloatType, &valf);
			propDict = CFDictionaryCreate(NULL, (const void **)&key, (const void **)&val, 1, 
				&kCFTypeDictionaryKeyCallBacks, &kCFTypeDictionaryValueCallBacks);
			/* 
			 * The docs say we should be able to set these properties like this,
			 * but that has no effect. We need to pass it as the 'properties' argument
			 * to CGImageDestinationAddImage(), below. 
			 * See <rdar://problem/5670614> 
			 */
			//CGImageDestinationSetProperties(imageDestRef, propDict);
			CFRelease(val);
			break;
		}
		default:
			break;
	}
	
	CGImageDestinationAddImage(imageDestRef, imageRef, propDict);
	
	/* Write the image to disk */
	if(!CGImageDestinationFinalize(imageDestRef)) {
		fprintf(stderr, "***BIWriteFile: Error on CGImageDestinationFinalize\n");
		ourRtn = 1;
	}
	if(propDict != NULL) {
		CFRelease(propDict);
	}
errOut:
	CFRELEASE(bitmapContextRef);
	CFRELEASE(imageRef);
	CFRELEASE(fileURL);
	CFRELEASE(imageDestRef);
	CFRELEASE(rgbColorSpaceRef);
	return ourRtn;
}

