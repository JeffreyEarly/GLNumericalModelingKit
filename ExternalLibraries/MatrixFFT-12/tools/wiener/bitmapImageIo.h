/*
 * bitmapImageIo.h - high-level support for reading and writing image files.
 *
 * Created 7/13/2007. 
 */
 
/* 
 * Linking against this module requires the following at load time:
 *
 * -framework CoreServices -framework ApplicationServices -framework CoreFoundation
 */
#ifndef	_BITMAP_IMAGE_H_
#define _BITMAP_IMAGE_H_

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * Input/output file type specification.
 */
typedef enum {
	BI_FT_Default = 0,			/* figure it out from file name */
	BI_FT_JPEG,					/* lossy JPEG */
	BI_FT_JPEG_Lossless,		/* lossless JPEG */
	BI_FT_TIFF,
	BI_FT_JPEG2000,				/* lossy JPEG 2000 */
	BI_FT_JPEG2000_Lossless,	/* lossless JPEG 2000 */
	BI_FT_Invalid				/* not used in interface */
} BIFileType;

/* 
 * Info about a bitmap image.
 * Note: color images are RGB, 4 bytes per pixel, no alpha (which takes up the last
 * byte of each 32 bit pixel, e.g. byte 3 of the pixel at bytes 0...3). 
 */
typedef struct {
	size_t imageWidth;			/* in pixels */
	size_t imageHeight;			/* in pixels */
	size_t bitsPerComponent;
	size_t bitsPerPixel;
	size_t bytesPerRow;			/* allows for a row using up more memory than the pixels in it occupy */
	size_t effectHeight;		/* includes extra mallocd data */
} BIImageInfo;

/* obtain a list of file type args available for use in BIFileTypeFromArg() */
extern const char *BIFileTypeArgString();

/*
 * Map a cmd line file type arg to BIFileType.
 * Returns BI_FT_Invalid if necessary.
 */
extern BIFileType BIFileTypeFromArg(
	const char			*arg);

/* 
 * Specification of padding for BIReadFile().
 */
typedef enum {
	PM_None,		/* no padding, mallocd array exactly fits image */
	PM_Pixels,		/* round up to specified pad size in pixels */
	PM_Bytes		/* round up to specified pad size in bytes */
} BIPadMode;

/* 
 * Read image file.
 * Returns nonzero on error.
 * Row size, and the total number of rows mallocd, are optionally rounded up to 
 * next padSize pixels or bytes per specified padMode. 
 */
extern int BIReadFile(
	const char			*fileName,
	BIFileType			fileType,
	BIPadMode			padMode,
	unsigned			padSize,
	BIImageInfo			*imageInfo,		/* RETURNED */
	unsigned char		**bitmap);		/* mallocd and RETURNED; caller must free */

/* 
 * Write image file.
 * Returns nonzero on error.
 */
extern int BIWriteFile(
	const char			*fileName,
	BIFileType			fileType,
	const BIImageInfo	*imageInfo,
	const unsigned char *bitmap);	

#ifdef __cplusplus
}
#endif

#endif	/* _BITMAP_IMAGE_H_ */

