/*
 * plplotUtils.h simplified interface to plplot.
 *
 * Created Dec. 19 2007. 
 */
 
#ifndef	_PLPLOT_UTILS_H_
#define _PLPLOT_UTILS_H_

/* 
 * plplot is available from http://plplot.sourceforge.net/index.html
 *
 * You need /usr/local/include/plplot in your header search list for this to build.
 */
#include <plplot.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int32_t minX;			/* x on left side of graph */
	int32_t maxX;			/* x on right side of graph */
	PLFLT minY;				/* y on bottom of graph */
	PLFLT maxY;				/* y at top of graph */
	const char *xAxisName;	/* optional */
	const char *yAxisName;	/* optional */
} PlplotSetup;

/* common point codes */
#define PLPOINT_NONE		0	/* e.g. when plotting lines */
#define PLPOINT_DOT			1	/* small dot */
#define PLPOINT_CIRCLE		4	/* small open circle */
#define PLPOINT_TRIANGLE	7
#define PLPOINT_TARGET		9	/* open circle with dot in center */
#define PLPOINT_LARGE_DOT	17	/* large dot */

/* our own enumerated colors */
#define PLCOLOR_BLACK		0
#define PLCOLOR_WHITE		1
#define PLCOLOR_RED			2
#define PLCOLOR_GREEN		3
#define PLCOLOR_BLUE		4
#define PLCOLOR_BROWN		5
#define PLCOLOR_CYAN		6
#define PLCOLOR_VIOLET		7
#define PLCOLOR_MAX			PLCOLOR_VIOLET

/*
 * Definition of one line - all share the same x array 
 */
typedef struct {
	PLFLT *fy;
	PLINT pointCode;		/* point glyph, only used if plPoints true */
	unsigned color;			/* PLCOLOR_BLACK, etc. */
} LineDef;
	
/* 
 * Plot specified x/y coordinates, x as int32_t and y as float, with an arbitrary number
 * of separate lines, as lines and/or points. 
 * Output is PDF file, location specified by 'fileName'.
 */
extern int plplotLines(
	PlplotSetup *setup,	
	uint32_t numSamples,	/* size of ix[] and fy[] arrays */
	uint32_t numLines,		/* size of lineDef[] array */
	int32_t *ix,			/* numSamples */
	LineDef *lineDef,		/* numLines */
	const char *graphName,
	const char *fileName,
	bool plotPoints,		/* true: plot points */
	bool plotLine,			/* true: plot line */
	bool skipTrailZeroes);	/* don't plot zero Y values at end of graph */

/* 
 * Plot a histogram showing the number of occurences of each possible 
 * value of 'samples'.
 * Output is PDF file, location specified by 'fileName'.
 */
extern int plplotHist(
	const int32_t *samples,
	uint32_t numSamples,
	const char *graphName,
	const char *fileName,
	const char *xAxisName,		/* optional */
	const char *yAxisName);		/* optional */
	
/* 
 * Plot a histogram of prebinned data. X values are int32_t's, and the corresponding
 * Y values - the counts for each X - are uint32_t's. 
 * Output is PDF file, location specified by 'fileName'.
 */
int plplotBins(
	uint32_t numBins,
	const int32_t *x,			/* numBins of X values, monotonically increasing */
	const uint32_t *y,			/* numBins of Y values for each associated X */
	const char *graphName,
	const char *fileName,
	const char *xAxisName,		/* optional */
	const char *yAxisName);		/* optional */

#ifdef __cplusplus
}
#endif

#endif	/* _PLPLOT_UTILS_H_ */

