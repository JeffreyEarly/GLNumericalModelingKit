/*
 * plplotUtils.cpp simplified interface to plplot.
 *
 * Created Dec. 19 2007. 
 */

#include "plplotUtils.h"
#include <unistd.h>

/* map our API's colors to RGB */
typedef struct {
	PLINT	r;
	PLINT	g;
	PLINT	b;
} RGB;

#define PLOT_WHITE	255,255,255
#define PLOT_BLACK	0,0,0

static const RGB plColors[] = 
{
	{ PLOT_BLACK },		/* PLCOLOR_BLACK */
	{ PLOT_WHITE },		/* PLCOLOR_WHITE */
	{ 255,0,0 },		/* PLCOLOR_RED */
	{ 0,128,0 },		/* PLCOLOR_GREEN */
	{ 0,0,255 },		/* PLCOLOR_BLUE */
	{ 255,255,0 },		/* PLCOLOR_BROWN ? */
	{ 0,255,255 },		/* PLCOLOR_CYAN ? */
	{ 255,0,255 },		/* PLCOLOR_VIOLET ? */
};

#define PSTOPDF		"/usr/bin/pstopdf"

/* 
 * Convert plplot's output file, which is in PostScript format, 
 * to PDF.
 */
static int psToPdf(
	const char *inFile,
	const char *outFile)
{
	char cmd[300];
	snprintf(cmd, 300, "%s %s -o %s", PSTOPDF, inFile, outFile);
	int rtn = system(cmd);
	if(rtn) {
		fprintf(stderr, "error return %d from '%s'\n", rtn, cmd);
		return -1;
	}
	return 0;
}

#define PL_TMP_FILE	"/tmp/plplot"

/* Cook up temp file name for ps output */
static void makeTempFile(
	char *fn,
	int fnLen)
{
	snprintf(fn, fnLen, "%s_%d.ps", PL_TMP_FILE, (int)getpid());
}

/* basic "plot one line" atom */
static void plotOneLine(
	uint32_t numSamples,
	PLFLT *fx,
	LineDef *lineDef,
	bool plotPoints,		/* true: plot points */
	bool plotLine)			/* true: plot line */
{
	/* set drawing color - set the value for color map 0[1], then set it */
	const RGB *rgb = &plColors[lineDef->color];
	plscol0(1, rgb->r, rgb->g, rgb->b);			
	plcol0(1);
	
	if(plotPoints) {
		plpoin(numSamples, fx, lineDef->fy, lineDef->pointCode);
	}
	if(plotLine) {
		plline(numSamples, fx, lineDef->fy);
	}
}

#define TEMP_FN_LEN		128

/* 
 * Plot specified x/y coordinates, x as int32_t and y as float, with an arbitrary number
 * of separate lines, as lines and/or points. 
 * Output is PDF file, location specified by 'fileName'.
 */
int plplotLines(
	PlplotSetup *setup,	
	uint32_t numSamples,	/* size of ix[] and fy[] arrays */
	uint32_t numLines,		/* size of lineDef[] array */
	int32_t *ix,			/* numSamples */
	LineDef *lineDef,		/* numLines */
	const char *graphName,
	const char *fileName,
	bool plotPoints,		/* true: plot points */
	bool plotLine,			/* true: plot line */
	bool skipTrailZeroes)	/* don't plot zero Y values at end of graph */
{
	if(setup == NULL) {
		printf("***plplotLine: setup required\n");
		return -1;
	}
	
	char tmpFile[TEMP_FN_LEN];
	makeTempFile(tmpFile, TEMP_FN_LEN);
	
	PLFLT fx[numSamples];
	PLFLT minX = setup->minX;
	PLFLT maxX = setup->maxX;
	PLFLT minY = setup->minY;
	PLFLT maxY = setup->maxY;
	
	for(uint32_t dex=0; dex<numSamples; dex++) {
		fx[dex] = ix[dex];
	}
	
	const char *yName = "";
	if(setup->yAxisName) {
		yName = setup->yAxisName;
	}
	const char *xName = "";
	if(setup->xAxisName) {
		xName = setup->xAxisName;
	}
	
	plsdev ("psc");
	
	/* standard: background white, foreground (axes, labels, etc.) black */
	plscolbg(PLOT_WHITE);	
	plscol0(1, PLOT_BLACK);			

	plsfnam(tmpFile);
	plsdiori(1.0);						// portrait 
	plinit();
	plenv(minX, maxX, minY, maxY, 0, 0);
	pllab(xName, yName, graphName);
	
	for(uint32_t dex=0; dex<numLines; dex++) {
		uint32_t thisSamples = numSamples;
		if(skipTrailZeroes) {
			while((lineDef[dex].fy[thisSamples-1] == 0.0) && (thisSamples > 0)) {
				thisSamples--;
			}
			if(thisSamples == 0) {
				printf("***plplotLines: Warning: line with all zeroes skipped\n");
				continue;
			}
		}
		plotOneLine(thisSamples, fx,  &lineDef[dex], plotPoints, plotLine);
	}
	
	plend();
	
	int ourRtn = psToPdf(tmpFile, fileName);
	unlink(tmpFile);
	return ourRtn;
}

#define HIST_COLOR		0

/* 
 * Plot a histogram showing the number of occurences of each possible 
 * value of 'samples'.
 */
int plplotHist(
	const int32_t *samples,
	uint32_t numSamples,
	const char *graphName,
	const char *fileName,
	const char *xAxisName,		/* optional */
	const char *yAxisName)		/* optional */
{
	char tmpFile[TEMP_FN_LEN];
	makeTempFile(tmpFile, TEMP_FN_LEN);

	/* First determine the range, i.e. the number of bins */
	int32_t minSamp = samples[0];
	int32_t maxSamp = samples[0];
	for(uint32_t dex=0; dex<numSamples; dex++) {
		int32_t s = samples[dex];
		if(s < minSamp) {
			minSamp = s;
		}
		if(s > maxSamp) {
			maxSamp = s;
		}
	}

	/* When we specify PL_BIN_CENTRED, the min and max values are half the normal width */
	minSamp--;
	maxSamp++;
	
	PLINT numBins = maxSamp - minSamp + 1;
	
	/* One array containing the sample values, x */
	PLFLT *x = (PLFLT *)malloc(numBins * sizeof(PLFLT));
	int32_t binNum = minSamp;
	for(uint32_t dex=0; dex<(uint32_t)numBins; dex++) {
		x[dex] = binNum++;
	}
	
	/* Now make and fill the bins proper */
	PLFLT *y = (PLFLT *)malloc(numBins * sizeof(PLFLT));
	for(uint32_t dex=0; dex<(uint32_t)numBins; dex++) {
		y[dex] = 0;
	}
	PLFLT maxY = 0.0;

	for(uint32_t dex=0; dex<numSamples; dex++) {
		int32_t s = samples[dex];
		PLFLT *yp = y + s - minSamp;
		*yp += 1.0;
		if(*yp > maxY) {
			maxY = *yp;
		} 
	}
	
	const char *yName = yAxisName ? yAxisName : "";
	const char *xName = xAxisName ? xAxisName : "";
	
	#if	HIST_COLOR
	plsdev ("psc");
	plscolor(1);
	plscolbg(255, 255, 255);			/* white background */
	#else
	plsdev ("ps");
	#endif
	plsfnam(tmpFile);
	plsdiori(1.0);						// portrait 
	plinit();
	plenv(minSamp, maxSamp, 0, maxY, 0, 0);
	
	#if	HIST_COLOR
	/* can we alter colors of lines and the spaces inside the histograms? */
	plscolbg(255, 0, 0);				/* red background */
	plscol0(1, 255, 0, 0);				/* red foreground - no effect */
	#endif
	
	pllab(xName, yName, graphName);
	
	plbin(numBins, x, y, PL_BIN_CENTRED);
	plend();
	
	free(x);
	free(y);

	int ourRtn = psToPdf(tmpFile, fileName);
	unlink(tmpFile);
	return ourRtn;
}

/* 
 * Plot a histogram of prebinned data. X values are int32_t's, and the corresponding
 * Y values - the counts for each X - are uint32_t's. 
 */
int plplotBins(
	uint32_t numBins,
	const int32_t *x,			/* numBins of X values, monotonically increasing */
	const uint32_t *y,			/* numBins of Y values for each associated X */
	const char *graphName,
	const char *fileName,
	const char *xAxisName,		/* optional */
	const char *yAxisName)		/* optional */
{
	char tmpFile[TEMP_FN_LEN];
	makeTempFile(tmpFile, TEMP_FN_LEN);
	
	PLINT totalBins = numBins + 2;
	
	/* PLFLT array of sample values */
	PLFLT *xf = (PLFLT *)malloc(totalBins * sizeof(PLFLT));
	
	/* these two will have Y values of zero */
	xf[0] = x[0] - 1;		
	xf[totalBins - 1] = x[numBins-1] + 1;
	
	const int32_t *ip = x;
	PLFLT *op = xf + 1;
	for(uint32_t dex=0; dex<numBins; dex++) {
		*op++ = *ip++;
	}
	
	/* PLFLT array of bins */
	PLFLT *yf = (PLFLT *)malloc(totalBins * sizeof(PLFLT));
	yf[0] = 0.0;
	yf[totalBins - 1]  = 0.0;
	const uint32_t *uip = y;
	op = yf + 1;
	for(uint32_t dex=0; dex<numBins; dex++) {
		*op++ = *uip++;
	}
	
	/* get max Y value */
	uint32_t maxY = 0;
	uip = y;
	for(uint32_t dex=0; dex<numBins; dex++) {
		uint32_t currY = *uip++;
		if(currY > maxY) {
			maxY = currY;
		}
	}
	
	const char *yName = yAxisName ? yAxisName : "";
	const char *xName = xAxisName ? xAxisName : "";
	
	#if	HIST_COLOR
	plsdev ("psc");
	plscolbg(255, 255, 255);			/* white background */
	#else
	plsdev ("ps");
	#endif
	plsfnam(tmpFile);
	plsdiori(1.0);						// portrait 
	plinit();
	plenv(xf[0], xf[totalBins - 1], 0, maxY, 0, 0);
	
	#if	HIST_COLOR
	/* can we alter colors of lines and the spaces inside the histograms? */
	plscolbg(255, 0, 0);				/* red background */
	plscol0(1, 255, 0, 0);				/* red foreground - no effect */
	#endif
	
	pllab(xName, yName, graphName);
	
	plbin(totalBins, xf, yf, PL_BIN_CENTRED);
	plend();
	
	free(xf);
	free(yf);
	
	int ourRtn = psToPdf(tmpFile, fileName);
	unlink(tmpFile);
	return ourRtn;
}

