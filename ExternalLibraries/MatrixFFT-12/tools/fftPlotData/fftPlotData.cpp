/*
 * compareFFT.cpp - convert the output of compareFFT to either a
 *                  a plplot-generated graph or TeX source. 
 *
 * Created 7/25/2008. 
 * Copyright 2008 by Apple, Inc. 
 */
 
/* 
 * The files we're parsing have one of two forms. 
 * 1-D outputs look like this:
 *
 * ----start of file ----
 * FFT type    : One-dimension complex
 * Precision   : Double
 * Signal type : Random
 * Threads     : 8
 * System      : Version 10.6 (Build 10A209)
 * CPU         : Intel(R) Xeon(R) CPU X5482 @ 3.20GHz
 * Built as    : 32 bit
 * Complex     : Interleaved
 * Test time   : 10:14:46  Jan 9 2009
 * 
 *   Samples |Matrix CTGs|vecLib CTGs| FFTW CTGs 
 *  ---------+-----------+-----------+-----------
 *     2^13  |   2.129   |   3.993   |   0.692   
 *     2^14  |   2.633   |   4.222   |   1.364   
 * --- etc... ---
 *
 * 2-D input files look like this:
 *
 * ----start of file ----
 * FFT type    : Two-dimension complex
 * Precision   : Double
 * Signal type : Random
 * Threads     : 8
 * System      : Version 10.6 (Build 10A209)
 * CPU         : Intel(R) Xeon(R) CPU X5482 @ 3.20GHz
 * Built as    : 32 bit
 * Complex     : Interleaved
 * Test time   : 10:27:22  Jan 9 2009
 * 
 *   Width | Height | Samples |Matrix CTGs|vecLib CTGs| FFTW CTGs 
 *  -------+--------+---------+-----------+-----------+-----------
 *   2^6   |  2^6   |   2^12  |    1.521  |    0.184  |    0.418  
 *   2^6   |  2^7   |   2^13  |    2.163  |    1.383  |    0.819  
 * --- etc... ---
 *
 *
 * Notes:
 * -----
 * -- X axis label is fixed, X_AXIS_NAME ("Samples")
 * -- Y axis is fixed, Y_AXIS_NAME ("Cooley-Tukey GFlops")
 * -- Each column - e.g. Matrix, vecLib, FFTW - is a separate line on the graph, 
 *    separate colors and point glyphs.
 * -- Title is optionally specified on the command line, or else gleaned from "FFT Type" 
 *    and "Precision" lines
 * -- X axis labeled in log2(samples); edit this with a PDF editing tool if that's 
 *    not good enough
 * -- 2-D input files can have duplicate entries for some sizes, e.g. 2048x8192 followed by
 *    8192x2048. When we encounter a dup, we skip the second one. 
 * -- Input files normally have 3 sets of data. This can be overridden by the -n numGraphs
 *    option. 
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <libMatrixFFT/src/TextParser.h>		/* private interface */
#include <libMatrixFFT/fftUtils.h>		
#include "plplotUtils.h"

#define X_AXIS_NAME "lg N"
#define Y_AXIS_NAME "CTGflop/sec"

/* 
 * 5 lines is the max we can handle.
 * Normally, with 3 FFT implementations:
 * MatrixFFT = Green, circle
 * Accelerate = Blue, triangle
 * FFTW = Red, large dot
 */
static const PLINT pointCodes[] = {PLPOINT_CIRCLE, PLPOINT_TRIANGLE, PLPOINT_LARGE_DOT, 
								   PLPOINT_TARGET, PLPOINT_CIRCLE };
static const unsigned colors[]  = {PLCOLOR_GREEN, PLCOLOR_BLUE, PLCOLOR_RED, 
								   PLCOLOR_BLACK, PLCOLOR_VIOLET};
#define NUM_POINTS_COLORS	(sizeof(colors) / sizeof(colors[0]))

#define FFT_PLOT_DEBUG		0
#if		FFT_PLOT_DEBUG
#define dprintf(args...)			printf(args)
#else
#define dprintf(args...)	
#endif

static void usage(char **argv)
{
	printf("usage: %s [options] inFile outFile\n", argv[0]);
	printf("Options:\n");
	printf("   -l label       -- graph label; default is inferred from 'FFT Type' line\n");
	printf("   -t lineLabel   -- TeX mode line label (default is graph mode)\n");
	printf("   -c             -- TeX/convolve mode (default is FFT graph mode)\n");
	printf("   -n numGraphs   -- number of graph lines; default is 3\n");
	printf("   -y minY        -- minimum value on Y axis; default is 0.0\n");
	printf("   -Y maxY        -- maximum value on Y axis; default is from data\n");
	printf("   -r             -- raw data mode, not compareFFT\n");
	/* etc. */
	exit(1);
}

/* ASCII float to float */
static float strToFloat(
	const char *str)
{
	return strtof(str, NULL);
}

/* ASCII integer, a power of 2, in either decimal or 2^n form, to log2(n) */
static unsigned strToLog2(
	const char *str)
{
	char *carat = strchr(str, '^');
	if(carat) {
		long pwr = strtol(carat+1, NULL, 0);
		return (unsigned)pwr;
	}
	
	size_t n = strtoll(str, NULL, 0);
	unsigned log2n;
	
	if(!fftIsPowerOfTwo(n, &log2n)) {
		/* 
		 * This program has to build 32 bit because plplot builds 32 bit...handle
		 * one special case here that's troublesome.
		 */
		if(!strcmp(str, "4294967296")) {
			return 32;
		}
		printf("***Found a sample count that is not a power of two (%llu)\n", (unsigned long long)n);
		return 0;
	}
	return log2n;
}

/* 
 * Emit (to stdout) a line in a TeX table for FFT measurements. 
 *
 * The TeX for a 1-D table entry looks like this:
 *
 *           Complex & $2^{19}$ & 0.81 & 1.36 & 1.58  \\
 * 
 * And for a 2-D looks like this:
 *
 *           Complex & $2^{14}\times 2^{11}$ &  0.08 &  0.30 & 2.03 \\
 */

static void emitTexFFT(
	const char *lineLabel,		/* e.g. "Complex" */
	bool twoDim,				/* inferred by caller from banner */
	unsigned numTokens,			/* line as parsed by TextParser */
	const char **tokens,
	unsigned sourceCols,		/* CTGs in source */
	unsigned destCols)			/* CTGs in output */
{
	unsigned ctgToken;
	
	/* print line label and size */
	if(twoDim) {
		printf("%s & $2^{%u}\\times 2^{%u}$ ", 
			lineLabel, strToLog2(tokens[0]), strToLog2(tokens[2]));
		ctgToken = 6;
	}
	else {
		printf("%s & $2^{%u}$ ", lineLabel, strToLog2(tokens[0]));
		ctgToken = 2;
	}
	
	/* print first two CTGs (or less) */
	unsigned n = (destCols < 2) ? destCols : 2;
	for(unsigned dex=0; dex<n; dex++) {
		float val = strToFloat(tokens[ctgToken]);
		if(val == 0.0) {
			printf("& --- ");
		}
		else {
			printf("& %.2f ", val);
		}
		ctgToken += 2;
	}
	if(destCols == 2) {
		printf("\\\\\n");
		return;
	}
	
	if((sourceCols == 4) && (destCols == 3)) {
		/* here's the weird case where we skip column 3, rowOP */
		ctgToken += 2;
	}
	
	/* now print the rest */
	for(unsigned dex=2; dex<destCols; dex++) {
		float val = strToFloat(tokens[ctgToken]);
		if(val == 0.0) {
			printf("& --- ");
		}
		else {
			printf("& %.2f ", val);
		}
		ctgToken += 2;
	}
	printf("\\\\\n");
}
	
/*
 * Emit (to stdout) a line in a TeX table for Convolution measurements. 
 *
 * The TeX for a convolution table entry looks like this:
 *
 *           $64 \times 64$ &  $5 \times 5$ &  0.00013 & 0.00019 \\
 *
 * In our input files, rows and height of the image are specified separately but there
 * is only one dimension for the kernel. Also note that in the input file the image
 * is described as rows x columns, the reverse of what we emit. 
 */
static void emitTexConv(
	unsigned numTokens,			/* line as parsed by TextParser */
	const char **tokens)
{
	if(numTokens != 9) {
		printf("***Badly formatted convolution table\n");
		return;
	}
	printf("$%s \\times %s$ & $%s \\times %s$ & %s & %s\\\\\n",
		tokens[2], tokens[0], tokens[4], tokens[4], tokens[6], tokens[8]);
}
	
int main(int argc, char **argv)
{
	const char *label = NULL;
	char genLabel[LINE_LENGTH_MAX];
	bool texMode = false;
	const char *texLineLabel = NULL;
	bool convolveMode = false;
	float minYSpec = 0.0;
	float maxYSpec = 0.0;
	bool minYSpecd = false;
	bool maxYSpecd = false;
	bool rawMode = false;		// raw data output, not compareFFT
	
	/*
	 * The number of lines to actually graph
	 */
	unsigned numGraphLines = 3;
	
	int arg;
	extern char *optarg;
	while ((arg = getopt(argc, argv, "l:t:cn:y:Y:rh")) != -1) {
		switch (arg) {
			case 'l':
				label = optarg;
				break;
			case 't':
				texLineLabel = optarg;
				texMode = true;
				break;
			case 'n':
				numGraphLines = atoi(optarg);
				if((numGraphLines > NUM_POINTS_COLORS) && (numGraphLines != 8)) {
					printf("We can handle a max of %u graph lines\n", (unsigned)NUM_POINTS_COLORS);
					exit(1);
				}
				break;
			case 'c':
				convolveMode = true;
				texMode = true;
				break;
			case 'y':
				minYSpec = strtof(optarg, NULL);
				minYSpecd = true;
				break;
			case 'Y':
				maxYSpec = strtof(optarg, NULL);
				maxYSpecd = true;
				break;
			case 'r':
				rawMode = true;
				numGraphLines = 1;
				break;
			default:
			case 'h':
				usage(argv);
		}
	}

	int fileArgs = texMode ? 1 : 2;
	if(optind != (argc - fileArgs)) {
		printf("Missing filename(s)\n");
		usage(argv);
	}
	
	const char *inFile = argv[argc - fileArgs];
	const char *outFile = NULL;
	if(!texMode) {
		outFile = argv[argc - fileArgs + 1];
	}
	
	/* open TextParser from input file - exists on error */
	TextParser parser(inFile);

	/* Get FFT type */
	char lineBuf[LINE_LENGTH_MAX];
	if(!parser.findLine("FFT type", lineBuf)) {
		printf("***%s does not appear to be a compareFFT output file (1). Aborting.\n", inFile);
		exit(1);
	}

	/* 1D/2D from FFT type */
	bool twoDim = strstr(lineBuf, "Two") ? true : false;
	
	/* index of "samples" token inferred from 1D/2D */
	unsigned samplesDex = twoDim ? 4 : 0;
	
	/* label from FFT Type and Precision if caller doesn't specify */
	if(!texMode && (label == NULL)) {
		const char *sigType = strstr(lineBuf, "complex") ? "complex" : "real";

		if(!parser.findLine("Precision", lineBuf)) {
			printf("***%s does not appear to be a compareFFT output file (2.5). Aborting.\n", inFile);
			exit(1);
		}
		const char *precStr = strstr(lineBuf, "Single") ? "single" : "double";
		
		/*
		 * Examples:
		 * "1-D complex, double-precision FFT"
		 * "2-D real, single-precision FFT"
		 */
		sprintf(genLabel, "%s %s, %s-precision FFT", 
			twoDim ? "2-D" : "1-D",
			sigType, precStr);
		label = genLabel;
	}
	
	/* There should only be one of these lines, right before the start of the actual data */
	if(!parser.findLine("---------+", lineBuf)) {
		printf("***%s does not appear to be a compareFFT output file (3). Aborting.\n", inFile);
		exit(1);
	}
	
	/* cursor is at first data line, save this to restart here in the main loop */
	size_t firstDataLine = parser.getCursor();
	
	/* 
	 * Get this data line, count how many CTGs columns there are.
	 * For each CTG column there are two tokens ("| num").
	 */
	unsigned numTokens;
	const char **tokens;
	if(parser.getTokens(numTokens, tokens)) {
		printf("***%s does not appear to be a compareFFT output file (4). Aborting.\n", inFile);
		exit(1);
	}
	parser.freeTokens(numTokens, tokens);
	if(numTokens < (samplesDex + 1)) {
		printf("***%s does not appear to be a compareFFT output file (5). Aborting.\n", inFile);
		exit(1);
	}
	
	/* remember this, ensure each line is the same */
	unsigned tokensPerLine = numTokens;
	
	unsigned remCols = numTokens - samplesDex - 1;
	if(remCols & 0x01) {
		/* odd number remaining, we're lost */
		printf("***%s does not appear to be a compareFFT output file (6). Aborting.\n", inFile);
		exit(1);
	}
	
	/* 
	 * the number of actual CTGs in the table - normally 3; user can override
	 */
	uint32_t numDataColumns = remCols / 2;
	if(!convolveMode) {
		uint32_t numImpl = rawMode ? 1 : numDataColumns;
		if(numImpl < numGraphLines) {
			printf("***%s has data from %u FFTs; we can't graph %u lines. Aborting.\n",
				inFile, numDataColumns, numGraphLines);
				exit(1);
		}
		printf("...processing data from %u FFT implementations\n", numImpl);
	}
	
	/* 
	 * Prepare to process data for numGraphLines. We just realloc along 
	 * the way since we don't know how many lines are in the table. 
	 * This is all a nop for TeX mode. 
	 */
	uint32_t numSamples = 0;
	int32_t *ix = NULL;
	LineDef lineDefs[numGraphLines];
	
	if(!texMode) {
		if(numGraphLines == 8) {
			/* 
			 * Special case, comparing all 4 thread results from FFTW and MatrixFFT
			 */
			for(unsigned dex=0; dex<4; dex++) {
				/* FFTW */
				LineDef *lineDef = &lineDefs[dex];
				lineDef->fy = NULL;
				lineDef->pointCode = pointCodes[dex];
				lineDef->color = PLCOLOR_BLUE;
			}
			for(unsigned dex=4; dex<8; dex++) {
				/* MatrixFFT */
				LineDef *lineDef = &lineDefs[dex];
				lineDef->fy = NULL;
				lineDef->pointCode = pointCodes[dex - 4];
				lineDef->color = PLCOLOR_RED;
			}
		}
		else {
			for(unsigned dex=0; dex<numGraphLines; dex++) {
				LineDef *lineDef = &lineDefs[dex];
				lineDef->fy = NULL;
				lineDef->pointCode = pointCodes[dex];
				lineDef->color = colors[dex];
			}
		}
	}
	
	/* Rewind to start of first data line and proceed with main loop */
	parser.setCursor(firstDataLine);
	
	PLFLT minY = 100000.0;
	PLFLT maxY = 0.0;
	
	/***
	 *** main loop, one iteration per data line (i.e. per N) 
	 ***/
	for(;;) {
		if(parser.getTokens(numTokens, tokens)) {
			/* EOF, done */
			break;
		}
		if(numTokens != tokensPerLine) {
			/* trash between end of data and EOF, ignore */
			break;
		}
		
		if(texMode) {
			if(convolveMode) {
				emitTexConv(numTokens, tokens);
			}
			else {
				emitTexFFT(texLineLabel, twoDim, numTokens, tokens,
					numDataColumns, numGraphLines);
			}
			parser.freeTokens(numTokens, tokens);
			continue;
		}
		
		/* ix := log2(samples) */
		unsigned log2Samples = strToLog2(tokens[samplesDex]);
		
		/* detect 2-D dups */
		if(twoDim && (numSamples > 0) && (log2Samples == ix[numSamples - 1])) {
			dprintf("...skipping 2-D dup %s\n", tokens[samplesDex]);
			continue;
		}
		
		/* record new 'x' value */
		numSamples++;
		ix = (int32_t *)realloc(ix, numSamples * sizeof(int32_t));
		ix[numSamples-1] = log2Samples;
		dprintf("x[%u] = %u\n", (unsigned)numSamples-1, log2Samples);
		
		/* each graph line has two tokens - "| ctgs" - here */
		for(unsigned dex=0; dex<numDataColumns; dex++) {
		
			/* 
			 * lineDex is the lineDefs index 
			 * dex is index into parsed text
			 */
			unsigned lineDex = dex;			
			unsigned ctgDex;
			if(rawMode) {
				ctgDex = tokensPerLine - 1;
			}
			else {
				ctgDex = samplesDex + (2 * (dex+1));
			}
			LineDef *lineDef = &lineDefs[lineDex];
			lineDef->fy = (PLFLT *)realloc(lineDef->fy, numSamples * sizeof(PLFLT));
			
			float ctg = strToFloat(tokens[ctgDex]);
			dprintf("y[%u][%u] = %f\n", dex, (unsigned)numSamples-1, ctg);

			/* detect possible new min/max Y; record new Y for this graph */
			if(ctg > maxY) {
				maxY = ctg;
				dprintf("maxY = %f\n", maxY);
			}
			if(ctg < minY) {
				minY = ctg;
				dprintf("minY = %f\n", minY);
			}
			lineDef->fy[numSamples - 1] = ctg;
		}
		parser.freeTokens(numTokens, tokens);
	}
	
	if(texMode) {
		/* all done */
		return 0;
	}
	
	if(numSamples == 0) {
		printf("***No data found in file %s. Aborting.\n", inFile);
		exit(1);
	}
	
	/* OK, we have all the data, write to outFile */
	PlplotSetup setup;
	setup.minX = ix[0];		// may want some extra room here...
	setup.maxX = ix[numSamples-1];
	if(minYSpecd) {
		setup.minY = minYSpec;
	}
	else {
		setup.minY = 0.0;		// was minY, but we always want 0
	}
	if(maxYSpecd) {
		setup.maxY = maxYSpec;
	}
	else {
		setup.maxY = maxY;
	}
	setup.xAxisName = X_AXIS_NAME;
	setup.yAxisName = Y_AXIS_NAME;
	
	if(plplotLines(&setup, numSamples, numGraphLines, ix, lineDefs, label, outFile, 
		true,		// plotPoints
		true,		// plotLine
		true)) {	// skipTrailZeroes - avoid plotting NULL CTGs resulting from untested sizes
		printf("***Error writing to %s. Aborting.\n", outFile);
		return 1;
	}
	else {
		printf("...wrote %lu samples from %lu FFTs to %s\n",
			(unsigned long)numSamples, (unsigned long)numGraphLines, outFile);
	}
	return 0;
}
