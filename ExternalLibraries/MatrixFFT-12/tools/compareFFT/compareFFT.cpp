/*
 * compareFFT.cpp - collate the outputs of a group of runs of mtime*.
 *
 * The outputs must match in "FFT Type", precision, signal type, and also to some extent 
 * the number of samples reported in each line (processing terminates on detection of a 
 * mismatch). 
 *
 * Created 8/16/2007. 
 * Copyright 2008 by Apple, Inc. 
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <libMatrixFFT/src/TextParser.h>		/* private interface */
#include <libMatrixFFT/fftUtils.h>		

static void usage(char **argv)
{
	printf("usage: %s [options] inFile1 inFile2 [...]\n", argv[0]);
	printf("Options:\n");
    printf("   -n     -- NUFFT mode - compare reference vs. optimized\n");
	/* etc. */
	exit(1);
}

/* 
 * Search for line containing lineTitle, return portion of the line following ":".
 * Returns nonzero if no "lineTitle" line found.
 * Line must contain at least 3 tokens - lineTitle, ":", and the string we return. 
 */
static int fineLineWithTitle(TextParser &parser, 
	const char *lineTytle,
	char *str)			// caller allocates, RETURNED
{
	char lineBuf[LINE_LENGTH_MAX];
	parser.setCursor(0);
	if(!parser.findLine(lineTytle, lineBuf)) {
		return 1;
	}
	
	unsigned numTokens;
	const char **tokens;
	numTokens = parser.parseLine(lineBuf, tokens);
	if(numTokens < 3) {
		return 1;
	}
	
	/* copy all tokens after ":" to caller's string */
	str[0] = '\0';
	for(unsigned dex=2; dex<numTokens; dex++) {
		strcat(str, tokens[dex]);
		if(dex < (numTokens - 1)) {
			strcat(str, " ");
		}
	}
	
	/* Convenience: "Accelerate" is too big. */
	if(strstr(str, "Accelerate")) {
		strcpy(str, "vecLib");
	}
	/* so is CLMatrix */
	else if(strstr(str, "CLMatrix")) {
		strcpy(str, "CLM");
	}
	/* ditto MatrixFFT */
	else if(strstr(str, "MatrixFFT")) {
		strcpy(str, "Matrix");
	}
	return 0;
}

/*
 * Determine if lines from two TextParsers containing *lineSpec are identical.
 * Returns nonzero if *lineSpec line not found in both files.
 */
static int compareLines(
	TextParser &parser1,
	TextParser &parser2,
	const char *lineSpec,			// e.g., "FFT Type" 
	bool &match,					// RETURNED, true if both parsers have the line and they match
	char *lineBufOut)				// RETURNED, one of the lines
{
	char lineBuf[LINE_LENGTH_MAX];

	parser1.setCursor(0);
	parser2.setCursor(0);
	if(!parser1.findLine(lineSpec, lineBuf)) {
		return 1;
	}
	if(!parser2.findLine(lineSpec, lineBufOut)) {
		return 1;
	}
	if(strcmp(lineBuf, lineBufOut)) {
		match = false;
	}
	else {
		match = true;
	}
	return 0;
}

/* 
 * Create a string with specified total length, containing "str1 str2" centered
 * within that string. Caller must free() the result.
 */
static char *centerStr(
	const char *str1, 
	const char *str2,
	unsigned totalLen)
{
	unsigned actLen = strlen(str1) + strlen(str2) + 1;
	unsigned leadSpace = (totalLen - actLen + 1) / 2;
	
	if(actLen > totalLen) {
		return strdup("BAD TOKEN");
	}
	char *outStr = (char *)malloc(totalLen + 1);
	memset(outStr, 0, totalLen + 1);
	for(unsigned dex=0; dex<leadSpace; dex++) {
		strcat(outStr, " ");
	}
	strcat(outStr, str1);
	strcat(outStr, " ");
	strcat(outStr, str2);
	
	if(actLen == totalLen) {
		return outStr;
	}
	
	unsigned padLen = totalLen - (actLen + 1);
	char *endp = outStr + actLen + leadSpace;
	for(unsigned dex=0; dex<padLen; dex++) {
		*endp++ = ' ';
	}
	outStr[totalLen] = '\0';
	return outStr;
}

int main(int argc, char **argv)
{
    bool nufftMode = false;
    
	int arg;
	while ((arg = getopt(argc, argv, "nh")) != -1) {
		switch (arg) {
            case 'n':
                nufftMode = true;
                break;
			default:
			case 'h':
				usage(argv);
		}
	}

	unsigned numFiles = argc - optind;
	if(numFiles < 2) {
		usage(argv);
	}
	
    if(nufftMode && (numFiles != 2)) {
        printf("***NUFFT mode requires two input files.\n");
        exit(1);
    }
    
	/* open TextParsers for each input file */
	char *fileName[numFiles];
	TextParser *parser[numFiles];
	for(unsigned dex=0; dex<numFiles; dex++) {
		fileName[dex] = argv[optind + dex];
		parser[dex] = new TextParser(fileName[dex]);
	}
	
	/*
	 * The header of all outputs looks like this:
	 *
	 * Library     : MatrixFFT
	 * FFT type    : Two-dimension complex
	 * Precision   : Single
	 * Signal type : Chirp
	 * Loops       : 10
	 * Threads     : 8
	 * Options     : Don't time first loop, In-place
	 * System      : Version 10.6 (Build 10A110)
	 * CPU         : Intel(R) Xeon(R) CPU X5482 @ 3.20GHz
	 * Built as    : 64 bit
	 * Complex     : Interleaved
	 * Test time   : 14:02:29  Oct 6 2008
	 *
	 * "Loops", "Threads", and "Options" lines are optional. 
	 * Do basic validation by getting "Library Type" for each file.
	 */
	char libType[numFiles][LINE_LENGTH_MAX];
	for(unsigned dex=0; dex<numFiles; dex++) {
		if(fineLineWithTitle(*parser[dex], "Library", libType[dex])) {
			printf("***%s is not an FFT output file\n", fileName[dex]);
			exit(1);
		}
	}
	
	/* 
	 * Make sure precisions, FFT types, and signal types match; 
	 * log them as well as log test time for the first file.
	 */
	char lineBuf[LINE_LENGTH_MAX];
	bool match;
	for(unsigned dex=1; dex<numFiles; dex++) {
		if(compareLines(*parser[0], *parser[dex], "FFT type", match, lineBuf)) {
			printf("***No FFT type line found\n");
			exit(1);
		}
		if(!match) {
			printf("***%s and %s measure different FFT types\n", fileName[0], fileName[dex]);
			exit(1);
		}
	}
	printf("%s\n", lineBuf);
	
	bool isConvolve = false;
	if(strstr(lineBuf, "Convolve")) {
		isConvolve = true;
	}
	
	for(unsigned dex=1; dex<numFiles; dex++) {
		if(compareLines(*parser[0], *parser[dex], "Precision", match, lineBuf)) {
			printf("***No Precision line found\n");
			exit(1);
		}
		if(!match) {
			printf("***%s and %s measure different Precision\n", fileName[0], fileName[dex]);
			exit(1);
		}
	}
	printf("%s\n", lineBuf);

	for(unsigned dex=1; dex<numFiles; dex++) {
		if(compareLines(*parser[0], *parser[dex], "Signal type", match, lineBuf)) {
			printf("***No Signal type line found\n");
			exit(1);
		}
		if(!match) {
			printf("***%s and %s measure different Signal types\n", fileName[0], fileName[dex]);
			exit(1);
		}
	}
	printf("%s\n", lineBuf);

	/* Optional "Threads" line - print it */
	if(parser[0]->findLine("Threads", lineBuf)) {
		printf("%s\n", lineBuf);
	}
	
	/* and rewind this one in any case; we may be past the Test time line */
	parser[0]->setCursor(0);
	
	/* 
	 * print following required fields:
	 *    system version (S/W)
	 *    CPU
	 *    Built as (32/64 bit)
	 *    Complex format
	 *    Test time 
	 *
	 * NOTE: the Complex Format line anly applies to MatrixFFT; we generally
	 * reuse data from one run of FFTW and vDSP for both Split and Interleaved
	 * format MFFT builds. So, the MFFT file really needs to be specified
	 * first to get the correct COmplex format in the collated output. 
	 */
	if(!parser[0]->findLine("System", lineBuf)) {
		printf("***No System found\n");
		exit(1);
	}
	printf("%s\n", lineBuf);
	if(!parser[0]->findLine("CPU", lineBuf)) {
		printf("***No CPU found\n");
		exit(1);
	}
	printf("%s\n", lineBuf);
	if(!parser[0]->findLine("Built as", lineBuf)) {
		printf("***No Built as found\n");
		exit(1);
	}
	printf("%s\n", lineBuf);
	if(!parser[0]->findLine("Complex", lineBuf)) {
		printf("***No Complex found\n");
		exit(1);
	}
	printf("%s\n", lineBuf);
	if(!parser[0]->findLine("Test time", lineBuf)) {
		printf("***No Test Time found\n");
		exit(1);
	}
	printf("%s\n\n", lineBuf);
	
	/* 
	 * The output we're parsing looks like this 2-D:
	 *
     *     Width  | Height |   Samples  | Wall time (s) |  CTGs 
     * -----------+--------+------------+---------------+---------
     *       1024 |   1024 |    1048576 |   0.606012    | 3.46058
	 *		   ...
	 *
	 * And for 1-D:
	 *
	 *     Samples | Wall time (s) |   CTGs 
	 *  -----------+---------------+--------
     *        1024 |   0.000013    | 3.86415
	 *
	 * For convolve:
	 *
	 *  Image Rows | Image Cols | Kernel size | Convolution time (s)
     *  -----------+------------+-------------+---------------------
     *
     * For NUFFT:
     *
     *   Samples | Wall time (s) | ms/loop 
     *  ---------+---------------+----------
     *      2^5  |     0.002     |  0.16
     * 
     *      1000   |     1000   |        3    |    0.07004
	 *
	 * So, we generally use the column header line to figure out what kind of output we're 
	 * dealing with.
	 */
	if(!isConvolve && !nufftMode) {
		if(!parser[0]->findLine("Samples", lineBuf)) {
			printf("***Malformed input file (no Samples line)\n");
			exit(1);
		}
	}
	
	unsigned totalTokens = 0;		/* expected tokens per line */
	unsigned ctgToken = 0;			/* index of CTGs/time token */
	unsigned matchTokens = 0;		/* % of tokens that must match exactly */
	bool is2Dim = false;
	
	const char **tokens = NULL;
	unsigned numTokens = parser[0]->parseLine(lineBuf, tokens);
	char *numStrings[3] = {NULL, NULL, NULL};
	
	if(isConvolve) {
		is2Dim = true;
		totalTokens = 7;
		ctgToken = 6;
		matchTokens = 6;
	}
    else if(nufftMode) {
        is2Dim = false;
        totalTokens = 5;
        matchTokens = 2;
        ctgToken = 4;       // not CTGs, but ms/loop, and we print it as a string 
    }
	else {
		switch(numTokens) {
			case 7:
				/* 1-dim */
				is2Dim = false;
				totalTokens = 5;
				ctgToken = 4;
				matchTokens = 2;
				break;
			case 11:
				/* 2-dim */
				is2Dim = true;
				totalTokens = 9;
				ctgToken = 8;
				matchTokens = 6;
				break;
			default:
				printf("***Malformed input file (3)\n");
				exit(1);
		}
	}
	parser[0]->freeTokens(numTokens, tokens);
	
	/* Skip ahead to first actual data in all files */
	for(unsigned dex=0; dex<numFiles; dex++) {
		if(!parser[dex]->findLine("-----------+", lineBuf)) {
			printf("***Malformed input file %s (1)\n", fileName[dex]);
			exit(1);
		}
	}
	
	if(isConvolve) {
		printf(" Image Rows | Image Cols | Kernel size ");
		for(unsigned dex=0; dex<numFiles; dex++) {
			printf("| %6s time (s) ", libType[dex]);
		}
		printf("\n");
		printf(" -----------+------------+-------------");
		for(unsigned dex=0; dex<numFiles; dex++) {
			printf("+-----------------");
		}
		printf("\n");
	}
    else if(nufftMode) {
        printf("             reference    optimized\n");
        printf("  Samples |   ms/loop   |  ms/loop    ratio\n");
        printf(" ---------+-------------+-----------+--------\n");
    }
	else if(is2Dim) {
		printf("  Width | Height | Samples ");
		for(unsigned dex=0; dex<numFiles; dex++) {
			char *title = centerStr(libType[dex], "CTGs", 11);
			printf("|%s", title);
			free(title);
		}
		printf("\n");
		printf(" -------+--------+---------");
		for(unsigned dex=0; dex<numFiles; dex++) {
			printf("+-----------");
		}
		printf("\n");
	}
	else {
		printf("  Samples ");
		for(unsigned dex=0; dex<numFiles; dex++) {
			char *title = centerStr(libType[dex], "CTGs", 11);
			printf("|%s", title);
			free(title);
		}
		printf("\n");
		printf(" ---------");
		for(unsigned dex=0; dex<numFiles; dex++) {
			printf("+-----------");
		}
		printf("\n");
	}
	
	while(1) {
		unsigned numTokensArr[numFiles];
		const char **tokensArr[numFiles];

		/* We quit when any file runs out of data */
		bool outOfData = false;
		for(unsigned dex=0; dex<numFiles; dex++) {
			if(parser[dex]->getTokens(numTokensArr[dex], tokensArr[dex])) {
				/* no more data */
				outOfData = true;
				break;
			}
			if(numTokensArr[dex] == 0) {
				/* empty line */
				outOfData = true;
				break;
			}
		}
		if(outOfData) {
			break;
		}
		
		/* verify all files have the right number of tokens */
		for(unsigned dex=0; dex<numFiles; dex++) {
			if(numTokensArr[dex] != totalTokens) {
				printf("***Malformed input file %s (2)\n", fileName[dex]);
				exit(1);
			}
		}
		
		/* the first matchTokens tokens must match exactly */
		for(unsigned fileDex=1; fileDex<numFiles; fileDex++) {
			for(unsigned tokenDex=0; tokenDex<matchTokens; tokenDex++) {
				if(strcmp(tokensArr[0][tokenDex], tokensArr[fileDex][tokenDex])) {
					/* maybe they are numbers printed with different representations.. */
					size_t num1 = fftParseStringRep(tokensArr[0][tokenDex]);
					size_t num2 = fftParseStringRep(tokensArr[fileDex][tokenDex]);
					if((num1 != 0) && (num1 != num2)) {
						printf("***Mismatched input files\n");
						exit(1);
					}
				}
			}
		}
		
		if(isConvolve) {
			printf("   %6s   |   %6s   |  %6s     ", tokensArr[0][0], tokensArr[0][2], tokensArr[0][4]);
			for(unsigned dex=0; dex<numFiles; dex++) {
				printf("|  %10s     ", tokensArr[dex][ctgToken]);
			}
			printf("\n");
		}
        else if(nufftMode) {
            /*
             * Though we print the ms/loop values as uninterpreted strings, we have to parse
             * them into floats in order to calculate the ratio.
             */
            float f0 = strtof(tokensArr[0][ctgToken], NULL);
            float f1 = strtof(tokensArr[1][ctgToken], NULL);
            float ratio = f0 / f1;
            printf(" %6s   |  %8s   | %8s  | %5.1f\n", 
                tokensArr[0][0],
                tokensArr[0][ctgToken], tokensArr[1][ctgToken], ratio);
            
        }
		else {
			/* 
			 * Convert all numbers to canonical form, no matter how they appear in
			 * input files
			 */
			size_t n = fftParseStringRep(tokensArr[0][0]);
			numStrings[0] = fftStringRepPow2(n);
			if(is2Dim) {
				n = fftParseStringRep(tokensArr[0][2]);
				numStrings[1] = fftStringRepPow2(n);
				n = fftParseStringRep(tokensArr[0][4]);
				numStrings[2] = fftStringRepPow2(n);
			}
			if(is2Dim) {
				printf(" %-6s | %-6s |  %-6s ", numStrings[0], numStrings[1], numStrings[2]);
				for(unsigned dex=0; dex<numFiles; dex++) {
					printf("| %8s  ", tokensArr[dex][ctgToken]);
				}
				printf("\n");
			}
			else {
				printf("   %-6s ", numStrings[0]);
				for(unsigned dex=0; dex<numFiles; dex++) {
					printf("| %7s   ", tokensArr[dex][ctgToken]);
				}
				printf("\n");
			}
			free(numStrings[0]);
			if(is2Dim) {
				free(numStrings[1]);
				free(numStrings[2]);
			}
		}
		for(unsigned dex=0; dex<numFiles; dex++) {
			parser[dex]->freeTokens(numTokensArr[dex], tokensArr[dex]);
		}
	}
	printf("\n");
	
	/* clean up */
	for(unsigned dex=0; dex<numFiles; dex++) {
		delete parser[dex];
	}

	return 0;
}
