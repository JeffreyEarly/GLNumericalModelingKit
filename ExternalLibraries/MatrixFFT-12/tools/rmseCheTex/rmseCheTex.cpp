/*
 * Copyright (c) 2009 Apple, Inc. All Rights Reserved.
 * 
 * rmseCheTex.cpp - convert the output of the runRmseChe script to TeX source.
 *
 * Created Aug. 12 2009. 
 *
 *
 * To use this, run the runRmseChe script twice, once for each of (single, double) precision. 
 * It doesn't matter what complex config - split or interleaved - is selected. 
 * Send the outputs of the two runRmseChe runs to two separate text files. 
 * Run this program, give it the names of those two files. It doesn't matter which 
 * order you specify the files in - this program figures out which one is Double and 
 * which is Single precision. 
 *
 * The output of this program can be copy&psted right into the TeX source in two chunks 
 * delineated by big 3-line TeX comments, also emitted here). 
 */
 
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <libMatrixFFT/src/TextParser.h>		/* private interface */
#include <libMatrixFFT/fftUtils.h>

typedef enum {
    PR_Double = 0,
    PR_Single
} Precision;

static void usage(char **argv)
{
	printf("usage: %s runRmseChe_outfile1 runRmseChe_outfile2\n", argv[0]);
	exit(1);
}

/* 
 * Convert a string containing an ASCII integer, a power of 2, in either decimal 
 * or 2^n form, to log2(n) 
 */
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
 * Convert a printf-style "2.8e-17" to TeX-style 
 * "$2.8 \times 10^{-17}$"
 */
static void printTexFloat(
    const char *inf)        /* e.g. "2.8e-17" */
{
    char mant[100];
    char *outp = mant;
    
    for(;;) {
        char c = *inf++;
        if(c == 'e') {
            break;
        }
        *outp++ = c;
    }
    
    /* NULL-terminate mantissa */
    *outp = '\0';

    /* 
     * *inf points to NULL-terminated exponent.
     * A simple way to strip off leading zeroes...
     */
    int expo = atoi(inf);
    
    printf("$%s \\times 10^{%d}$", mant, expo);
}

/*
  1-D real. 
  Input (cursor is currently on the "underline" line, after the column headings):
 
   samples  | RMSE round | MXE round | RMSE fwd |  MXE fwd  
 -----------+------------+-----------+----------+-----------
    2^10        2.8e-17     3.3e-16     1.2e-13    2.1e-12
    2^20        3.2e-18     1.8e-15     1.4e-13    4.7e-12
    2^30        3.2e-18     1.8e-15     1.4e-13    4.7e-12


  Output:
                N         end-end RMSE            end-end MXE              fwd RMSE               fwd MXE
 1D Double & $2^{10}$ & $2.8 \times 10^{-17}$ & $3.3 \times 10^{-16}$ & $1.2 \times 10^{-13}$ & $2.1 \times 10^{-12}$ & --- \\
 1D Double & $2^{20}$ & $1.5 \times 10^{-18}$ & $4.4 \times 10^{-16}$ & $1.4 \times 10^{-13}$ & $4.7 \times 10^{-12}$ & --- \\
 1D Double & $2^{31}$ & $7.1 \times 10^{-20}$ & $1.8 \times 10^{-15}$ & $1.4 \times 10^{-13}$ & $5.2 \times 10^{-12}$ & --- \\
*/
static int process1DReal(
    TextParser *parser,
    Precision prec)
{
    unsigned numTokens = 0;
    const char **tokens = NULL;
    const char *precStr = (prec == PR_Double) ? "Double" : "Single";
    
    /*
     * Cursor at underline; skip 
     */
    if(parser->skipLine()) {
        printf("***Unexpected EOF at %s\n", parser->fileName());
        return -1;
    }
    for(;;) {
        if(parser->getTokens(numTokens, tokens)) {
            /* EOF */
            break;
        }
        if(numTokens != 5) {
            /* done */
            break;
        }
        
        /* line header */
        printf("1D %s & ", precStr);
        
        /* N */
        unsigned log2N = strToLog2(tokens[0]);
        printf("$2^{%u}$ & ", log2N);
        
        /* meat & potatoes, same order in input and output */
        printTexFloat(tokens[1]);
        printf(" & ");
        printTexFloat(tokens[2]);
        printf(" & ");
        printTexFloat(tokens[3]);
        printf(" & ");
        printTexFloat(tokens[4]);
        printf(" & --- \\\\\n");

        parser->freeTokens(numTokens, tokens);
    }   
    if(tokens != NULL) {
        parser->freeTokens(numTokens, tokens);
    }
    return 0;
}

/*
   Complex, 1-D and 2-D.
  
   Input, 1-D (cursor is currently on the "underline" line, after the column headings):
   
   samples  | RMSE round | MXE round |   CHE    
 -----------+------------+-----------+----------
    2^10        1.5e-07     4.9e-07    3.3e-07
    2^20        2.7e-07     1.4e-06    7.2e-07
    2^30        2.7e-07     1.4e-06    7.2e-07

   Input, 2-D (cursor is currently on the "underline" line, after the column headings):
   
  rows | cols | samples | RMSE round | MXE round |   CHE   
 ------+------+---------+------------+-----------+----------
  2^5    2^5      2^10     1.1e-07      1.9e-07    1.2e-07
  2^10   2^10     2^20     1.8e-07      7.1e-07    3.9e-07
  2^15   2^15     2^30     1.8e-07      7.1e-07    3.9e-07
  
    Output: 
                N         end-end RMSE            end-end MXE                          CHE   
 1D Double & $2^{10}$ & $3.0 \times 10^{-16}$ & $8.0 \times 10^{-16}$ & --- & --- & $1.4 \times 10^{-13}$ \\
 1D Double & $2^{20}$ & $7.0 \times 10^{-16}$ & $3.0 \times 10^{-15}$ & --- & --- & $2.3 \times 10^{-10}$ \\
 1D Double & $2^{30}$ & $7.8 \times 10^{-16}$ & $5.7 \times 10^{-15}$ & --- & --- & $3.6 \times 10^{-7}$ \\

*/

static int processComplex(
    TextParser *parser,
    Precision prec,
    bool is2D)
{
    unsigned numTokens = 0;
    const char **tokens = NULL;
    const char *precStr = (prec == PR_Double) ? "Double" : "Single";
    const char *dimStr = is2D ? "2D" : "1D";
    
    /*
     * Cursor at underline; skip 
     */
    if(parser->skipLine()) {
        printf("***Unexpected EOF at %s\n", parser->fileName());
        return -1;
    }

    unsigned NDex = 0;
    unsigned RMSEDex = 0;
    unsigned MXEDex = 0;
    unsigned CHEDex = 0;
    unsigned expectTokens = 0;
    
    if(is2D) {
        expectTokens = 6;
        NDex = 2;
        RMSEDex = 3;
        MXEDex = 4;
        CHEDex = 5;
    }
    else {
        expectTokens = 4;
        NDex = 0;
        RMSEDex = 1;
        MXEDex = 2;
        CHEDex = 3;
    }
    
    for(;;) {
        if(parser->getTokens(numTokens, tokens)) {
            /* EOF */
            break;
        }
        if(numTokens != expectTokens) {
            /* done */
            break;
        }
        
        /* line header */
        printf("%s %s & ", dimStr, precStr);
        
        /* N */
        unsigned log2N = strToLog2(tokens[NDex]);
        printf("$2^{%u}$ & ", log2N);
        
        printTexFloat(tokens[RMSEDex]);
        printf(" & ");
        printTexFloat(tokens[MXEDex]);
        printf(" & --- & --- & ");
        printTexFloat(tokens[CHEDex]);
        printf("  \\\\\n");

        parser->freeTokens(numTokens, tokens);
    }   
    if(tokens != NULL) {
        parser->freeTokens(numTokens, tokens);
    }
    return 0;

}

int main(int argc, char **argv)
{
    if(argc != 3) {
        usage(argv);
    }
    
    /* these exit on error */
    TextParser parser1(argv[1]);
    TextParser parser2(argv[2]);
    
    /* 
     * Figure out which input file is single precision.
     * The precision line looks this this:
     *
     * Precision   : Double
     */
    char lineBuf[LINE_LENGTH_MAX];
    TextParser *singleParser = NULL;
    TextParser *doubleParser = NULL;
    
    if(!parser1.findLine("Precision", lineBuf)) {
        printf("***%s is not a runRmseChe output file (1). Aborting.\n", argv[1]);
        exit(1);
    }
    
    const char **tokens = NULL;
    unsigned numTokens;
    numTokens = parser1.parseLine(lineBuf, tokens);
    if(numTokens != 3) {
        printf("***%s is not a runRmseChe output file (2). Aborting.\n", argv[1]);
        exit(1);
    }
    if(!strcmp(tokens[2], "Double")) {
        doubleParser = &parser1;
        singleParser = &parser2;
    }
    else {
        doubleParser = &parser2;
        singleParser = &parser1;
    }
    parser1.freeTokens(numTokens, tokens);
    
    /*
     * 1-D real: double then single
     */
    printf("%%%%%%\n");
    printf("%%%%%% 1-D real\n");
    printf("%%%%%%\n");

    if(!doubleParser->findLine("samples", lineBuf)) {
        printf("***%s is not a runRmseChe output file (3). Aborting.\n", doubleParser->fileName());
        exit(1);
    }
    if(process1DReal(doubleParser, PR_Double)) {
        exit(1);
    }
    if(!singleParser->findLine("samples", lineBuf)) {
        printf("***%s is not a runRmseChe output file (4). Aborting.\n", singleParser->fileName());
        exit(1);
    }
    if(process1DReal(singleParser, PR_Single)) {
        exit(1);
    }
    
    /*
     * 1D Double
     */
    printf("%%%%%%\n");
    printf("%%%%%% Complex\n");
    printf("%%%%%%\n");
    
    if(!doubleParser->findLine("samples", lineBuf)) {
        printf("***%s is not a runRmseChe output file (5). Aborting.\n", doubleParser->fileName());
        exit(1);
    }
    if(processComplex(doubleParser, PR_Double, false)) {
        exit(1);
    }
    
    /*
     * 1D single
     */
    printf("\\hline\n");
    if(!singleParser->findLine("samples", lineBuf)) {
        printf("***%s is not a runRmseChe output file (6). Aborting.\n", singleParser->fileName());
        exit(1);
    }
    if(processComplex(singleParser, PR_Single, false)) {
        exit(1);
    }

    /*
     * 2D double
     */
    printf("\\hline\n");
    if(!doubleParser->findLine("rows", lineBuf)) {
        printf("***%s is not a runRmseChe output file (7). Aborting.\n", doubleParser->fileName());
        exit(1);
    }
    if(processComplex(doubleParser, PR_Double, true)) {
        exit(1);
    }

    /*
     * 2D single
     */
    printf("\\hline\n");
    if(!singleParser->findLine("rows", lineBuf)) {
        printf("***%s is not a runRmseChe output file (7). Aborting.\n", singleParser->fileName());
        exit(1);
    }
    if(processComplex(singleParser, PR_Single, true)) {
        exit(1);
    }

}
