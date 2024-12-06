//
//  mexPlan.c
//  
//
//  Created by Jeffrey Early on 12/5/24.
//

#include <stdio.h>
#include "fftw3.h"
#include <mex.h>

void  mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, j, bw, bw2_1, size, size2_1, nrow, ncol;
    int data_is_real;
    int cutoff;
    int rank, howmany_rank;
    double *rresult, *iresult, *rdata, *idata;
    double *workspace, *weights;

    fftw_plan dctPlan;
    fftw_plan fftPlan;
    fftw_iodim dims[1], howmany_dims[1];

    bw = 2;
    weights = (double *)malloc(sizeof(double) * 4 * bw);
    rdata = (double *)malloc(sizeof(double) * 5 * bw);
    dctPlan = fftw_plan_r2r_1d(2 * bw, weights, rdata,
        FFTW_REDFT10, FFTW_ESTIMATE);
}
