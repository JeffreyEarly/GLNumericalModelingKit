//
//  simple_add.c
//  
//
//  Created by Jeffrey Early on 12/8/24.
//

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    // Check for proper number of input and output arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:mexFunction:nargin", "Two input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgIdAndTxt("MyToolbox:mexFunction:nargout", "Too many output arguments.");
    }

    // Get the dimensions of the input matrices
    const mwSize *dims1 = mxGetDimensions(prhs[0]);
    const mwSize *dims2 = mxGetDimensions(prhs[1]);

    // Check if the input matrices have the same dimensions
    if (dims1[0] != dims2[0] || dims1[1] != dims2[1]) {
        mexErrMsgIdAndTxt("MyToolbox:mexFunction:inputDimensions", "Input matrices must have the same dimensions.");
    }

    // Create the output matrix
    plhs[0] = mxCreateNumericArray(2, dims1, mxDOUBLE_CLASS, mxREAL);

    // Get pointers to the input and output data
    double *A = mxGetPr(prhs[0]);
    double *B = mxGetPr(prhs[1]);
    double *C = mxGetPr(plhs[0]);

    // Perform the matrix addition
    for (mwIndex i = 0; i < dims1[0] * dims1[1]; ++i) {
        C[i] = A[i] + B[i];
    }
}
