//
//  GLSpectralMatrices.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 3/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef GLNumericalModelingKit_GLSpectralMatrices_h
#define GLNumericalModelingKit_GLSpectralMatrices_h

#import <GLNumericalModelingKit/Precision.h>

// diff will be populated with a differentiation matrix. This matrix can be element-wise
// multiplied with a standard order fft matrix to functionally differentiate the matrix.
// diffReals is for the real part and diffImag is for the imaginary part. shouldSwap returns
// whether or not the real part and the imaginary part need to be swapped (before applying the matrix).
void spectralDifferentiationMatrixComplex( GLSplitComplex *diff, int diffX, int diffY, int nx, int ny, GLFloat xdomain, GLFloat ydomain );


void spectralVanishingViscoityMatrix( GLFloat *q, GLFloat sigma, int nx, int ny );

#endif
