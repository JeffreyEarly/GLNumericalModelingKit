//
//  GLSpectralMatrices.c
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 3/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include "GLSpectralMatrices.h"

/************************************************/
/*		Direct Matrix Creation					*/
/************************************************/

#pragma mark -
#pragma mark Direct Matrix Creation
#pragma mark

GLFloat **matrixPointerToVector( GLFloat *vec, int nx, int ny ) {
	GLFloat **matrix = (GLFloat **) malloc( nx * sizeof( GLFloat *) );
	int i;
	for ( i=0; i<nx; i++ ) {
		matrix[i] = &vec[i*ny];
	}
	
	return matrix;
}

void spectralDifferentiationMatrixComplex( GLSplitComplex *diff, int diffX, int diffY, int nx, int ny, GLFloat xdomain, GLFloat ydomain )
{	
	// For convinience
	GLFloat **realMatrix = matrixPointerToVector( diff->realp, nx, ny);
	GLFloat **imagMatrix = matrixPointerToVector( diff->imagp, nx, ny);
	
	// Determine the appropriate signs. Careful bookkeeping required!!!
	// These first two are for positive frequencies in both X and Y directions
	// a*i^0 (+ x + i y) = (+ a + i 0) (+ x + i y) = ( +a x +a i y)
	// a*i^1 (+ x + i y) = (+ 0 + i a) (+ x + i y) = ( -a y +a i x)
	// a*i^2 (+ x + i y) = (- a + i 0) (+ x + i y) = ( -a x -a i y)
	// a*i^3 (+ x + i y) = (+ 0 - i a) (+ x + i y) = ( +a y -a i x)
	
	GLFloat realSign = ((diffX+diffY) % 2) == 1 ? 0.0 : 1.0; // Odd powers result in no real part
	realSign *= ((diffX+diffY) % 4) == 2 ? -1.0 : 1.0; // Multiples of 2 result in a negative real part.
	
	GLFloat imagSign = ((diffX+diffY) % 2) == 0 ? 0.0 : 1.0; // Even powers result in no imaginary part
	imagSign *= ((diffX+diffY) % 4) == 3 ? -1.0 : 1.0; // Multiples of 3 result in a negative imaginary part.
	
	// The additional sign change when dealing with negative frequencies.
	GLFloat negXSign = ( (diffX % 2) == 1 ? -1.0 : 1.0 );
	GLFloat negYSign = ( (diffY % 2) == 1 ? -1.0 : 1.0 );
	
	GLFloat xexponent = 2.0 * M_PI / ( xdomain );
	GLFloat yexponent = 2.0 * M_PI / ( ydomain );
	
	int i,j;
	for(i=1; i < nx/2; i++) {
		GLFloat xFactor = pow( xexponent * ( (GLFloat) i ), (GLFloat) diffX);
		for(j=1; j < ny/2; j++) {
			GLFloat xyFactor = xFactor * pow( yexponent * ( (GLFloat) j ), (GLFloat) diffY);
			// positive frequency values for X, positive frequency values for Y
			realMatrix[i][j] = realSign * xyFactor;
			imagMatrix[i][j] = imagSign * xyFactor;
			
			// negative frequency values for X, positive frequency values for Y
			realMatrix[nx-i][j] = negXSign * realSign * xyFactor;
			imagMatrix[nx-i][j] = negXSign * imagSign * xyFactor;
			
			// positive frequency values for X, negative frequency values for Y
			realMatrix[i][ny-j] = negYSign * realSign * xyFactor;
			imagMatrix[i][ny-j] = negYSign * imagSign * xyFactor;
			
			// negative frequency values for X, negative frequency values for Y
			realMatrix[nx-i][ny-j] = negXSign * negYSign * realSign * xyFactor;
			imagMatrix[nx-i][ny-j] = negXSign * negYSign * imagSign * xyFactor;
		}
	}
	
	//Now deal with the 0th and N/2 terms	
	for(j=1; j < ny/2; j++) {
		GLFloat xFactor = (diffX == 0 ? 1.0 : 0.0);
		GLFloat yFactor = pow( yexponent * ( (GLFloat) j ), (GLFloat) diffY);
		realMatrix[0][j] = realSign * yFactor * xFactor;
		imagMatrix[0][j] = imagSign * yFactor * xFactor;
		realMatrix[0][ny-j] = negYSign * realSign * yFactor * xFactor;
		imagMatrix[0][ny-j] = negYSign * imagSign * yFactor * xFactor;
		
		// IOW, it's zero for odd differentiations, even otherwise
		xFactor = ((GLFloat) ((diffX + 1)%2)) * pow( xexponent * ( (GLFloat) nx/2 ), (GLFloat) diffX);
		realMatrix[nx/2][j] = negXSign * realSign * yFactor * xFactor;
		imagMatrix[nx/2][j] = negXSign * imagSign * yFactor * xFactor;
		realMatrix[nx/2][ny-j] = negXSign * negYSign * realSign * yFactor * xFactor;
		imagMatrix[nx/2][ny-j] = negXSign * negYSign * imagSign * yFactor * xFactor;
	}
	
	for(i=1; i < nx/2; i++) {
		GLFloat xFactor = pow( xexponent * ( (GLFloat) i ), (GLFloat) diffX);
		GLFloat yFactor = (diffY == 0 ? 1.0 : 0.0);
		realMatrix[i][0] = realSign * xFactor * yFactor;
		imagMatrix[i][0] = imagSign * xFactor * yFactor;
		realMatrix[nx-i][0] = negXSign * realSign * xFactor * yFactor;
		imagMatrix[nx-i][0] = negXSign * imagSign * xFactor * yFactor;
		
		yFactor = ((GLFloat) ((diffY + 1)%2)) * pow( yexponent * ( (GLFloat) ny/2 ), (GLFloat) diffY);
		realMatrix[i][ny/2] = negYSign * realSign * xFactor * yFactor;
		imagMatrix[i][ny/2] = negYSign * imagSign * xFactor * yFactor;
		realMatrix[nx-i][ny/2] = negXSign * negYSign * realSign * xFactor * yFactor;
		imagMatrix[nx-i][ny/2] = negXSign * negYSign * imagSign * xFactor * yFactor;
	}
	
	realMatrix[0][0] = ( (diffX + diffY > 0 ) ? 0.0 : 1.0 );
	imagMatrix[0][0] = ( (diffX + diffY > 0 ) ? 0.0 : 1.0 );
	
	GLFloat xFactor = ((GLFloat) ((diffX + 1)%2)) * pow( xexponent * ((GLFloat) nx/2), (GLFloat) diffX);
	GLFloat yFactor = ( (diffY > 0 ) ? 0.0 : 1.0 );
	realMatrix[nx/2][0] = negXSign * realSign * xFactor * yFactor;
	imagMatrix[nx/2][0] = negXSign * imagSign * xFactor * yFactor;
	
	xFactor = ( (diffX > 0 ) ? 0.0 : 1.0 );
	yFactor = ((GLFloat) ((diffY + 1)%2)) * pow( yexponent * ( (GLFloat) ny/2 ), (GLFloat) diffY);
	realMatrix[0][ny/2] = negYSign * realSign * xFactor * yFactor;
	imagMatrix[0][ny/2] = negYSign * imagSign * xFactor * yFactor;
	
	xFactor = ((GLFloat) ((diffX + 1)%2)) * pow(xexponent * ((GLFloat) nx/2), (GLFloat) diffX);
	yFactor = ((GLFloat) ((diffY + 1)%2)) * pow( yexponent * ( (GLFloat) ny/2 ), (GLFloat) diffY);	
	realMatrix[nx/2][ny/2] = negXSign * negYSign * realSign * xFactor * yFactor;
	imagMatrix[nx/2][ny/2] = negXSign * negYSign * imagSign * xFactor * yFactor;
	
	free(realMatrix);
	free(imagMatrix);
	
	return;
}

void spectralVanishingViscoityMatrix( GLFloat *q, GLFloat sigma, int nx, int ny )
{
	GLFloat **qMatrix = matrixPointerToVector( q, nx, ny);
	GLFloat alpha = 1.0;
	GLFloat p = 2.0;
	GLFloat max = 1.0; //sigma + 0.2;
	
	int i,j;
	for(i=1; i < nx/2; i++) {
		for(j=1; j < ny/2; j++) {
			GLFloat k = sqrt( ( 2.0*((GLFloat)i)/((GLFloat)nx) )*( 2.0*((GLFloat)i)/((GLFloat)nx) ) + ( 2.0*((GLFloat)j)/((GLFloat)ny) )*( 2.0*((GLFloat)j)/((GLFloat)ny) ));
			qMatrix[i][j] = qMatrix[nx-i][j] = qMatrix[i][ny-j] = qMatrix[nx-i][ny-j] = k < sigma ? 0.0 : expf( -alpha * powf( (k-max)/(k-sigma), p) );
			if (k > max) qMatrix[i][j] = qMatrix[nx-i][j] = qMatrix[i][ny-j] = qMatrix[nx-i][ny-j] = 1.0;
		}
	}
	
	//Now deal with the 0th and N/2 terms	
	for(j=1; j < ny/2; j++) {
		GLFloat k = ( 2.0*((GLFloat)j)/((GLFloat)ny) );
		qMatrix[0][j] = qMatrix[0][ny-j] = k < sigma ? 0.0 : expf( -alpha * powf( (k-max)/(k-sigma), p) );
		if (k > max) qMatrix[0][j] = qMatrix[0][ny-j] = 1.0;
		
		k = sqrt( 1.0 + ( 2.0*((GLFloat)j)/((GLFloat)ny) )*( 2.0*((GLFloat)j)/((GLFloat)ny) ));
		qMatrix[nx/2][j] = qMatrix[nx/2][ny-j] = k < sigma ? 0.0 : expf( -alpha * powf( (k-max)/(k-sigma), p) );
		if (k > max) qMatrix[nx/2][j] = qMatrix[nx/2][ny-j] = 1.0;
	}
	
	for(i=1; i < nx/2; i++) {
		GLFloat k = ( 2.0*((GLFloat)i)/((GLFloat)nx) );
		qMatrix[i][0] = qMatrix[nx-i][0] = k < sigma ? 0.0 : expf( -alpha * powf( (k-max)/(k-sigma), p) );
		if (k > max) qMatrix[i][0] = qMatrix[nx-i][0] = 1.0;
		
		k = sqrt( ( 2.0*((GLFloat)i)/((GLFloat)nx) )*( 2.0*((GLFloat)i)/((GLFloat)nx) ) + 1.0 );
		qMatrix[i][ny/2] = qMatrix[nx-i][ny/2] = k < sigma ? 0.0 : expf( -alpha * powf( (k-max)/(k-sigma), p) );
		if (k > max) qMatrix[i][ny/2] = qMatrix[nx-i][ny/2] = 1.0;
	}
	
	qMatrix[0][0] = 0.0;
	qMatrix[nx/2][0] = 1.0;
	qMatrix[0][ny/2] = 1.0;
	qMatrix[nx/2][ny/2] = 1.0;
	free(qMatrix);
}