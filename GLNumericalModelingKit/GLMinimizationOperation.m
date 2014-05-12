//
//  GLMinimizationOperation.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 11/14/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import "GLMinimizationOperation.h"
#import "GLOperationOptimizer.h"

@implementation GLMinimizationOperation

- (GLMinimizationOperation *) initAtPoint: (NSArray *) startingPoints withDeltas: (NSArray *) deltas forFunction: (yFromX) functionBlock
{
	GLOperationOptimizer *functionOptimizer = [[GLOperationOptimizer alloc] initWithTopVariables: startingPoints bottomVariables: @[functionBlock(startingPoints)]];
	variableOperation functionOperation = functionOptimizer.operationBlock;
	
	NSUInteger nDims = startingPoints.count;
	NSUInteger nVertices = nDims+1;
	GLFloat ftol = 1e-5;
	
	// resultArray contains nDim GLScalar variables with the minimum point, followed by a single GLScalar variable with the value at that point.
	NSMutableArray *resultPrototypes = [NSMutableArray array];
	for (GLVariable *variable in startingPoints) {
		if (variable.isComplex || variable.rank != 0) {
			[NSException raise: @"InvalidFormat" format: @"Input must be a real scalar"];
		}
		[resultPrototypes addObject: [GLVariable variableWithPrototype: variable]];
	}
	[resultPrototypes addObject: [GLVariable variableWithPrototype: startingPoints[0]]];
	
	// operandArray contains nDim GLScalar variables defining the starting point, followed by nDim GLScalar variables defining the linearly independent directions
	NSMutableArray *operandPrototypes = [NSMutableArray arrayWithArray: startingPoints];
	[operandPrototypes addObjectsFromArray: deltas];
	
	// bufferArray contains nDim*nVertices GLScalar variables defining the vertices of the simplex.
	// followed by nVertices GLScalar variables defining y values of each of those vertices
	// followed by the buffers for the functionOperation
	NSMutableArray *buffers = [NSMutableArray array];
	for (NSUInteger i=0; i<(nDims+1)*nVertices; i++) {
		[buffers addObject: [[GLBuffer alloc] initWithLength: [startingPoints[0] dataBytes]]];
	}
	
	
	[buffers addObjectsFromArray: functionOptimizer.internalDataBuffers];
	
	variableOperation operation = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
		
		// This comes from Numerical Recipes, Third Edition, page 505-507.
        
		NSMutableArray *functionOperationBuffers = [NSMutableArray array];
		for (NSUInteger i=(nDims+1)*nVertices; i<bufferArray.count; i++) {
			[functionOperationBuffers addObject: bufferArray[i]];
		}
		
		// vertexBuffers buffers indexes the vertex (row), followed by the point (column)
		// vertices points the data in the buffers with iVertex*nDims + iDim
		GLFloat **vertices = malloc(nDims*nVertices*sizeof(GLFloat*));
		NSMutableArray *vertexBuffers = [NSMutableArray array];
		for (NSUInteger iVertex=0; iVertex < nVertices; iVertex++) {
			NSMutableArray *pointBuffers = [NSMutableArray array];
			for (NSUInteger iDim=0; iDim < nDims; iDim++) {
				NSData *theBuffer = bufferArray[iVertex*nDims + iDim];
				GLFloat *startingPoint = (GLFloat *) [operandArray[iDim] bytes];
				pointBuffers[iDim] = theBuffer;
				vertices[iVertex*nDims + iDim] = (GLFloat *) theBuffer.bytes;
				vertices[iVertex*nDims + iDim][0] = *startingPoint;
			}
			vertexBuffers[iVertex] = pointBuffers;
			if (iVertex != 0) {
				GLFloat *delta = (GLFloat *) [operandArray[nDims + iVertex-1] bytes];
				vertices[iVertex*nDims + (iVertex-1)][0] += *delta;
			}
		}
		
		// yBuffers indexes the data object that contains the y value at the vertex
		// y points to the data in that buffer
		GLFloat **y = malloc(nDims*nVertices*sizeof(GLFloat*));
		NSMutableArray *yBuffers = [NSMutableArray array];
		for (NSUInteger iVertex=0; iVertex < nVertices; iVertex++) {
			NSData *theBuffer = bufferArray[nVertices*nDims + iVertex];
			yBuffers[iVertex] = theBuffer;
			y[iVertex] = (GLFloat *) theBuffer.bytes;
			
			functionOperation(@[yBuffers[iVertex]], vertexBuffers[iVertex], functionOperationBuffers);
		}
		
		// And finally, these point to the result buffers
		GLFloat **resultVertex = malloc(nDims*nVertices*sizeof(GLFloat*));
		NSMutableArray *resultVertexBuffers = [NSMutableArray array];
		for (NSUInteger iDim=0; iDim < nDims; iDim++) {
			resultVertexBuffers[iDim] = resultArray[iDim];
			resultVertex[iDim] = (GLFloat *) [resultArray[iDim] bytes];
		}
		GLFloat *resultY = (GLFloat *) [resultArray[nDims] bytes];
		
		
		// This block should capture the memory address
		GLFloat *psum = malloc(nDims*sizeof(GLFloat));
		void (^compute_psum) (void) = ^(void) {
			for (NSUInteger iDim=0; iDim < nDims; iDim++) {
				psum[iDim] = 0.0;
				for (NSUInteger iVertex=0; iVertex < nVertices; iVertex++) {
					psum[iDim] += vertices[iVertex*nDims + iDim][0];
				}
			}
		};
		
		// This block is used to adjust the vertex points
		GLFloat (^amotry) (NSUInteger , GLFloat ) = ^(NSUInteger iHighVertex, GLFloat fraction ) {
			GLFloat fraction1 = (1.0 - fraction)/ ((GLFloat)nDims);
			GLFloat fraction2 = fraction1-fraction;
			// We are using the resultVertex memory buffers temporarily here
			for (NSUInteger iDim=0; iDim < nDims; iDim++) {
				resultVertex[iDim][0] = psum[iDim]*fraction1 - vertices[iHighVertex*nDims + iDim][0]*fraction2;
			}
			// Is this really how we want to capture this function block?
			functionOperation(@[resultArray[nDims]], resultVertexBuffers, functionOperationBuffers);
			
			// Replace the the old worst vertex with our new vertex
			if ( *resultY < y[iHighVertex][0]) {
				y[iHighVertex][0] = *resultY;
				for (NSUInteger iDim=0; iDim < nDims; iDim++) {
					psum[iDim] += resultVertex[iDim][0] - vertices[iHighVertex*nDims + iDim][0];
					vertices[iHighVertex*nDims + iDim][0] = resultVertex[iDim][0];
				}
			}
			
			return *resultY;
		};

//		GLFloat *y0 = y[0];
//		GLFloat *y1 = y[1];
//		GLFloat *y2 = y[2];
//		
//		GLFloat *p0_x = vertices[0];
//		GLFloat *p0_y = vertices[1];
//		GLFloat *p1_x = vertices[2];
//		GLFloat *p1_y = vertices[3];
//		GLFloat *p2_x = vertices[4];
//		GLFloat *p2_y = vertices[5];
        
//        GLFloat *y0 = y[0];
//		GLFloat *y1 = y[1];
//		GLFloat *y2 = y[2];
//        GLFloat *y3 = y[3];
//		
//		GLFloat *p0_x = vertices[0];
//		GLFloat *p0_y = vertices[1];
//        GLFloat *p0_z = vertices[2];
//		GLFloat *p1_x = vertices[3];
//		GLFloat *p1_y = vertices[4];
//        GLFloat *p1_z = vertices[5];
//		GLFloat *p2_x = vertices[6];
//		GLFloat *p2_y = vertices[7];
//        GLFloat *p2_z = vertices[8];
//        GLFloat *p3_x = vertices[9];
//		GLFloat *p3_y = vertices[10];
//        GLFloat *p3_z = vertices[11];
		
        GLFloat previousMin = y[0][0];
        
		compute_psum();
		NSUInteger nFunctionEvaluations = 0;
		while (1)
		{
			// First we must determine which point is the higest, next-highest, and lowest by looping over the points in the simplex
			NSUInteger ihi, inhi, ilo=0; // index of the: highest (worst), next-highest, lowest (best) vertex
			if (y[0][0] > y[1][0]) {
				ihi = 0; inhi = 1;
			} else {
				ihi = 1; inhi = 0;
			}
			for (NSUInteger iVertex=0; iVertex < nVertices; iVertex++) {
				if ( y[iVertex][0] <= y[ilo][0] ) ilo = iVertex;
				if ( y[iVertex][0] > y[ihi][0] ) {
					inhi=ihi;
					ihi=iVertex;
				}
				else if ( y[iVertex][0] > y[inhi][0] && iVertex != ihi) {
					inhi=iVertex;
				}
			}
			
			// Compute the fractional range from highest to lowest vertex
			GLFloat rtol = 2.0*fabs( y[ihi][0] - y[ilo][0] )/(fabs(y[ihi][0]) + fabs(y[ilo][0])+1.0e-10);
			
			// and return if satisfactory
			if (rtol<ftol) {
				for (NSUInteger iDim=0; iDim < nDims; iDim++) {
					resultVertex[iDim][0] = vertices[ilo*nDims + iDim][0];
				}
				*resultY = y[ilo][0];
				//NSLog(@"Total function evaluations: %lu", nFunctionEvaluations);
				
                GLFloat restarttol = 2.0*fabs( *resultY - previousMin )/(fabs(*resultY) + fabs(previousMin)+1.0e-10);
                if (restarttol > ftol) {
                    previousMin = *resultY;
                    for (NSUInteger iVertex=0; iVertex < nVertices; iVertex++) {
                        for (NSUInteger iDim=0; iDim < nDims; iDim++) {
                            vertices[iVertex*nDims + iDim][0] = resultVertex[iDim][0];
                        }
                        if (iVertex != 0) {
                            GLFloat *delta = (GLFloat *) [operandArray[nDims + iVertex-1] bytes];
                            vertices[iVertex*nDims + (iVertex-1)][0] += *delta;
                        }
                        functionOperation(@[yBuffers[iVertex]], vertexBuffers[iVertex], functionOperationBuffers);
                    }
                    compute_psum();
                    //NSLog(@"Automatic restart at %lu function evaluations.", nFunctionEvaluations);
                    continue;
                }
                else {
                    //NSLog(@"Total function evaluations: %lu", nFunctionEvaluations);
                    free(vertices); free(y); free(psum); free(resultVertex);
                    break;
                }
			}
			
			if (nFunctionEvaluations > 15000) {
				NSLog(@"Exceed the maximum number of function evaluations");
				
				free(vertices); free(y); free(psum); free(resultVertex);
				break;
			}
			nFunctionEvaluations+=2;
			
			// Begin a new iteration. First extrapolate by a factor -1 through the face of the sinplex across from the high point, i.e., reflect the simplex from the high point.
			GLFloat ytry = amotry( ihi, -1.0);
			if (ytry <= y[ilo][0]) {
				// Gives a result better than the best point, so try an additional extrapolation by a factor of 2.
				ytry = amotry( ihi, 2.0);
			} else if (ytry >= y[inhi][0]) {
				// the reflected point is worse than the second-highest, so look for an intermediate lower point, i.e., do a one-dimensional contraction
				GLFloat ysave = y[ihi][0];
				ytry = amotry( ihi, 0.5);
				if (ytry >= ysave) { // Can't seem to get rid of that high point
					// Better contract around the lowest (best) point.
					for (NSUInteger iVertex=0; iVertex < nVertices; iVertex++) {
						if (iVertex != ilo) {
							for (NSUInteger iDim=0; iDim < nDims; iDim++) {
								vertices[iVertex*nDims + iDim][0] = 0.5*( vertices[iVertex*nDims + iDim][0] + vertices[ilo*nDims + iDim][0] );
								functionOperation(@[yBuffers[iVertex]], vertexBuffers[iVertex], functionOperationBuffers);
							}
						}
					}
					nFunctionEvaluations += nDims;
					compute_psum();
				}
			} else {
				nFunctionEvaluations--;
			}
			
		}
	};
	
	if ((self = [super initWithResult: resultPrototypes operand: operandPrototypes buffers: buffers operation:operation])) {
		
	}
	
	return self;
}

@end
