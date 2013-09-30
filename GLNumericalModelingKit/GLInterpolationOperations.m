//
//  GLInterpolationOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 5/2/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLInterpolationOperations.h"
#import "GLMemoryPool.h"

BOOL noWrap = NO;
BOOL xWrap = NO;
BOOL allWrap = YES;

// This operation has the simple task of converting points along a dimension, to fractional indices.
@interface GLLinearInterpolationIndexOperation : GLUnaryOperation
- (id) initWithOperand:(GLVariable *)positionVar dimension: (GLDimension *)	dimension;
@end

@implementation GLLinearInterpolationIndexOperation

NSInteger indexBelow( GLFloat *monotonicallyIncreasingValues, GLFloat value, NSInteger imin, NSInteger imax )
{
	if (imin+1 == imax) {
		return imin;
	}
	else
	{
		NSInteger imid=imin+(imax-imin)/2;
		if (monotonicallyIncreasingValues[imid] > value) {
			return indexBelow(monotonicallyIncreasingValues, value, imin, imid);
		}
		else {
			return indexBelow(monotonicallyIncreasingValues, value, imid, imax);
		}
	}
}

- (id) initWithOperand:(GLVariable *)positionVar dimension: (GLDimension *)	dimension
{
	if ((self=[super initWithOperand: positionVar]))
	{
		GLFloat minIndex = 0.0;
		GLFloat maxIndexX = dimension.nPoints-1;
		GLFloat nx = dimension.nPoints;
		
		NSUInteger numInterpPoints = positionVar.nDataPoints;
		
		if (dimension.isEvenlySampled)
		{
			GLFloat xDimInvInterval = 1.0/dimension.sampleInterval;
			GLFloat xDimOffset = - dimension.domainMin / dimension.sampleInterval;
			
			if (noWrap)
			{
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					vGL_vsmsa( (GLFloat *)operand.bytes, 1, (GLFloat *) &xDimInvInterval, (GLFloat *) &xDimOffset, result.mutableBytes, 1, numInterpPoints);
					vGL_vclip( result.mutableBytes, 1, (GLFloat *)&minIndex,(GLFloat *) &maxIndexX, result.mutableBytes, 1, numInterpPoints);
				};
			}
			else if (xWrap || allWrap)
			{
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					vGL_vsmsa( (GLFloat *)operand.bytes, 1, (GLFloat *) &xDimInvInterval, (GLFloat *) &xDimOffset, result.mutableBytes, 1, numInterpPoints);
					
					GLFloat *f = (GLFloat *) result.mutableBytes;
					for (NSUInteger i=0; i<numInterpPoints; i++) {
						if (f[i]>nx) {
							f[i]=fmodf(f[i], nx);
						} else if (f[i]<0) {
							f[i]=nx+fmodf(f[i], nx);
						}
					}
				};
			}
		}
		else
		{
			NSData *dimData = dimension.data;
			dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
			
			GLFloat domainMin = dimension.domainMin;
			GLFloat domainLength = dimension.domainLength;
			
			if (noWrap)
			{
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLFloat *monotonicallyIncreasingValues = (GLFloat *) dimData.bytes;
					GLFloat *values = (GLFloat *) operand.bytes;
					GLFloat *fractionalIndex = (GLFloat *) result.bytes;
					
					dispatch_apply(numInterpPoints, globalQueue, ^(size_t i) {
						NSInteger index = indexBelow(monotonicallyIncreasingValues, values[i], 0, nx);
						if (index == nx) {
							fractionalIndex[i] = index;
						} else {
							fractionalIndex[i] = index + (values[i] - monotonicallyIncreasingValues[index])/(monotonicallyIncreasingValues[index+1] - monotonicallyIncreasingValues[index]);
						}
					});
				};
			}
			else if (xWrap || allWrap)
			{
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLFloat *monotonicallyIncreasingValues = (GLFloat *) dimData.bytes;
					GLFloat *values = (GLFloat *) operand.bytes;
					GLFloat *fractionalIndex = (GLFloat *) result.bytes;
					
					for (NSUInteger i=0; i<numInterpPoints; i++) {
						if (values[i]>domainMin+domainLength) {
							values[i]=domainMin+fmodf(values[i]-domainMin, domainLength);
						} else if (values[i]<domainMin) {
							values[i]=domainMin+domainLength+fmodf(values[i]-domainMin, domainLength);
						}
					}
		
					dispatch_apply(numInterpPoints, globalQueue, ^(size_t i) {
						NSInteger index = indexBelow(monotonicallyIncreasingValues, values[i], 0, nx);
						if (index == nx) {
							fractionalIndex[i] = index;
						} else {
							fractionalIndex[i] = index + (values[i] - monotonicallyIncreasingValues[index])/(monotonicallyIncreasingValues[index+1] - monotonicallyIncreasingValues[index]);
						}
					});
				};
			}
		}
	}
	return self;
}



@end

/************************************************/
/*		GLInterpolationOperation				*/
/************************************************/

// Given a vector of x positions and y positions, as well as the x dimension and y dimension,
// this function returns four variables:
// 1) a vector of indices i*ny+j (incl. j fraction) corresponding to the nearest "floor"/lower-left/truncated index in x, but not j
// 2) a vector of the (i+1)*ny+j (incl. j fraction) indices
// 3) a vector, xFrac, of the relative distance the position is between the i and i+1 index. Value is between (0,1).
// 4) a vector with 1-xFrac.
// This is designed as the input for the bilinear interpolation algorithm.

@interface GLBilinearInterpolationIndexOperation : GLUnaryVectorOperation
@end

@implementation GLBilinearInterpolationIndexOperation

- (id) initWithOperand:(NSArray *)operandArray dimensions: (NSArray *) dimensions
{	
	if (dimensions.count != 2 && operandArray.count !=2) {
		[NSException raise: @"MethodNotImplemented" format: @"Interpolation can only be done in two dimensions."];
	}
	
	GLVariable *xPositionVar = [operandArray objectAtIndex: 0];
	GLVariable *firstInterpIndicesVar = [GLVariable variableOfRealTypeWithDimensions:xPositionVar.dimensions forEquation:xPositionVar.equation];
	GLVariable *secondInterpIndicesVar = [GLVariable variableOfRealTypeWithDimensions:xPositionVar.dimensions forEquation:xPositionVar.equation];
	GLVariable *xFracVar = [GLVariable variableOfRealTypeWithDimensions:xPositionVar.dimensions forEquation:xPositionVar.equation];
	GLVariable *xFracOneMinusVar = [GLVariable variableOfRealTypeWithDimensions:xPositionVar.dimensions forEquation:xPositionVar.equation];
	
	NSArray *resultArray = @[firstInterpIndicesVar, secondInterpIndicesVar, xFracVar, xFracOneMinusVar];
	
	if ((self=[super initWithResult: resultArray operand: operandArray]))
	{
		GLDimension *xDim = [dimensions objectAtIndex: 0];
		GLDimension *yDim = [dimensions objectAtIndex: 1];
		
		GLFloat xDimInvInterval = 1.0/xDim.sampleInterval;
		GLFloat xDimOffset = - xDim.domainMin / xDim.sampleInterval;
		
		GLFloat yDimInvInterval = 1.0/yDim.sampleInterval;
		GLFloat yDimOffset = - yDim.domainMin / yDim.sampleInterval;
		
		GLFloat nx = xDim.nPoints;
		GLFloat ny = yDim.nPoints;
		GLFloat minIndex = 0.0;
		GLFloat maxIndexX = xDim.nPoints-1;
		GLFloat maxIndexY = yDim.nPoints-1;
		
		NSUInteger numInterpPoints = xPositionVar.nDataPoints;
		
		dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
		
		if (noWrap)
		{
			self.blockOperation = ^(NSArray *result, NSArray *operand)
			{
				NSMutableData *xPositions = [operand objectAtIndex: 0];
				NSMutableData *yPositions = [operand objectAtIndex: 1];
				
				NSMutableData *firstInterpIndices = [result objectAtIndex: 0];
				NSMutableData *secondInterpIndices = [result objectAtIndex: 1];
				NSMutableData *xFrac = [result objectAtIndex: 2];
				NSMutableData *xFracOneMinus = [result objectAtIndex: 3];
				
				dispatch_group_t my_group = dispatch_group_create();
				dispatch_group_t finalGroup = dispatch_group_create();

				dispatch_group_enter(finalGroup);
				dispatch_group_enter(finalGroup);
				
				// Convert from a position to an index. index = (position - domainMin)/sampleInterval
				// and ignore peridiocity -- simply clip the upper and lower indices.
				dispatch_group_async( my_group, globalQueue, ^{
					// D = A*b + c
					vGL_vsmsa( xPositions.mutableBytes, 1, (GLFloat *) &xDimInvInterval, (GLFloat *) &xDimOffset, firstInterpIndices.mutableBytes, 1, numInterpPoints);
					vGL_vclip( firstInterpIndices.mutableBytes, 1, (GLFloat *)&minIndex,(GLFloat *) &maxIndexX, firstInterpIndices.mutableBytes, 1, numInterpPoints);
					
					vGL_vfrac( firstInterpIndices.mutableBytes, 1, xFrac.mutableBytes, 1, numInterpPoints);
					// C = B - A
					vGL_vsub( xFrac.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, numInterpPoints);
									
					// firstInterpIndices now contains an integer index (i) and xFrac contains the remainder,
					// so it's safe to go compute the xFrac scaling factor we need on our own.
					dispatch_async( globalQueue, ^{
						// xfrac = 1-xfrac
						GLFloat one = 1;
						vGL_vneg( xFrac.mutableBytes, 1, xFracOneMinus.mutableBytes, 1, numInterpPoints);
						vGL_vsadd( xFracOneMinus.mutableBytes, 1, &one, xFracOneMinus.mutableBytes, 1, numInterpPoints );
						
						dispatch_group_leave(finalGroup);
					});
					
					// firstInterpIndices = ny*i
					vGL_vsmul( firstInterpIndices.mutableBytes, 1,  &ny, firstInterpIndices.mutableBytes, 1, numInterpPoints );
				});
				
				// Go compute the y index value
				dispatch_group_async( my_group, globalQueue, ^{
					vGL_vsmsa( yPositions.mutableBytes, 1,(GLFloat *) &yDimInvInterval,(GLFloat *) &yDimOffset, secondInterpIndices.mutableBytes, 1, numInterpPoints);
					vGL_vclip( secondInterpIndices.mutableBytes, 1, (GLFloat *)&minIndex, (GLFloat *)&maxIndexY, secondInterpIndices.mutableBytes, 1, numInterpPoints);
				});	
				
				dispatch_group_notify(my_group, globalQueue, ^{
					// firstInterpIndices = ny*i + j (including j fraction)
					vGL_vadd( firstInterpIndices.mutableBytes, 1, secondInterpIndices.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, numInterpPoints);
					// secondInterpIndices = ny*(i+1) + j (including j fraction)
					vGL_vsadd( firstInterpIndices.mutableBytes, 1, (void *) &ny, secondInterpIndices.mutableBytes, 1, numInterpPoints);
					
					dispatch_group_leave(finalGroup);
				});
				
				dispatch_group_wait(finalGroup, DISPATCH_TIME_FOREVER);
			};
		}
		else if (xWrap)
		{
			maxIndexX = xDim.nPoints;
			maxIndexY = yDim.nPoints-1;
			self.blockOperation = ^(NSArray *result, NSArray *operand)
			{
				NSMutableData *xPositions = [operand objectAtIndex: 0];
				NSMutableData *yPositions = [operand objectAtIndex: 1];
				
				NSMutableData *firstInterpIndices = [result objectAtIndex: 0];
				NSMutableData *secondInterpIndices = [result objectAtIndex: 1];
				NSMutableData *xFrac = [result objectAtIndex: 2];
				NSMutableData *xFracOneMinus = [result objectAtIndex: 3];
				
				// Convert from a position to an index. index = (position - domainMin)/sampleInterval
				// D = A*b + c
				vGL_vsmsa( xPositions.mutableBytes, 1, (GLFloat *) &xDimInvInterval, (GLFloat *) &xDimOffset, firstInterpIndices.mutableBytes, 1, numInterpPoints);
				
				GLFloat *f = (GLFloat *) firstInterpIndices.mutableBytes;
				for (NSUInteger i=0; i<numInterpPoints; i++) {
					if (f[i]>nx) {
						f[i]=fmodf(f[i], nx);
					} else if (f[i]<0) {
						f[i]=nx+fmodf(f[i], nx);
					}
				}
				
				NSMutableData *buffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: firstInterpIndices.length];
				GLFloat one = 1.0;
				// secondInterpIndices = i+1
				vGL_vsadd( firstInterpIndices.mutableBytes, 1, (void *) &one, buffer.mutableBytes, 1, numInterpPoints);
				GLFloat *s = (GLFloat *) buffer.mutableBytes;
				for (NSUInteger i=0; i<numInterpPoints; i++) {
					if (s[i]>nx) {
						s[i]=fmodf(s[i], nx);
					} else if (s[i]<0) {
						s[i]=nx+fmodf(s[i], nx);
					}
				}
				
				
					vGL_vfrac( firstInterpIndices.mutableBytes, 1, xFrac.mutableBytes, 1, numInterpPoints);
					// C = B - A
					vGL_vsub( xFrac.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, numInterpPoints);
					vGL_vsub( xFrac.mutableBytes, 1, buffer.mutableBytes, 1, buffer.mutableBytes, 1, numInterpPoints);
					
					// xfrac = 1-xfrac
					vGL_vneg( xFrac.mutableBytes, 1, xFracOneMinus.mutableBytes, 1, numInterpPoints);
					vGL_vsadd( xFracOneMinus.mutableBytes, 1, &one, xFracOneMinus.mutableBytes, 1, numInterpPoints );
					
					
					
					
					// firstInterpIndices = ny*i
					vGL_vsmul( firstInterpIndices.mutableBytes, 1,  &ny, firstInterpIndices.mutableBytes, 1, numInterpPoints );
					vGL_vsmul( buffer.mutableBytes, 1,  &ny, buffer.mutableBytes, 1, numInterpPoints );
					
					
					// Go compute the y index value
					vGL_vsmsa( yPositions.mutableBytes, 1,(GLFloat *) &yDimInvInterval,(GLFloat *) &yDimOffset, secondInterpIndices.mutableBytes, 1, numInterpPoints);
					vGL_vclip( secondInterpIndices.mutableBytes, 1, (GLFloat *)&minIndex, (GLFloat *)&maxIndexY, secondInterpIndices.mutableBytes, 1, numInterpPoints);
					
					// firstInterpIndices = ny*i + j (including j fraction)
					vGL_vadd( firstInterpIndices.mutableBytes, 1, secondInterpIndices.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, numInterpPoints);
					// secondInterpIndices = ny*(i+1) + j (including j fraction)
					vGL_vadd( buffer.mutableBytes, 1, secondInterpIndices.mutableBytes, 1, secondInterpIndices.mutableBytes, 1, numInterpPoints);
					[[GLMemoryPool sharedMemoryPool] returnData: buffer];
			};
		}
		else if (allWrap)
		{
			maxIndexX = xDim.nPoints;
			maxIndexY = yDim.nPoints-1;
			self.blockOperation = ^(NSArray *result, NSArray *operand)
			{
				NSMutableData *xPositions = [operand objectAtIndex: 0];
				NSMutableData *yPositions = [operand objectAtIndex: 1];
				
				NSMutableData *firstInterpIndices = [result objectAtIndex: 0];
				NSMutableData *secondInterpIndices = [result objectAtIndex: 1];
				NSMutableData *xFrac = [result objectAtIndex: 2];
				NSMutableData *xFracOneMinus = [result objectAtIndex: 3];
				
				// Convert from a position to an index. index = (position - domainMin)/sampleInterval
				// D = A*b + c
				vGL_vsmsa( xPositions.mutableBytes, 1, (GLFloat *) &xDimInvInterval, (GLFloat *) &xDimOffset, firstInterpIndices.mutableBytes, 1, numInterpPoints);
				
				GLFloat *f = (GLFloat *) firstInterpIndices.mutableBytes;
				for (NSUInteger i=0; i<numInterpPoints; i++) {
					if (f[i]>nx) {
						f[i]=fmodf(f[i], nx);
					} else if (f[i]<0) {
						f[i]=nx+fmodf(f[i], nx);
					}
				}
				
				NSMutableData *buffer = [[GLMemoryPool sharedMemoryPool] dataWithLength: firstInterpIndices.length];
				GLFloat one = 1.0;
				// secondInterpIndices = i+1
				vGL_vsadd( firstInterpIndices.mutableBytes, 1, (void *) &one, buffer.mutableBytes, 1, numInterpPoints);
				GLFloat *s = (GLFloat *) buffer.mutableBytes;
				for (NSUInteger i=0; i<numInterpPoints; i++) {
					if (s[i]>nx) {
						s[i]=fmodf(s[i], nx);
					} else if (s[i]<0) {
						s[i]=nx+fmodf(s[i], nx);
					}
				}
				
				
				vGL_vfrac( firstInterpIndices.mutableBytes, 1, xFrac.mutableBytes, 1, numInterpPoints);
				// C = B - A
				vGL_vsub( xFrac.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, numInterpPoints);
				vGL_vsub( xFrac.mutableBytes, 1, buffer.mutableBytes, 1, buffer.mutableBytes, 1, numInterpPoints);
				
				// xfrac = 1-xfrac
				vGL_vneg( xFrac.mutableBytes, 1, xFracOneMinus.mutableBytes, 1, numInterpPoints);
				vGL_vsadd( xFracOneMinus.mutableBytes, 1, &one, xFracOneMinus.mutableBytes, 1, numInterpPoints );
				
				
				
				
				// firstInterpIndices = ny*i
				vGL_vsmul( firstInterpIndices.mutableBytes, 1,  &ny, firstInterpIndices.mutableBytes, 1, numInterpPoints );
				vGL_vsmul( buffer.mutableBytes, 1,  &ny, buffer.mutableBytes, 1, numInterpPoints );
				
				
				// Go compute the y index value
				vGL_vsmsa( yPositions.mutableBytes, 1,(GLFloat *) &yDimInvInterval,(GLFloat *) &yDimOffset, secondInterpIndices.mutableBytes, 1, numInterpPoints);
				
				f = (GLFloat *) secondInterpIndices.mutableBytes;
				for (NSUInteger i=0; i<numInterpPoints; i++) {
					if (f[i]>ny) {
						f[i]=fmodf(f[i], ny);
					} else if (f[i]<0) {
						f[i]=ny+fmodf(f[i], ny);
					}
				}
				
				// firstInterpIndices = ny*i + j (including j fraction)
				vGL_vadd( firstInterpIndices.mutableBytes, 1, secondInterpIndices.mutableBytes, 1, firstInterpIndices.mutableBytes, 1, numInterpPoints);
				// secondInterpIndices = ny*(i+1) + j (including j fraction)
				vGL_vadd( buffer.mutableBytes, 1, secondInterpIndices.mutableBytes, 1, secondInterpIndices.mutableBytes, 1, numInterpPoints);
				[[GLMemoryPool sharedMemoryPool] returnData: buffer];
			};

		}
		
	}
	return self;
}

- (unaryVectorOperation) blockOperation {
	return [super blockOperation];
}

@end

@interface GLLinearInterpolationOperation : GLBinaryOperation
@end

@implementation GLLinearInterpolationOperation

// The first operand is the function we're approximating, the second function is the (fractional) indices that we're approximating it at.
- (id) initWithFirstOperand:(GLVariable *)fOperand secondOperand:(GLVariable *)sOperand
{
	GLVariable *resultVariable = [GLVariable variableOfRealTypeWithDimensions: sOperand.dimensions forEquation: sOperand.equation];
	if (( self = [super init] )) {
		self.firstOperand = fOperand;
		self.secondOperand = sOperand;
        self.result = resultVariable;
        
        [self performSelector:@selector(setupDependencies)];
        
		NSUInteger numInterpPoints = sOperand.nDataPoints;
		NSUInteger numElements = fOperand.nDataElements;
		self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
			// C = A_b + alpha*(A_{b+1} - A_b)
			// b = trunc(B)
			// alpha = B-float(B)
			vGL_vlint( (GLFloat *) fOperand.bytes,(GLFloat *) sOperand.bytes, 1, result.mutableBytes, 1, numInterpPoints, numElements);
		};
		
//		self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
//			
//			GLFloat *f = (GLFloat *) fOperand.bytes;
//			GLFloat *ind = (GLFloat *) sOperand.bytes;
//			GLFloat *g = (GLFloat *) result.mutableBytes;
//			
//			for (NSUInteger i=0; i<numInterpPoints; i++)
//			{
//				NSUInteger idx = floor(ind[i]);
//				GLFloat frac = ind[i]-(GLFloat)idx;
//				
//				NSLog(@"Neighbors: (%f, %f)", f[idx], f[idx+1]);
//				
//				g[i] = f[idx]+frac*(f[idx+1]-f[idx]);
//				g[i] = (1-frac)*f[idx]+frac*f[idx+1];
//			}
//			
//		};
	}
	return self;
}

@end

@implementation GLInterpolationOperation

- (id) initWithFirstOperand: (NSArray *) fOperand secondOperand: (NSArray *) sOperand
{	
	for (GLVariable *fVariable in fOperand) {	
		if (fVariable.dimensions.count != sOperand.count) {
			[NSException raise: @"MethodNotImplemented" format: @"Position vector must have the same dimensionality as the functions being interpolated."];
		}
	}
	
	if (sOperand.count > 2) {
		[NSException raise: @"MethodNotImplemented" format: @"Interpolation is not yet implemented for three or more dimensions."];
	}
	
	for (GLVariable * fOperandVariable in fOperand) {
		if ( fOperandVariable.isComplex ) {
			[NSException raise: @"MethodNotImplemented" format: @"Interpolation can only be done at real points."];
		}
	}
	
	NSMutableArray *resultArray = [[NSMutableArray alloc] init];
	NSArray *dimensions = [[fOperand lastObject] dimensions];
	
	for (GLVariable *fVariable in fOperand) {
		if ( fVariable.isComplex ) {
			[NSException raise: @"MethodNotImplemented" format: @"Interpolation can only be done at real points."];
		} else if ( ![fVariable.dimensions isEqualToArray: dimensions] ) {
			[NSException raise: @"IllegalDimensions" format: @"All variables being interpolated must have the same dimensions."];
		} else {
			// The resulting variables will have the same length as the requested interpolation points.
			GLVariable *resultVariable = [GLVariable variableOfRealTypeWithDimensions:[[sOperand lastObject] dimensions] forEquation:[[sOperand lastObject] equation]];
			[resultArray addObject: resultVariable];
		}
	}
	
	if ( sOperand.count == 1)
	{
		GLLinearInterpolationIndexOperation *indexOperation = [[GLLinearInterpolationIndexOperation alloc] initWithOperand:sOperand.lastObject dimension: dimensions.lastObject];
		GLVariable *interpIndices = indexOperation.result;
		
		NSMutableArray *partialInterpolations = [[NSMutableArray alloc] init];
		for (GLVariable *fVariable in fOperand) {
			GLLinearInterpolationOperation *interp = [[GLLinearInterpolationOperation alloc] initWithFirstOperand: fVariable secondOperand: interpIndices];
			[partialInterpolations addObject: interp.result];
		}
		
		// This is a hollow operation. It does nothing.
		if (( self = [super init] ))
		{
			self.result = partialInterpolations;
			self.blockOperation = ^(NSArray *result, NSArray *fOperand, NSArray *sOperand) {};
		}
	}
	else if ( sOperand.count == 2)
	{
		GLBilinearInterpolationIndexOperation *indexOperation = [[GLBilinearInterpolationIndexOperation alloc] initWithOperand:sOperand dimensions: dimensions];
		GLVariable *firstInterpIndices = [indexOperation.result objectAtIndex: 0];
		GLVariable *secondInterpIndices = [indexOperation.result objectAtIndex: 1];
		GLVariable *xFrac = [indexOperation.result objectAtIndex: 2];
		GLVariable *xFracOneMinus = [indexOperation.result objectAtIndex: 3];
		
		NSMutableArray *partialInterpolations = [[NSMutableArray alloc] init];
		for (GLVariable *fVariable in fOperand) {
			GLLinearInterpolationOperation *firstInterp = [[GLLinearInterpolationOperation alloc] initWithFirstOperand: fVariable secondOperand: firstInterpIndices];
			GLLinearInterpolationOperation *secondInterp = [[GLLinearInterpolationOperation alloc] initWithFirstOperand: fVariable secondOperand: secondInterpIndices];
			
			[partialInterpolations addObject: firstInterp.result];
			[partialInterpolations addObject: secondInterp.result];
		}
		
		NSArray *fractions = [NSArray arrayWithObjects: xFrac, xFracOneMinus, nil];
		
		// What we *actually* do is initialize the second operand with the results of the index operation.
		if ( (self = [self initWithResult: resultArray firstOperand: partialInterpolations secondOperand: fractions]))
		{		
			dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
			NSUInteger numInterpPoints = firstInterpIndices.nDataPoints;
			NSUInteger numInterpFunctions = fOperand.count;
			
			self.blockOperation = ^(NSArray *result, NSArray *fOperand, NSArray *sOperand)
			{
				NSMutableData *xFracData = [sOperand objectAtIndex: 0];
				NSMutableData *xFracOneMinusData = [sOperand objectAtIndex: 1];
				
				dispatch_apply(numInterpFunctions, globalQueue, ^(size_t i){
					NSMutableData *firstInterpData = [fOperand objectAtIndex: 2*i];
					NSMutableData *secondInterpData = [fOperand objectAtIndex: 2*i+1];
					NSMutableData *resultData = [result objectAtIndex: i];
					
					// (1-frac)*firstInterp + frac*secondInterp
					vGL_vmma( xFracOneMinusData.mutableBytes, 1, firstInterpData.mutableBytes, 1, xFracData.mutableBytes, 1, secondInterpData.mutableBytes, 1, resultData.mutableBytes, 1, numInterpPoints);
				});
			};
		}
	}
	
	return self;
}

@end
