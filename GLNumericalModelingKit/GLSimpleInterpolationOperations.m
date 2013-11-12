//
//  GLSimpleInterpolationOperations.m
//  GLNumericalModelingKitCopy
//
//  Created by Jeffrey J. Early on 1/28/13.
//
//

#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import "GLSimpleInterpolationOperations.h"

/************************************************/
/*		GLSimpleLinearInterpolationIndexOperation				*/
/************************************************/

// Given 1) a vector of positions, 2) the dimension associated with those positions and 3) end point behavior
// this function returns four variables:
// 1) a vector of fractional indices corresponding to the nearest "floor"/"lower" index
// 2) a vector of fractional indices corresponding to the nearest "ceil"/"next" index
// 3) a vector, fraction, of the relative distance the position is between the above two indices. Value is between (0,1).
// 4) a vector with 1-fraction.
// This is designed as the input for the bilinear interpolation algorithm.

@interface GLDimensionalPositionOperation : GLVariableOperation
- (id) initWithPositionVector:(GLFunction *)positionVar dimension: (GLDimension *) dimension endPointPointBehavior: (GLInterpolationEndpointBehavior) behavior;
@property(strong) GLDimension *dimension;
@end

@implementation GLDimensionalPositionOperation

// This recursive function takes an array of monotonically increasing values, and locates the index of the value
// immediately below the argument value. The index returned will be between imin and imax-1.
NSInteger indexBelow2( GLFloat *monotonicallyIncreasingValues, GLFloat value, NSInteger imin, NSInteger imax )
{
	if (imin+1 == imax) {
		return imin;
	}
	else
	{
		NSInteger imid=imin+(imax-imin)/2;
		if (monotonicallyIncreasingValues[imid] > value) {
			return indexBelow2(monotonicallyIncreasingValues, value, imin, imid);
		}
		else {
			return indexBelow2(monotonicallyIncreasingValues, value, imid, imax);
		}
	}
}

- (id) initWithPositionVector:(GLFunction *) position dimension: (GLDimension *) dimension endPointPointBehavior: (GLInterpolationEndpointBehavior) behavior
{
	GLFunction *lowerIndicesVar = [GLFunction variableOfRealTypeWithDimensions:position.dimensions forEquation:position.equation];
	GLFunction *upperIndicesVar = [GLFunction variableOfRealTypeWithDimensions:position.dimensions forEquation:position.equation];
	GLFunction *fractionVar = [GLFunction variableOfRealTypeWithDimensions:position.dimensions forEquation:position.equation];
	GLFunction *oneMinusFractionVar = [GLFunction variableOfRealTypeWithDimensions:position.dimensions forEquation:position.equation];
	
	if ((self=[super initWithResult: @[lowerIndicesVar, upperIndicesVar, fractionVar, oneMinusFractionVar] operand: @[position]]))
	{
        self.dimension = dimension;
//		GLFloat minIndex = 0.0;
//		GLFloat maxIndexX = dimension.nPoints-1;
		GLFloat nx = dimension.nPoints;
        NSInteger n = dimension.nPoints;
		
		NSUInteger numInterpPoints = position.nDataPoints;
        
        // This stupid dimension looks like this: 0,..,fc-1,-fc,..,-1
        // But, this algorithm assume it's going like: -fc,..,-1,0,..,fc-1
        BOOL dimensionIsWrapped = (dimension.basisFunction == kGLExponentialBasis && dimension.isStrictlyPositive == NO);
		
		if (dimension.isEvenlySampled)
		{
			GLFloat xDimInvInterval = 1.0/dimension.sampleInterval;
			GLFloat xDimOffset = - dimension.domainMin / dimension.sampleInterval;
			
			if (behavior == kGLTruncationBehavior)
			{
				self.operation = ^(NSArray *result, NSArray *operands, NSArray *bufferArray) {
                    GLFloat *positions = (GLFloat *) [operands[0] bytes];
                    
                    GLFloat *lowerIndices = (GLFloat *) [result[0] mutableBytes];
                    GLFloat *upperIndices = (GLFloat *) [result[1] mutableBytes];
                    GLFloat *fraction = (GLFloat *) [result[2] mutableBytes];
                    GLFloat *oneMinusFraction = (GLFloat *) [result[3] mutableBytes];
					
                    // fraction contains the fractional, lower index.
					vGL_vsmsa( positions, 1, (GLFloat *) &xDimInvInterval, (GLFloat *) &xDimOffset, fraction, 1, numInterpPoints);
					
					for (NSUInteger i=0; i<numInterpPoints; i++) {
                        lowerIndices[i] = floor(fraction[i]);
                        fraction[i] = fraction[i] - lowerIndices[i];
                        oneMinusFraction[i] = 1.0 - fraction[i];
                        
                        NSInteger a = lowerIndices[i];
                        NSInteger b = a+1;
                        
                        lowerIndices[i] = a < 0 ? 0 : ( a >= n ? n-1 : a);
                        upperIndices[i] = b < 0 ? 0 : ( b >= n ? n-1 : b);
					}
                    
                    // Not efficient, but good enough for now.
                    if (dimensionIsWrapped) {
                        for (NSUInteger i=0; i<numInterpPoints; i++) {
                            if (lowerIndices[i] < nx/2) {
                                lowerIndices[i] += nx/2;
                            } else {
                                lowerIndices[i] -= nx/2;
                            }
                            
                            if (upperIndices[i] < nx/2) {
                                upperIndices[i] += nx/2;
                            } else {
                                upperIndices[i] -= nx/2;
                            }
                        }
                    }
				};
			}
			else if (behavior == kGLPeriodicBehavior)
			{
                // one minus frac weights the lower index, frac weights the upper index
                //
                // Say position = -0.25.
                // We want lower index to have nx-1, and upper index to have 0
                // We want one minus frac to have 0.25 and frac to have 0.75.
                //
                // Say position = 0.25.
                // We want lower index to have 0, and upper index to have 1
                // We want one minus frac to have 0.75 and frac to have 0.25.
                //
                // Say position = -1.25.
                // We want lower index to have nx-2, and upper index to have nx-1
                // We want one minus frac to have 0.25 and frac to have 0.75.
				self.operation = ^(NSArray *result, NSArray *operands, NSArray *bufferArray) {
                    GLFloat *positions = (GLFloat *) [operands[0] bytes];
                    
                    GLFloat *lowerIndices = (GLFloat *) [result[0] mutableBytes];
                    GLFloat *upperIndices = (GLFloat *) [result[1] mutableBytes];
                    GLFloat *fraction = (GLFloat *) [result[2] mutableBytes];
                    GLFloat *oneMinusFraction = (GLFloat *) [result[3] mutableBytes];
					
                    // fraction contains the fractional, lower index.
					vGL_vsmsa( positions, 1, (GLFloat *) &xDimInvInterval, (GLFloat *) &xDimOffset, fraction, 1, numInterpPoints);
					
					for (NSUInteger i=0; i<numInterpPoints; i++) {
                        lowerIndices[i] = floor(fraction[i]);
                        fraction[i] = fraction[i] - lowerIndices[i];
                        oneMinusFraction[i] = 1.0 - fraction[i];
                        
                        NSInteger a = lowerIndices[i];
                        NSInteger b = a+1;                        
                        
                        lowerIndices[i] = (((a < 0) ? ((a % n) + n) : a) % n);
                        upperIndices[i] = (((b < 0) ? ((b % n) + n) : b) % n);                      
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
			
			if (behavior == kGLTruncationBehavior)
			{
				self.operation = ^(NSArray *result, NSArray *operands, NSArray *bufferArray) {
					GLFloat *monotonicallyIncreasingValues = (GLFloat *) dimData.bytes;
                    GLFloat *positions = (GLFloat *) [operands[0] bytes];
                    
                    GLFloat *lowerIndices = (GLFloat *) [result[0] mutableBytes];
                    GLFloat *upperIndices = (GLFloat *) [result[1] mutableBytes];
                    GLFloat *fraction = (GLFloat *) [result[2] mutableBytes];
                    GLFloat *oneMinusFraction = (GLFloat *) [result[3] mutableBytes];
					
					dispatch_apply(numInterpPoints, globalQueue, ^(size_t i) {
						lowerIndices[i] = indexBelow2(monotonicallyIncreasingValues, positions[i], 0, nx);
                        
                        if (lowerIndices[i] == nx-1) {
                            upperIndices[i] = nx-1;
                            fraction[i] = 0.0;
                            oneMinusFraction[i] = 1.0;
                        } else {
                            upperIndices[i] = lowerIndices[i] + 1;
                            fraction[i] = (positions[i] - monotonicallyIncreasingValues[(NSUInteger)lowerIndices[i]])/(monotonicallyIncreasingValues[(NSUInteger)upperIndices[i]] - monotonicallyIncreasingValues[(NSUInteger)lowerIndices[i]]);
                            oneMinusFraction[i] = 1.0 - fraction[i];
                        }
					});
				};
			}
			else if (behavior == kGLPeriodicBehavior)
			{
				self.operation = ^(NSArray *result, NSArray *operands, NSArray *bufferArray) {
					GLFloat *monotonicallyIncreasingValues = (GLFloat *) dimData.bytes;
                    GLFloat *positions = (GLFloat *) [operands[0] bytes];
                    
                    GLFloat *lowerIndices = (GLFloat *) [result[0] mutableBytes];
                    GLFloat *upperIndices = (GLFloat *) [result[1] mutableBytes];
                    GLFloat *fraction = (GLFloat *) [result[2] mutableBytes];
                    GLFloat *oneMinusFraction = (GLFloat *) [result[3] mutableBytes];
					
					dispatch_apply(numInterpPoints, globalQueue, ^(size_t i) {
                        GLFloat modPosition;
                        if (positions[i]>domainMin+domainLength) {
							modPosition=domainMin+fmodf(positions[i]-domainMin, domainLength);
						} else if (positions[i]<domainMin) {
							modPosition=domainMin+domainLength+fmodf(positions[i]-domainMin, domainLength);
						} else {
                            modPosition = positions[i];
                        }
                        
						lowerIndices[i] = indexBelow2(monotonicallyIncreasingValues, modPosition, 0, nx);
                        
                        if (lowerIndices[i] == nx-1) {
                            upperIndices[i] = 0;
                        } else {
                            upperIndices[i] = lowerIndices[i] + 1;
                        }
                        
                        fraction[i] = (modPosition - monotonicallyIncreasingValues[(NSUInteger)lowerIndices[i]])/(monotonicallyIncreasingValues[(NSUInteger)upperIndices[i]] - monotonicallyIncreasingValues[(NSUInteger)lowerIndices[i]]);
                        oneMinusFraction[i] = 1.0 - fraction[i];
                        
					});
				};
			}
		}
	}
	return self;
}

- (BOOL) isEqualToOperation: (id) otherOperation
{
    if ( [super isEqualToOperation: otherOperation]) {
        if ( self.dimension == [(GLDimensionalPositionOperation *)otherOperation dimension] ) {
            return YES;
        }
    }
    
    return NO;
}

@end


@interface GLOneDimLinearInterpolationOperation : GLVariableOperation
@end

@implementation GLOneDimLinearInterpolationOperation

// The first operand is the function we're approximating, the second function is the (fractional) indices that we're approximating it at.
- (id) initWithFunction:(GLFunction *)fOperand lowerIndices: (GLFunction *) lIndices upperIndices: (GLFunction *) uIndices fraction: (GLFunction *) frac
{
	GLFunction *resultVariable = [GLFunction variableOfRealTypeWithDimensions: frac.dimensions forEquation: frac.equation];
	if (( self = [super initWithResult: @[resultVariable] operand: @[fOperand, lIndices, uIndices, frac]] )) {
		
        NSUInteger numInterpPoints = frac.nDataPoints;
		self.operation = ^(NSArray *result, NSArray *operands, NSArray *bufferArray) {
            GLFloat *f = (GLFloat *) [operands[0] bytes];
            GLFloat *lowerIndices = (GLFloat *) [operands[1] bytes];
            GLFloat *upperIndices = (GLFloat *) [operands[2] bytes];
            GLFloat *fraction = (GLFloat *) [operands[3] bytes];
            
            GLFloat *g = (GLFloat *) [result[0] mutableBytes];

			for (NSUInteger i=0; i<numInterpPoints; i++)
			{
				g[i] = (1.0 - fraction[i])*f[(NSUInteger)lowerIndices[i]] + fraction[i]*f[(NSUInteger)upperIndices[i]];
			}
			
		};
	}
	return self;
}
@end

@interface GLWeightedAverageOperation : GLVariableOperation
- (id) initWithLeftVariables: (NSArray *) leftVars rightVariables: (NSArray *) rightVars leftWeighting: (GLFunction *) leftWeight rightWeighting: (GLFunction *) rightWeight;
@end

@implementation GLWeightedAverageOperation

- (id) initWithLeftVariables: (NSArray *) leftVars rightVariables: (NSArray *) rightVars leftWeighting: (GLFunction *) leftWeight rightWeighting: (GLFunction *) rightWeight;
{
    
    NSMutableArray *resultArray = [[NSMutableArray alloc] init];
    NSMutableArray *inputArray = [[NSMutableArray alloc] init];
	for (NSUInteger i=0; i<leftVars.count; i++) {
        GLFunction *resultVariable = [GLFunction variableOfRealTypeWithDimensions: [leftVars[i] dimensions] forEquation: [leftVars[i] equation]];
        [resultArray addObject: resultVariable];
        inputArray[2*i] = leftVars[i];
        inputArray[2*i+1] = rightVars[i];
	}
    
	NSMutableArray *operandArray = [NSMutableArray arrayWithArray:@[leftWeight, rightWeight]];
	[operandArray addObjectsFromArray: inputArray];
	
	if ( (self = [self initWithResult: resultArray operand: operandArray]))
    {
        dispatch_queue_t globalQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
        NSUInteger numInterpPoints = leftWeight.nDataPoints;
        NSUInteger numInterpFunctions = leftVars.count;
        
		self.operation = ^(NSArray *result, NSArray *operands, NSArray *bufferArray)
        {
            NSMutableData *left = [operands objectAtIndex: 0];
            NSMutableData *right = [operands objectAtIndex: 1];
            
            dispatch_apply(numInterpFunctions, globalQueue, ^(size_t i){
                NSMutableData *firstInterpData = [operands objectAtIndex: 2*(i+1)];
                NSMutableData *secondInterpData = [operands objectAtIndex: 2*(i+1)+1];
                NSMutableData *resultData = [result objectAtIndex: i];
                
                // (1-frac)*firstInterp + frac*secondInterp
                vGL_vmma( left.mutableBytes, 1, firstInterpData.mutableBytes, 1, right.mutableBytes, 1, secondInterpData.mutableBytes, 1, resultData.mutableBytes, 1, numInterpPoints);
            });
        };
    }
    
    return self;
}

@end

@implementation GLSimpleInterpolationOperation

- (id) initWithFirstOperand: (NSArray *) fOperand secondOperand: (NSArray *) sOperand
{
	for (GLFunction *fVariable in fOperand) {
		if (fVariable.dimensions.count != sOperand.count) {
			[NSException raise: @"MethodNotImplemented" format: @"Position vector must have the same dimensionality as the functions being interpolated."];
		}
	}
	
	if (sOperand.count > 3) {
		[NSException raise: @"MethodNotImplemented" format: @"Interpolation is not yet implemented for four or more dimensions."];
	}
	
	NSMutableArray *resultArray = [[NSMutableArray alloc] init];
	NSArray *dimensions = [[fOperand lastObject] dimensions];
	
	for (GLFunction *fVariable in fOperand) {
		if ( fVariable.isComplex ) {
			[NSException raise: @"MethodNotImplemented" format: @"Interpolation can only be done at real points."];
		} else if ( ![fVariable.dimensions isEqualToArray: dimensions] ) {
			[NSException raise: @"IllegalDimensions" format: @"All variables being interpolated must have the same dimensions."];
		} else {
			// The resulting variables will have the same length as the requested interpolation points.
			GLFunction *resultVariable = [GLFunction variableOfRealTypeWithDimensions:[[sOperand lastObject] dimensions] forEquation:[[sOperand lastObject] equation]];
			[resultArray addObject: resultVariable];
		}
	}
	
	if ( sOperand.count == 1)
	{
		GLDimensionalPositionOperation *indexOperation = [[GLDimensionalPositionOperation alloc] initWithPositionVector: sOperand.lastObject dimension: dimensions.lastObject endPointPointBehavior: [dimensions.lastObject isPeriodic] ? kGLPeriodicBehavior : kGLTruncationBehavior];
		GLFunction *lowerIndices = indexOperation.result[0];
        GLFunction *upperIndices = indexOperation.result[1];
        GLFunction *fraction = indexOperation.result[2];
		
		NSMutableArray *partialInterpolations = [[NSMutableArray alloc] init];
		for (GLFunction *fVariable in fOperand) {
			GLOneDimLinearInterpolationOperation *interp = [[GLOneDimLinearInterpolationOperation alloc] initWithFunction: fVariable lowerIndices:lowerIndices upperIndices: upperIndices fraction: fraction];
			[partialInterpolations addObject: interp.result[0]];
		}
		
		// This is a hollow operation. It does nothing.
		if (( self = [super init] ))
		{
			self.result = partialInterpolations;
			self.operation = ^(NSArray *result, NSArray *fOperand, NSArray *sOperand) {};
		}
	}
	else if ( sOperand.count == 2)
	{
        GLMatrixDescription *matrixDescription = [fOperand[0] matrixDescription];
        GLDataStride *strides = matrixDescription.strides;
        
        GLDimensionalPositionOperation *indexOperation0 = [[GLDimensionalPositionOperation alloc] initWithPositionVector: sOperand[0] dimension: dimensions[0] endPointPointBehavior: [dimensions[0] isPeriodic] ? kGLPeriodicBehavior : kGLTruncationBehavior];
		GLFunction *lowerIndices0 = indexOperation0.result[0];
        GLFunction *upperIndices0 = indexOperation0.result[1];
        GLFunction *fraction0 = indexOperation0.result[2];
        GLFunction *oneMinusfraction0 = indexOperation0.result[3];
        
        GLDimensionalPositionOperation *indexOperation1 = [[GLDimensionalPositionOperation alloc] initWithPositionVector: sOperand[1] dimension: dimensions[1] endPointPointBehavior: [dimensions[1] isPeriodic] ? kGLPeriodicBehavior : kGLTruncationBehavior];
		GLFunction *lowerIndices1 = indexOperation1.result[0];
        GLFunction *upperIndices1 = indexOperation1.result[1];
        GLFunction *fraction1 = indexOperation1.result[2];
        
        GLFunction *lowerLower = [[lowerIndices0 times: @(strides[0].stride)] plus: lowerIndices1];
        GLFunction *lowerUpper = [[lowerIndices0 times: @(strides[0].stride)] plus: upperIndices1];
        
        GLFunction *upperLower = [[upperIndices0 times: @(strides[0].stride)] plus: lowerIndices1];
        GLFunction *upperUpper = [[upperIndices0 times: @(strides[0].stride)] plus: upperIndices1];
        
		NSMutableArray *leftVariables = [[NSMutableArray alloc] init];
        NSMutableArray *rightVariables = [[NSMutableArray alloc] init];
		for (GLFunction *fVariable in fOperand) {
			GLOneDimLinearInterpolationOperation *firstInterp = [[GLOneDimLinearInterpolationOperation alloc] initWithFunction: fVariable lowerIndices:lowerLower upperIndices: lowerUpper fraction: fraction1];
			GLOneDimLinearInterpolationOperation *secondInterp = [[GLOneDimLinearInterpolationOperation alloc] initWithFunction: fVariable lowerIndices:upperLower upperIndices: upperUpper fraction: fraction1];
			
			[leftVariables addObject: firstInterp.result[0]];
			[rightVariables addObject: secondInterp.result[0]];
		}
        
        GLWeightedAverageOperation *avgOp = [[GLWeightedAverageOperation alloc] initWithLeftVariables: leftVariables rightVariables: rightVariables leftWeighting: oneMinusfraction0 rightWeighting:fraction0];
				
		// This is a hollow operation. It does nothing.
		if (( self = [super init] ))
		{
			self.result = avgOp.result;
			self.operation = ^(NSArray *result, NSArray *fOperand, NSArray *sOperand) {};
		}
	}
    else if ( sOperand.count == 3)
	{
        GLMatrixDescription *matrixDescription = [fOperand[0] matrixDescription];
        GLDataStride *strides = matrixDescription.strides;
        
        GLDimensionalPositionOperation *indexOperation0 = [[GLDimensionalPositionOperation alloc] initWithPositionVector: sOperand[0] dimension: dimensions[0] endPointPointBehavior: [dimensions[0] isPeriodic] ? kGLPeriodicBehavior : kGLTruncationBehavior];
		GLFunction *lowerIndices0 = indexOperation0.result[0];
        GLFunction *upperIndices0 = indexOperation0.result[1];
        GLFunction *fraction0 = indexOperation0.result[2];
        GLFunction *oneMinusfraction0 = indexOperation0.result[3];
        
        GLDimensionalPositionOperation *indexOperation1 = [[GLDimensionalPositionOperation alloc] initWithPositionVector: sOperand[1] dimension: dimensions[1] endPointPointBehavior: [dimensions[1] isPeriodic] ? kGLPeriodicBehavior : kGLTruncationBehavior];
		GLFunction *lowerIndices1 = indexOperation1.result[0];
        GLFunction *upperIndices1 = indexOperation1.result[1];
        GLFunction *fraction1 = indexOperation1.result[2];
        GLFunction *oneMinusfraction1 = indexOperation1.result[3];
        
        GLDimensionalPositionOperation *indexOperation2 = [[GLDimensionalPositionOperation alloc] initWithPositionVector: sOperand[2] dimension: dimensions[2] endPointPointBehavior: [dimensions[2] isPeriodic] ? kGLPeriodicBehavior : kGLTruncationBehavior];
		GLFunction *lowerIndices2 = indexOperation2.result[0];
        GLFunction *upperIndices2 = indexOperation2.result[1];
        GLFunction *fraction2 = indexOperation2.result[2];
        
        GLFunction *lowerLowerLower = [[[[lowerIndices0 scalarMultiply: strides[1].nPoints] plus: lowerIndices1] scalarMultiply: strides[2].nPoints] plus: lowerIndices2];
        GLFunction *lowerLowerUpper = [[[[lowerIndices0 scalarMultiply: strides[1].nPoints] plus: lowerIndices1] scalarMultiply: strides[2].nPoints] plus: upperIndices2];
        
        GLFunction *lowerUpperLower = [[[[lowerIndices0 scalarMultiply: strides[1].nPoints] plus: upperIndices1] scalarMultiply: strides[2].nPoints] plus: lowerIndices2];
        GLFunction *lowerUpperUpper = [[[[lowerIndices0 scalarMultiply: strides[1].nPoints] plus: upperIndices1] scalarMultiply: strides[2].nPoints] plus: upperIndices2];
        
        GLFunction *upperLowerLower = [[[[upperIndices0 scalarMultiply: strides[1].nPoints] plus: lowerIndices1] scalarMultiply: strides[2].nPoints] plus: lowerIndices2];
        GLFunction *upperLowerUpper = [[[[upperIndices0 scalarMultiply: strides[1].nPoints] plus: lowerIndices1] scalarMultiply: strides[2].nPoints] plus: upperIndices2];
        
        GLFunction *upperUpperLower = [[[[upperIndices0 scalarMultiply: strides[1].nPoints] plus: upperIndices1] scalarMultiply: strides[2].nPoints] plus: lowerIndices2];
        GLFunction *upperUpperUpper = [[[[upperIndices0 scalarMultiply: strides[1].nPoints] plus: upperIndices1] scalarMultiply: strides[2].nPoints] plus: upperIndices2];
        
		NSMutableArray *leftVariables = [[NSMutableArray alloc] init];
        NSMutableArray *rightVariables = [[NSMutableArray alloc] init];
		for (GLFunction *fVariable in fOperand) {
			GLOneDimLinearInterpolationOperation *firstInterp = [[GLOneDimLinearInterpolationOperation alloc] initWithFunction: fVariable lowerIndices:lowerLowerLower upperIndices: lowerLowerUpper fraction: fraction2];
			GLOneDimLinearInterpolationOperation *secondInterp = [[GLOneDimLinearInterpolationOperation alloc] initWithFunction: fVariable lowerIndices:lowerUpperLower upperIndices: lowerUpperUpper fraction: fraction2];
			
			[leftVariables addObject: firstInterp.result[0]];
			[rightVariables addObject: secondInterp.result[0]];
            
             firstInterp = [[GLOneDimLinearInterpolationOperation alloc] initWithFunction: fVariable lowerIndices:upperLowerLower upperIndices: upperLowerUpper fraction: fraction2];
			 secondInterp = [[GLOneDimLinearInterpolationOperation alloc] initWithFunction: fVariable lowerIndices:upperUpperLower upperIndices: upperUpperUpper fraction: fraction2];
			
			[leftVariables addObject: firstInterp.result[0]];
			[rightVariables addObject: secondInterp.result[0]];
		}
        
        // The results should contain half the number of variables now.
        GLWeightedAverageOperation *avgOp = [[GLWeightedAverageOperation alloc] initWithLeftVariables: leftVariables rightVariables: rightVariables leftWeighting: oneMinusfraction1 rightWeighting:fraction1];
        leftVariables = [[NSMutableArray alloc] init];
        rightVariables = [[NSMutableArray alloc] init];
        for (NSUInteger i=0; i<fOperand.count; i++) {
            leftVariables[i] = avgOp.result[2*i];
            rightVariables[i] = avgOp.result[2*i+1];
        }
        
        avgOp = [[GLWeightedAverageOperation alloc] initWithLeftVariables: leftVariables rightVariables: rightVariables leftWeighting: oneMinusfraction0 rightWeighting:fraction0];
                
		// This is a hollow operation. It does nothing.
		if (( self = [super init] ))
		{
			self.result = avgOp.result;
			self.operation = ^(NSArray *result, NSArray *fOperand, NSArray *sOperand) {};
		}
	}
	
	return self;
}


@end
