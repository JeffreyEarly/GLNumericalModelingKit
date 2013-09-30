//
//  GLDimensionalOperations.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey Early on 1/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#import "GLDimensionalOperations.h"
#import "GLDimension.h"

/************************************************/
/*		GLAddDimensionOperation                 */
/************************************************/

@implementation GLAddDimensionOperation

- (id) initWithOperand: (GLVariable *) variable dimension: (GLDimension *) dim
{
	if ([variable.dimensions containsObject: dim]) {
		return nil;
	}
	
	NSMutableArray *newDimensions = [NSMutableArray arrayWithObject: dim];
	[newDimensions addObjectsFromArray: variable.dimensions];
	
	GLVariable *resultVariable = [GLVariable variableOfType: variable.dataFormat withDimensions: newDimensions forEquation: variable.equation];
	resultVariable.name = variable.name;
	resultVariable.units = variable.units;
	
	if (( self = [super initWithResult: resultVariable operand: variable] ))
	{		
		self.theDimension = dim;
		self.nPoints = dim.nPoints;
		
		NSUInteger n = self.operand.nDataPoints;
		NSUInteger numPoints = dim.nPoints;
		
		if ( self.result.isComplex )
		{
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLSplitComplex toSplit = splitComplexFromData(result);
				GLSplitComplex fromSplit = splitComplexFromData(operand);
				
				// We extended by making copies of the data. That's the only sensible thing to do.
				for (NSUInteger i=0; i < numPoints; i++)
				{
					memcpy( &(toSplit.imagp[i*n]), fromSplit.imagp, n*sizeof(GLFloat));
					memcpy( &(toSplit.realp[i*n]), fromSplit.realp, n*sizeof(GLFloat));
				}
			};
		}
		else {
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLFloat * toData = result.mutableBytes;
				const GLFloat * fromData = operand.bytes;
				
				// We extended by making copies of the data. That's the only sensible thing to do.
				for (NSUInteger i=0; i < numPoints; i++)
				{
					memcpy( &toData[i*n], fromData, n*sizeof(GLFloat));
				}
			};
		}
	}
	
    return self;
}

@synthesize theDimension;
@synthesize nPoints;

@end

/************************************************/
/*		GLSubdomainOperation					*/
/************************************************/

// Returns a variable with only the elements indicated by by the array of ranges.
// The size of the ranges array must match the number of dimensions.
// If the -shouldFlatten flag is set, this will cause any dimensions of length 1 to be eliminated.
@implementation GLSubdomainOperation

- (id) initWithOperand: (GLVariable *) variable indexRange: (NSArray *) ranges flatten: (BOOL) aFlag
{
	self.shouldFlatten = aFlag;
	
	// Sanity check first
	if ( variable.dimensions.count != ranges.count) return nil;
	
	// What do the reduced dimensions look like?
	NSMutableArray *newDimensions = [[NSMutableArray alloc] init];
	BOOL invalidRange = NO;
	for ( GLDimension *dim in variable.dimensions ) {
		NSRange range = [[ranges objectAtIndex: [variable.dimensions indexOfObject: dim]] rangeValue];
		GLDimension *aNewDimension = [dim subdimensionWithRange: range];
		if (aNewDimension) {
			if ( !(self.shouldFlatten && aNewDimension.nPoints == 1) ) {
				[newDimensions addObject: aNewDimension];
			}
		} else {
			invalidRange = YES;
			NSLog(@"GLSubdomain Operation reports index out of bounds.");
		}
	}
	if (invalidRange) return nil;
	
	GLVariable *resultVariable = [GLVariable variableOfType: variable.dataFormat withDimensions: newDimensions forEquation: variable.equation];
	
	if (( self = [super initWithResult:resultVariable operand:variable] ))
	{
		self.theRanges = ranges;
		NSArray *dimensions = self.operand.dimensions;
        NSMutableString *rangeDescrip = [[NSMutableString alloc] initWithString: @"("];
        for (NSValue *aValue in self.theRanges) {
            [rangeDescrip appendFormat: @"%ld:%ld,", (unsigned long)aValue.rangeValue.location, aValue.rangeValue.location+aValue.rangeValue.length-1];
        }
        [rangeDescrip deleteCharactersInRange: NSMakeRange(rangeDescrip.length-1, 1)];
        [rangeDescrip appendFormat: @")"];
		
		if ( self.result.isComplex )
		{
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLSplitComplex toSplit = splitComplexFromData(result);
				GLSplitComplex fromSplit = splitComplexFromData(operand);
				
				CopySubmatrix( dimensions, ranges, 0, 0, fromSplit.realp, toSplit.realp);
				CopySubmatrix( dimensions, ranges, 0, 0, fromSplit.imagp, toSplit.imagp);
			};
            self.graphvisDescription = [NSString stringWithFormat: @"complex subdomain %@", rangeDescrip];
		}
		else {
			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
				GLFloat * toData = result.mutableBytes;
				const GLFloat * fromData = operand.bytes;
				
				CopySubmatrix( dimensions, ranges, 0, 0, fromData, toData);
			};
            self.graphvisDescription = [NSString stringWithFormat: @"real subdomain %@", rangeDescrip];
		}
	}
	
    return self;
}

@synthesize theRanges;
@synthesize shouldFlatten;

void CopySubmatrix( NSArray *dimensions, NSArray *ranges, NSUInteger matrixIndex, NSUInteger submatrixIndex, const GLFloat *fromData, GLFloat *toData );

// The fromData matrix must be the product of all dimensions.nPoints
// The toData matrix must be the product of all ranges.lengths
// The three functions below this are explict examples of what this is optimizing.

// The "matrix index" is the starting index in the array---an initial offset that just needs to respected from some higher dimensional variables.
void CopySubmatrix( NSArray *dimensions, NSArray *ranges, NSUInteger matrixIndex, NSUInteger submatrixIndex, const GLFloat *fromData, GLFloat *toData )
{
	if ( dimensions.count == 0 )
	{
		return;
	}
	else
	{
		NSRange range = [ranges[0] rangeValue];
		GLDimension *dim = dimensions[0];
		
		matrixIndex = matrixIndex * dim.nPoints + range.location;
		submatrixIndex = submatrixIndex * range.length;
		
		if ( dimensions.count == 1 )
		{
			memcpy( &toData[submatrixIndex], &fromData[matrixIndex], range.length*sizeof(GLFloat));
		}
		else if ( dimensions.count == 2 )
		{
			NSRange fastRange = [ranges[1] rangeValue];
			GLDimension *fastDim = dimensions[1];
            
            matrixIndex = matrixIndex * fastDim.nPoints + fastRange.location;
            submatrixIndex = submatrixIndex * fastRange.length;
			
			vGL_mmov( (void *) &(fromData[matrixIndex]), &(toData[submatrixIndex]), fastRange.length, range.length, fastDim.nPoints, fastRange.length);
		}
		else
		{
			NSArray *poppedDimensions = [dimensions subarrayWithRange: NSMakeRange(1, dimensions.count-1)];
			NSArray *poppedRanges = [ranges subarrayWithRange: NSMakeRange(1, ranges.count-1)];
			
			for (NSUInteger i=0; i < range.length; i++)
			{
				NSUInteger newMatrixIndex = matrixIndex + i;
				NSUInteger newSubmatrixIndex = submatrixIndex + i;
				
				CopySubmatrix( poppedDimensions, poppedRanges, newMatrixIndex, newSubmatrixIndex, fromData, toData);
			}
		}
	}
}

// vec[(i*ny+j)*nz+k] = matrix[i][j][k]
- (void) threeLoop3DCopy
{
	GLFloat *fromData = self.operand.pointerValue;
	GLFloat *toData = self.result.pointerValue;
	
	NSUInteger ny = [[self.operand.dimensions objectAtIndex: 1] nPoints];
	NSUInteger nz = [[self.operand.dimensions objectAtIndex: 2] nPoints];
	
	NSRange xRange = [[self.theRanges objectAtIndex: 0] rangeValue];
	NSRange yRange = [[self.theRanges objectAtIndex: 1] rangeValue];
	NSRange zRange = [[self.theRanges objectAtIndex: 2] rangeValue];
	
	for (NSUInteger i=0; i < xRange.length; i++)
	{
		for (NSUInteger j=0; j < yRange.length; j++)
		{
			for (NSUInteger k=0; k < zRange.length; k++)
			{
				// When k=zRange.location, this is the start of a contiguous block of memory to be copied, of length zRange.length
				NSUInteger submatrixIndex = (i*yRange.length+j)*zRange.length+k;
				NSUInteger matrixIndex = ((i+xRange.location)*ny+(j+yRange.location))*nz+(k+zRange.location);
				toData[submatrixIndex] = fromData[matrixIndex];
			}
		}
	}
}

- (void) twoLoop3DCopy
{
	GLFloat *fromData = self.operand.pointerValue;
	GLFloat *toData = self.result.pointerValue;
	
	NSUInteger ny = [[self.operand.dimensions objectAtIndex: 1] nPoints];
	NSUInteger nz = [[self.operand.dimensions objectAtIndex: 2] nPoints];
	
	NSRange xRange = [[self.theRanges objectAtIndex: 0] rangeValue];
	NSRange yRange = [[self.theRanges objectAtIndex: 1] rangeValue];
	NSRange zRange = [[self.theRanges objectAtIndex: 2] rangeValue];
	
	for (NSUInteger i=0; i < xRange.length; i++)
	{
		for (NSUInteger j=0; j < yRange.length; j++)
		{
			NSUInteger submatrixIndex = (i*yRange.length+j)*zRange.length;
			NSUInteger matrixIndex = ((i+xRange.location)*ny+(j+yRange.location))*nz+zRange.location;
			memcpy( &(toData[submatrixIndex]), &(fromData[matrixIndex]), zRange.length);
		}
	}
}

- (void) oneLoop3DCopy
{
	GLFloat *fromData = self.operand.pointerValue;
	GLFloat *toData = self.result.pointerValue;
	
	NSUInteger ny = [[self.operand.dimensions objectAtIndex: 1] nPoints];
	NSUInteger nz = [[self.operand.dimensions objectAtIndex: 2] nPoints];
	
	NSRange xRange = [[self.theRanges objectAtIndex: 0] rangeValue];
	NSRange yRange = [[self.theRanges objectAtIndex: 1] rangeValue];
	NSRange zRange = [[self.theRanges objectAtIndex: 2] rangeValue];
	
	for (NSUInteger i=0; i < xRange.length; i++)
	{
		NSUInteger submatrixIndex = i*yRange.length;
		NSUInteger matrixIndex = (i+xRange.location)*ny + yRange.location;
		
		// The number of columns, is the length of the fastest index
		// The number of rows, is the length the slowest index
		vGL_mmov( &(fromData[matrixIndex]), &(toData[submatrixIndex]), zRange.length, yRange.length, nz, zRange.length);
	}
}

@end



/************************************************/
/*		GLConcatenationOperation                 */
/************************************************/

@implementation GLExistingDimensionConcatenationOperation

// If we concatenate on index 0, then...
// An array of size mxn concatenated with an array of size pxn produces an array of sized (m+p)xn
// An array of size mxn concatenated with an array of size n produces an array of sized (m+1)xn

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand dimensionIndex: (NSUInteger) dimIndex
{
	NSMutableArray *resultDimensions = [NSMutableArray arrayWithArray: fOperand.dimensions];
	
	if ( fOperand.dimensions.count == sOperand.dimensions.count)
	{
		for (NSUInteger i=0; i<fOperand.dimensions.count; i++)
		{
			GLDimension *leftDim = [fOperand.dimensions objectAtIndex: i];
			GLDimension *rightDim = [sOperand.dimensions objectAtIndex: i];
			
			if (i==dimIndex)
			{	// We only require that the left operand have an evenly sampled dimension.
				if (!leftDim.isEvenlySampled) {
					[NSException raise: @"DimensionsException" format: @"The dimension of concatenation must be evenly sampled."];
				}
				GLDimension *newDimension = [[GLDimension alloc] initPeriodicDimension:leftDim.isPeriodic nPoints:leftDim.nPoints+rightDim.nPoints domainMin:leftDim.domainMin sampleInterval:leftDim.sampleInterval];
				[resultDimensions replaceObjectAtIndex: dimIndex withObject: newDimension];
			}
			else
			{	// We only require the number of points match.
				if (leftDim.nPoints != rightDim.nPoints) {
					[NSException raise: @"DimensionsException" format: @"The dimensional lengths do not match."];
				}
			}
		}	
	}
	else if ( fOperand.dimensions.count == sOperand.dimensions.count+1)
	{
		NSMutableArray *reducedDimensions = [NSMutableArray arrayWithArray: fOperand.dimensions];
		GLDimension *oldDimension = [reducedDimensions objectAtIndex: dimIndex];
		if (!oldDimension.isEvenlySampled) {
			[NSException raise: @"DimensionsException" format: @"The dimension of concatenation must be evenly sampled."];
		}
		GLDimension *newDimension = [[GLDimension alloc] initPeriodicDimension:oldDimension.isPeriodic nPoints:oldDimension.nPoints+1 domainMin:oldDimension.domainMin sampleInterval:oldDimension.sampleInterval];
		[resultDimensions replaceObjectAtIndex: dimIndex withObject: newDimension];
		[reducedDimensions removeObjectAtIndex: dimIndex];
		
		for (NSUInteger i=0; i<reducedDimensions.count; i++)
		{
			GLDimension *leftDim = [reducedDimensions objectAtIndex: i];
			GLDimension *rightDim = [sOperand.dimensions objectAtIndex: i];
			
			// We only require the number of points match.
			if (leftDim.nPoints != rightDim.nPoints) {
				[NSException raise: @"DimensionsException" format: @"The dimensional lengths do not match."];
			}
		}
	}
	else {
		[NSException raise: @"DimensionsException" format: @"Incompatible number of dimensions for concatenation."];
		return nil;
	}
	
	if ( fOperand.isComplex != sOperand.isComplex ) {
		[NSException raise: @"MethodNotYetImplemented" format: @"Cannot perform a binary operation on two variables that are not either both complex or both real."]; return nil;
	}
	
	GLVariable *resultVariable = [GLVariable variableOfType: fOperand.dataFormat withDimensions: resultDimensions forEquation:fOperand.equation];
	if (( self = [super initWithResult: resultVariable firstOperand:fOperand secondOperand:sOperand] )) {
		self.blockOperation = ^(NSMutableData *result, NSData *operand1, NSData *operand2) {
			NSLog(@"Oh crap! You just called a function (GLExistingDimensionConcatenationOperation) that isn't yet implemented!");
		};
	}
    
    return self;
}

@synthesize theDimension;

- (void) main
{	
	NSLog(@"Oh crap! You just called a function (GLExistingDimensionConcatenationOperation) that isn't yet implemented!");
//	NSUInteger n = self.firstOperand.nDataPoints;
//	size_t nBytes = self.firstOperand.nDataPoints*sizeof(GLFloat);
//	
//	if ( self.result.isComplex )
//	{
//		GLSplitComplex firstData = self.firstOperand.splitComplex;
//		GLSplitComplex secondData = self.secondOperand.splitComplex;
//		GLSplitComplex toData = self.result.splitComplex;
//		
//		memcpy( toData.imagp, firstData.imagp, nBytes);
//		memcpy( toData.realp, firstData.realp, nBytes);
//		memcpy( &(toData.imagp[n]), secondData.imagp, nBytes);
//		memcpy( &(toData.realp[n]), secondData.realp, nBytes);
//	}
//	else
//	{
//		// If the dimension has 0 points, this should be alright.
//		GLFloat *firstData = self.firstOperand.pointerValue;
//		GLFloat *secondData = self.secondOperand.pointerValue;
//		GLFloat *toData = self.result.pointerValue;
//		
//		if (dimIndex == 0)
//		{
//			memcpy( toData, firstData, nBytes);
//			memcpy( &(toData[n]), secondData, nBytes);
//		}
//	}
}

//// vec[(i*ny+j)*nz+k] = matrix[i][j][k]
//- (void) threeLoop3DCopy
//{
//	GLFloat *firstData = self.firstOperand.pointerValue;
//	GLFloat *secondData = self.secondOperand.pointerValue;
//	GLFloat *toData = self.result.pointerValue;
//	
//	NSUInteger nx = [[self.result.dimensions objectAtIndex: 0] nPoints];
//	NSUInteger ny = [[self.result.dimensions objectAtIndex: 1] nPoints];
//	NSUInteger nz = [[self.result.dimensions objectAtIndex: 2] nPoints];
//	
//	NSUInteger catIndex = 0;
//	NSUInteger maxFirst = [[self.firstOperand.dimensions objectAtIndex: catIndex] nPoints];
//	NSUInteger maxSecond = [[self.secondOperand.dimensions objectAtIndex: catIndex] nPoints];
//	
//	for (NSUInteger i=0; i < nx; i++)
//	{
//		for (NSUInteger j=0; j < ny; j++)
//		{
//			for (NSUInteger k=0; k < nz; k++)
//			{
//				// When k=zRange.location, this is the start of a contiguous block of memory to be copied, of length zRange.length
//				NSUInteger submatrixIndex = (i*ny+j)*nz+k;
//				NSUInteger matrixIndex = ((i+xRange.location)*ny+(j+yRange.location))*nz+(k+zRange.location);
//				toData[submatrixIndex] = fromData[matrixIndex];
//			}
//		}
//	}
//}

@end

@implementation GLNewDimensionConcatenationOperation

- (id) initWithFirstOperand: (GLVariable *) fOperand secondOperand: (GLVariable *) sOperand dimension: (GLDimension *) dim
{	
	if ( ![fOperand.dimensions isEqualToArray: sOperand.dimensions] ) {
        [NSException raise: @"DimensionsNotEqualException" format: @"Cannot add two variables of different dimensions"]; return nil;
    }
	if ( fOperand.isComplex != sOperand.isComplex ) {
		[NSException raise: @"MethodNotYetImplemented" format: @"Cannot perform a binary operation on two variables that are not either both complex or both real."]; return nil;
	}
    
	GLDimension *theDimension = dim.nPoints==2 ? dim : [[GLDimension alloc] initPeriodicDimension: dim.isPeriodic nPoints: 2 domainMin: dim.domainMin sampleInterval:dim.sampleInterval];
	NSMutableArray *resultDimensions = [NSMutableArray arrayWithObject: theDimension];
	[resultDimensions addObjectsFromArray: fOperand.dimensions];
	
	GLVariable *resultVariable = [GLVariable variableOfType: fOperand.dataFormat withDimensions: resultDimensions forEquation:fOperand.equation];
	
	if (( self = [super initWithResult:resultVariable firstOperand:fOperand secondOperand:sOperand] )) {
		
		NSUInteger n = self.firstOperand.nDataPoints;
		size_t nBytes = self.firstOperand.nDataPoints*sizeof(GLFloat);
		
		if ( self.result.isComplex )
		{
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				GLSplitComplex firstData = splitComplexFromData( fOperand );
				GLSplitComplex secondData = splitComplexFromData( sOperand );
				GLSplitComplex toData = splitComplexFromData( result );
				
				memcpy( toData.imagp, firstData.imagp, nBytes);
				memcpy( toData.realp, firstData.realp, nBytes);
				memcpy( &(toData.imagp[n]), secondData.imagp, nBytes);
				memcpy( &(toData.realp[n]), secondData.realp, nBytes);
			};
		}
		else {
			self.blockOperation = ^(NSMutableData *result, NSData *fOperand, NSData *sOperand) {
				// If the dimension has 0 points, this should be alright.
				const GLFloat *firstData = fOperand.bytes;
				const GLFloat *secondData = sOperand.bytes;
				GLFloat *toData = result.mutableBytes;
				
				memcpy( toData, firstData, nBytes);
				memcpy( &(toData[n]), secondData, nBytes);
			};
		}

	}
    
    return self;
}

@end

/************************************************/
/*		GLHalfToSplitComplexOperation			*/
/************************************************/

@implementation GLHalfToSplitComplexOperation


//- (id) initWithOperand: (GLVariable *) variable
//{
//	if (variable.dataFormat != kGLHalfComplexDataFormat) {
//		return nil;
//	}
//	
//	GLVariable *resultVariable = [GLVariable variableOfType: kGLSplitComplexDataFormat withDimensions: variable.dimensions forEquation:variable.equation];
//	if (( self = [super initWithResult: resultVariable operand: variable] ))
//	{
//		NSUInteger nDataElements = resultVariable.nDataElements;
//		
//		// We create a fake set of dimensions where anything in half complex format is treated as if all of its elements are points.
//		NSMutableArray *dummyDimensions = [NSMutableArray array];
//		NSMutableArray *realRanges = [NSMutableArray array];
//		NSMutableArray *imaginaryRanges = [NSMutableArray array];
//		for (GLDimension *aDim in variable.dimensions) {
//			GLDataFormat dataFormat = [[variable.dataFormats objectAtIndex: [variable.dimensions indexOfObject: aDim]] unsignedIntegerValue];
//			if (dataFormat == kGLHalfComplexDataFormat) {
//				[dummyDimensions addObject: [[GLDimension alloc] initPeriodicDimension: NO nPoints: 2*(aDim.nPoints-1) domainMin: 0 length: 1.0]];
//				[realRanges addObject: [NSValue valueWithRange: NSMakeRange(0, aDim.nPoints)]];
//				[imaginaryRanges addObject: [NSValue valueWithRange: NSMakeRange(aDim.nPoints, aDim.nPoints-2)]];
//			} else {
//				[dummyDimensions addObject: aDim];
//				[realRanges addObject: [NSValue valueWithRange: NSMakeRange(0, aDim.nPoints)]];
//				[imaginaryRanges addObject: [NSValue valueWithRange: NSMakeRange(0, aDim.nPoints)]];
//			}
//		}
//		
//		if (resultVariable.dimensions.count == 1)
//		{
//			NSUInteger nPoints = [variable.dimensions.lastObject nPoints];
//			self.blockOperation = ^(NSMutableData *result, NSData *operand) {
//				
//				GLSplitComplex toSplit = splitComplexFromData(result);
//				const GLFloat * fromData = operand.bytes;
//				
//				memcpy(toSplit.realp, fromData, nPoints*sizeof(GLFloat));
//				toSplit.imagp[0] = 0;
//				memcpy(&(toSplit.imagp[1]), &(fromData[nPoints]), (nPoints-2)*sizeof(GLFloat));
//				vGL_vrvrs(&(toSplit.imagp[1]), 1, nPoints-2);
//				toSplit.imagp[nPoints-1] = 0;
//			};
//		}
//		
////		self.blockOperation = ^(NSMutableData *result, NSData *operand) {
////			
////			GLSplitComplex toSplit = splitComplexFromData(result);
////			const GLFloat * fromData = operand.bytes;
////			
////			vGL_vclr( result.mutableBytes, 1, nDataElements);
////			CopySubmatrix( dummyDimensions, realRanges, 0, 0, fromData, toSplit.realp);
////			CopySubmatrix( dummyDimensions, realRanges, 0, 0, fromData, toSplit.imagp);
////			
////		};
//	}
//	return self;
//}

@end

/************************************************/
/*		GLZeroPadOperation                      */
/************************************************/

@implementation GLZeroPadOperation

- (id) initWithOperand: (GLVariable *) variable newDimensions: (NSArray *) newDimensions basis: (NSArray *) basis
{
    if (variable.dimensions.count != newDimensions.count) {
        [NSException raise: @"DimensionsException" format: @"You must have the same number of dimensions."];
    }
    
    GLVariable *lowResTransformedVariable = [variable transformToBasis: basis];
    
    NSMutableArray *highResDimensions = [NSMutableArray array];
    for (NSUInteger i=0; i<lowResTransformedVariable.dimensions.count; i++) {
        GLDimension *dim = [[GLDimension alloc] initAsDimension: newDimensions[i] transformedToBasis: [basis[i] unsignedIntegerValue] strictlyPositive: [lowResTransformedVariable.dimensions[i] isStrictlyPositive]];
        [highResDimensions addObject: dim];
    }
    
	GLVariable *resultVariable = [GLVariable variableOfType: lowResTransformedVariable.dataFormat withDimensions: highResDimensions forEquation: variable.equation];
	
	if (( self = [super initWithResult:resultVariable operand:lowResTransformedVariable] ))
	{		
		if ( highResDimensions.count == 2 )
		{
            NSUInteger nDataElements = self.result.nDataElements;
            NSUInteger nKPoints = [self.operand.dimensions[0] nPoints];
            NSUInteger nLPoints = [self.operand.dimensions[1] nPoints];
			NSUInteger n2KPoints = [self.result.dimensions[0] nPoints];
            NSUInteger n2LPoints = [self.result.dimensions[1] nPoints];
			
			NSUInteger diff = n2KPoints - nKPoints;
			
			if (resultVariable.isComplex)
			{
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLSplitComplex toSplit = splitComplexFromData(result);
					GLSplitComplex fromSplit = splitComplexFromData(operand);
					
					vGL_vclr( result.mutableBytes, 1, nDataElements);
					
					
					for (NSUInteger i=0; i<nKPoints; i++) {
						for (NSUInteger j=0; j<nLPoints; j++) {
							if (i<nKPoints/2) { // positive wavenumbers
								toSplit.realp[i*n2LPoints+j] = fromSplit.realp[i*nLPoints+j];
								toSplit.imagp[i*n2LPoints+j] = fromSplit.imagp[i*nLPoints+j];
							} else { // negative wavenumbers
								toSplit.realp[(i+diff)*n2LPoints+j] = fromSplit.realp[i*nLPoints+j];
								toSplit.imagp[(i+diff)*n2LPoints+j] = fromSplit.imagp[i*nLPoints+j];
							}
						}
					}
				};
			}
			else
			{
				self.blockOperation = ^(NSMutableData *result, NSData *operand) {
					GLFloat *to = result.mutableBytes;
					GLFloat *from = (GLFloat *) operand.bytes;
					
					vGL_vclr( result.mutableBytes, 1, nDataElements);

					for (NSUInteger i=0; i<nKPoints; i++) {
						for (NSUInteger j=0; j<nLPoints; j++) {
							if (i<nKPoints/2) { // positive wavenumbers
								to[i*n2LPoints+j] = from[i*nLPoints+j];
							} else { // negative wavenumbers
								to[(i+diff)*n2LPoints+j] = from[i*nLPoints+j];
							}
						}
					}
				};
			}
            self.graphvisDescription = [NSString stringWithFormat: @"zero pad"];
		}
		else {
			[NSException raise: @"NotYetImplemented" format: @"This case has not yet been implemented."];
		}
	}
	
    return self;
}

@end

