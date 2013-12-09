//
//  GLNumericalModelingKitTests.m
//  GLNumericalModelingKitTests
//
//  Created by Jeffrey J. Early on 9/30/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <XCTest/XCTest.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface GLNumericalModelingKitTests : XCTestCase

@end

// First check if their difference is less than precision, then do the same, but scaled by the magnitude.
#define fequal(a,b) ((fabs((a) - (b)) < 10*FLT_EPSILON) || (fabs(((a) - (b))/a) < 10*FLT_EPSILON))

// Set your own precision
#define fequalprec(a,b,prec) ((fabs((a) - (b)) < prec) || (fabs(((a) - (b))/a) < prec))

@implementation GLNumericalModelingKitTests

- (void)setUp
{
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown
{
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

/************************************************/
/*		Grid Tests								*/
/************************************************/

#pragma mark -
#pragma mark Grid Tests
#pragma mark

- (void) testEndpointGrid
{
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 4 domainMin: 0.0 length: 12.0];
	GLFloat expected[4] = {0., 4., 8., 12.};
	
	for (int i=0; i<4; i++) {
		if ( !fequal([dim valueAtIndex: i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], [dim valueAtIndex: i]);
		}
	}
}

- (void) testInteriorGrid
{
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLInteriorGrid nPoints: 4 domainMin: 0.0 length: 12.0];
	GLFloat expected[4] = {1.5, 4.5, 7.5, 10.5};
	
	for (int i=0; i<4; i++) {
		if ( !fequal([dim valueAtIndex: i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], [dim valueAtIndex: i]);
		}
	}
}

- (void) testPeriodicGrid
{
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 4 domainMin: 0.0 length: 12.0];
	GLFloat expected[4] = {0., 3., 6., 9.};
	
	for (int i=0; i<4; i++) {
		if ( !fequal([dim valueAtIndex: i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], [dim valueAtIndex: i]);
		}
	}
}

- (void) testChebyshevEndpointGrid
{
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: 4 domainMin: 0.0 length: 12.0];
	GLFloat expected[4] = {12., 9., 3., 0.};
	
	for (int i=0; i<4; i++) {
		if ( !fequal([dim valueAtIndex: i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], [dim valueAtIndex: i]);
		}
	}
}

- (void) testChebyshevInteriorGrid
{
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevInteriorGrid nPoints: 4 domainMin: 0.0 length: 12.0];
	GLFloat expected[4] = {11.543277195067720, 8.296100594190539, 3.703899405809461, 0.456722804932280};
	
	for (int i=0; i<4; i++) {
		if ( !fequal([dim valueAtIndex: i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], [dim valueAtIndex: i]);
		}
	}
}

/************************************************/
/*		Format Tests							*/
/************************************************/

#pragma mark -
#pragma mark Format Tests
#pragma mark

// Start with split format, check it, go to interleaved, check it, then go back to split format.
- (void) testSplitAndInterleavedFormatConversion
{
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *dim = [GLDimension dimensionXWithNPoints: 4 length: 4];
    GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
    GLFunction *result = [[var plus: @(1)] plus: [var swapComplex]];
    // result should contain (1 +I0, 2+I1, 3+I2, 4+I3)
    
    GLFloat expected_realp[4] = {1.0, 2.0, 3.0, 4.0};
	GLFloat expected_imagp[4] = {0.0, 1.0, 2.0, 3.0};
    
    GLSplitComplex output = result.splitComplex;
    for (int i=0; i<4; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
    
    result = [result interleavedFormat];
    
    GLFloatComplex *output2 = result.floatComplex;
    GLFloatComplex expected[4] = {1.0+I*0, 2.0+I*1.0, 3.0+I*2.0, 4.0+I*3.0};
    
    for (int i=0; i<4; i++) {
		if ( !fequal(creal(output2[i]), creal(expected[i])) || !fequal(cimag(output2[i]), cimag(expected[i])) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", creal(expected[i]), cimag(expected[i]), creal(output2[i]), cimag(output2[i]));
		}
	}

    result = [result splitFormat];
    
    output = result.splitComplex;
    for (int i=0; i<4; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
}

- (void) testDataTranspose
{
	GLDimension *xDim = [GLDimension dimensionXWithNPoints:5 length: 5.0];
	GLEquation *equation = [[GLEquation alloc] init];
	GLLinearTransform *A = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
	
	GLFloat a[5*5] =
	{	-1.01,   0.86,  -4.60,   3.31,  -4.81,
		3.98,   0.53,  -7.04,   5.29,   3.55,
		3.30,   8.26,  -3.89,   8.20,  -1.51,
		4.43,   4.96,  -7.66,  -7.33,   6.18,
		7.31,  -6.43,  -6.16,   2.47,   5.58};
	memcpy(A.pointerValue, a, A.nDataElements*sizeof(GLFloat));
	
	GLFloat b[5*5] = {
		-1.01f,  3.98f,  3.30f,  4.43f,  7.31f,
		0.86f,  0.53f,  8.26f,  4.96f, -6.43f,
		-4.60f, -7.04f, -3.89f, -7.66f, -6.16f,
		3.31f,  5.29f,  8.20f, -7.33f,  2.47f,
		-4.81f,  3.55f, -1.51f,  6.18f,  5.58f
	};
	
	GLLinearTransform *B = [A columnMajorOrdered];
	
	for (int i=0; i<5*5; i++) {
		if ( !fequalprec(B.pointerValue[i], b[i],1e-6) ) {
			XCTFail(@"Expected %f, found %f.", b[i], B.pointerValue[i]);
		}
	}
}

/************************************************/
/*		Algebriac Vector-Scalar Tests			*/
/************************************************/

#pragma mark -
#pragma mark Algebriac Vector-Scalar Tests
#pragma mark

- (void)testVectorScalarMultiplication
{
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:0.0 length:4.0];
	GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	
	GLFunction *result = [var times: @(2.0)];
	[result solve];
	GLFloat *output = result.pointerValue;
	
	GLFloat expected[4] = {0.0, 2.0, 4.0, 6.0};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

- (void)testVectorScalarAddition
{
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:0.0 length:4.0];
	GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	
	GLFunction *result = [var plus: @(1.0)];
	[result solve];
	GLFloat *output = result.pointerValue;
	
	GLFloat expected[4] = {1.0, 2.0, 3.0, 4.0};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

- (void)testVectorScalarDivisionCaseI
{
	// Case I: real vector (real variable type). Testing 1.0/(1, 2, 3, 4)
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:1.0 length:4.0];
	GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	
	GLFunction *result = [var scalarDivide: 1.0];
	[result solve];
	GLFloat *output = result.pointerValue;
	
	GLFloat expected[4] = {1.0/1.0, 1.0/2.0, 1.0/3.0, 1.0/4.0};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

- (void)testVectorScalarDivisionCaseII
{
	// Case II: real vector (complex variable type). Testing 1.0/(1, 2, 3, 4)
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:1.0 length:4.0];
	GLFunction *var = [GLFunction functionOfComplexTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	
	GLFunction *result = [var scalarDivide: 1.0];
	[result solve];
	GLSplitComplex output = result.splitComplex;
	
	GLFloat expected_realp[4] = {1.0/1.0, 1.0/2.0, 1.0/3.0, 1.0/4.0};
	GLFloat expected_imagp[4] = {0.0, 0.0, 0.0, 0.0};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
}

- (void)testVectorScalarDivisionCaseIII
{
	// Case III: imaginary vector (complex variable type). Testing 1.0/(I*(1, 2, 3, 4))
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:1.0 length:4.0];
	GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	
	GLFunction *result = [[var swapComplex] scalarDivide: 1.0];
	[result solve];
	GLSplitComplex output = result.splitComplex;
	
	GLFloat expected_realp[4] = {0.0, 0.0, 0.0, 0.0};
	GLFloat expected_imagp[4] = {-1.0/1.0, -1.0/2.0, -1.0/3.0, -1.0/4.0};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
}

- (void)testVectorScalarDivisionCaseIV
{
	// Case III: imaginary vector (complex variable type). Testing 1.0/(1+I*1, 2+I*2, 3+I*3, 4+I*4)
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:1.0 length:4.0];
	GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	
	GLFunction *result = [[var plus: [var swapComplex]] scalarDivide: 1.0];
	[result solve];
	GLSplitComplex output = result.splitComplex;
	
	GLFloat expected_realp[4] = {1.0/(2.0*1.0), 1.0/(2.0*2.0), 1.0/(2.0*3.0), 1.0/(2.0*4.0)};
	GLFloat expected_imagp[4] = {-1.0/(2.0*1.0), -1.0/(2.0*2.0), -1.0/(2.0*3.0), -1.0/(2.0*4.0)};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
}

- (void)testVectorPowerOperationCaseI
{
	// Case I: real vector (real variable type). Testing (1, 2, 3, 4)^2.0
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:1.0 length:4.0];
	GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	
	GLFunction *result = [var pow: 2.0];
	[result solve];
	GLFloat *output = result.pointerValue;
	
	GLFloat expected[4] = {pow(1.0, 2.0), pow(2.0, 2.0), pow(3.0, 2.0), pow(4.0, 2.0)};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

/************************************************/
/*		Interpolation Tests                     */
/************************************************/

#pragma mark -
#pragma mark Interpolation Tests
#pragma mark

- (void) test1DEndpointInterpolation
{
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints:5 domainMin:0.0 length:4.0];
    GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
    GLFunction *xpos = [x plus: @(-1.75)];
    GLFunction *interpolatedFunction = [x interpolateAtPoints: @[xpos]];
    
    // 0.0, 1.0, 2.0, 3.0, 4.0
    GLFloat expected[5] = {0.0, 0.0, 0.25, 1.25, 2.25};
    [interpolatedFunction solve];
    GLFloat *output = interpolatedFunction.pointerValue;
    for (int i=0; i<5; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
    
    xpos = [x plus: @(1.75)];
    interpolatedFunction = [x interpolateAtPoints: @[xpos]];
    
    // 0.0, 1.0, 2.0, 3.0, 4.0
    GLFloat expected2[5] = {1.75, 2.75, 3.75, 4.0, 4.0};
    [interpolatedFunction solve];
    output = interpolatedFunction.pointerValue;
    for (int i=0; i<5; i++) {
		if ( !fequal(output[i], expected2[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected2[i], output[i]);
		}
	}
}

- (void) test1DPeriodicInterpolation
{
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:5 domainMin:0.0 length:5.0];
    GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
    GLFunction *xpos = [x plus: @(-1.75)];
    GLFunction *interpolatedFunction = [x interpolateAtPoints: @[xpos]];
    
    // 0.0, 1.0, 2.0, 3.0, 4.0
    GLFloat expected[5] = {3.25, 3.00, 0.25, 1.25, 2.25};
    [interpolatedFunction solve];
    GLFloat *output = interpolatedFunction.pointerValue;
    for (int i=0; i<5; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
    
    xpos = [x plus: @(1.75)];
    interpolatedFunction = [x interpolateAtPoints: @[xpos]];
    
    // 0.0, 1.0, 2.0, 3.0, 4.0
    GLFloat expected2[5] = {1.75, 2.75, 3.75, 1.0, 0.75};
    [interpolatedFunction solve];
    output = interpolatedFunction.pointerValue;
    for (int i=0; i<5; i++) {
		if ( !fequal(output[i], expected2[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected2[i], output[i]);
		}
	}
}

- (void) test2DPeriodicInterpolation
{
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:5 domainMin:0.0 length:5.0];
	GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:5 domainMin:0.0 length:5.0];
    GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim,yDim] forEquation: equation];
	GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: @[xDim,yDim] forEquation: equation];
	
	// There are really 8 cases to check in a doubly periodic domain.
	GLDimension *interpDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints:8 domainMin:0.0 length:7.0];
	GLFunction *xInterp = [GLFunction functionOfRealTypeFromDimension: interpDim withDimensions: @[interpDim] forEquation: equation];
	GLFunction *yInterp = [GLFunction functionOfRealTypeFromDimension: interpDim withDimensions: @[interpDim] forEquation: equation];
	
	// How many extra times we want to wrap. Shouldn't make a difference.
	GLFloat n = 0.0;
	
	// Scenario 1
//	GLFloat DeltaMinus = -1.25 - n*5.0; // Should correspond to 3.75
//	GLFloat DeltaPlus = 5.25 + n*5.0;	// Should correspond to 0.25
//	GLFloat Middle = 2.5 + n*5.0;
//	
//	GLFloat DeltaMinusResult = 3.75;
//	GLFloat DeltaPlusResult = 0.25;
//	GLFloat MiddleResult = 2.5;
	
	// Scenario 2
	GLFloat DeltaMinus = -.25 - n*5.0; // Should correspond to 3.75
	GLFloat DeltaPlus = 4.25 + n*5.0;	// Should correspond to 0.25
	GLFloat Middle = 2.5 + n*5.0;
	
	GLFloat DeltaMinusResult = 1.0;
	GLFloat DeltaPlusResult = 3.0;
	GLFloat MiddleResult = 2.5;
	
	xInterp.pointerValue[0] = DeltaMinus; yInterp.pointerValue[0] = Middle;
	xInterp.pointerValue[1] = DeltaMinus; yInterp.pointerValue[1] = DeltaPlus;
	xInterp.pointerValue[2] = Middle; yInterp.pointerValue[2] = DeltaPlus;
	xInterp.pointerValue[3] = DeltaPlus; yInterp.pointerValue[3] = DeltaPlus;
	xInterp.pointerValue[4] = DeltaPlus; yInterp.pointerValue[4] = Middle;
	xInterp.pointerValue[5] = DeltaPlus; yInterp.pointerValue[5] = DeltaMinus;
	xInterp.pointerValue[6] = Middle; yInterp.pointerValue[6] = DeltaMinus;
	xInterp.pointerValue[7] = DeltaMinus; yInterp.pointerValue[7] = DeltaMinus;
	
	GLFloat expectedX[8] = {DeltaMinusResult, DeltaMinusResult, MiddleResult, DeltaPlusResult, DeltaPlusResult, DeltaPlusResult, MiddleResult, DeltaMinusResult};
	GLFloat expectedY[8] = {MiddleResult, DeltaPlusResult, DeltaPlusResult, DeltaPlusResult, MiddleResult, DeltaMinusResult, DeltaMinusResult, DeltaMinusResult};
	
    GLFunction *interpolatedX = [x interpolateAtPoints: @[xInterp, yInterp]];
	GLFunction *interpolatedY = [y interpolateAtPoints: @[xInterp, yInterp]];
    
    GLFloat *xout = interpolatedX.pointerValue;
	GLFloat *yout = interpolatedY.pointerValue;
    for (int i=0; i<8; i++) {
		if ( !fequal(xout[i], expectedX[i]) || !fequal(yout[i], expectedY[i]) ) {
			XCTFail(@"Case %d failed. Expected (%f, %f), found (%f, %f).", i, expectedX[i], expectedY[i], xout[i], yout[i]);
		}
	}
}

/************************************************/
/*		Minimization Tests						*/
/************************************************/

#pragma mark -
#pragma mark Minimization Tests
#pragma mark

- (void) test1DMinimization
{
	GLEquation *equation = [[GLEquation alloc] init];
	GLScalar *start = [GLScalar scalarWithValue: 1.0 forEquation: equation];
	// Starting the direction with 1.0 or 0.5 causes the algorithm to fail to converge!.
	GLScalar *direction = [GLScalar scalarWithValue: 0.3513 forEquation: equation];
	GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[start] withDeltas: @[direction] forFunction: ^(NSArray *xArray) {
		GLScalar *x = xArray[0];
		return [[x plus: @(-5)] times: [x plus: @(-5)]];
	}];
	
	GLScalar *minPoint = minimizer.result[0];
	GLScalar *minValue = minimizer.result[1];
		
	if ( !fequal(minPoint.pointerValue[0], 5.0) ) {
		XCTFail(@"Expected %f, found %f.", 5.0, minPoint.pointerValue[0]);
	}
	
	if ( !fequal(minValue.pointerValue[0], 0.0) ) {
		XCTFail(@"Expected %f, found %f.", 0.0, minValue.pointerValue[0]);
	}
}

- (void) test2DMinimization
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	NSMutableArray *start = [NSMutableArray array];
	start[0] = [GLScalar scalarWithValue: 3.0 forEquation: equation];
	start[1] = [GLScalar scalarWithValue: 10.0 forEquation: equation];
	
	NSMutableArray *direction = [NSMutableArray array];
	direction[0] = [GLScalar scalarWithValue: 0.351 forEquation: equation];
	direction[1] = [GLScalar scalarWithValue: 0.426 forEquation: equation];
	
	GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: start withDeltas: direction forFunction: ^(NSArray *xArray) {
		GLScalar *x = xArray[0];
		GLScalar *y = xArray[1];
		return [[[x plus: @(-1)] times: [x plus: @(-1)]] plus: [[y plus: @(-2)] times: [y plus: @(-2)]]];
	}];
	
	GLScalar *minXPoint = minimizer.result[0];
	GLScalar *minYPoint = minimizer.result[1];
	GLScalar *minValue = minimizer.result[2];
	
	if ( !fequal(minXPoint.pointerValue[0], 1.0) || !fequal(minYPoint.pointerValue[0], 2.0) ) {
		XCTFail(@"Expected (%f, %f), found (%f, %f).", 1.0, 2.0, minXPoint.pointerValue[0], minYPoint.pointerValue[0]);
	}
	
	if ( !fequal(minValue.pointerValue[0], 0.0) ) {
		XCTFail(@"Expected %f, found %f.", 0.0, minValue.pointerValue[0]);
	}
}

/************************************************/
/*		Transform Tests							*/
/************************************************/

#pragma mark -
#pragma mark Transform Tests
#pragma mark

- (void) test1DDiscreteFourierTransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 8 domainMin: 0.0 length: 1.0];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[[x times: @(2*2*M_PI)] sin] times: @(3.0)] plus: @(1.0)];
    GLLinearTransform *matrix = [GLLinearTransform discreteTransformFromDimension: f.dimensions[0] toBasis: kGLExponentialBasis forEquation:equation];
	GLFunction *f_tilde = [matrix transform: f];
	[f_tilde solve];
	GLSplitComplex output = f_tilde.splitComplex;
	
	GLFunction *k = [GLFunction functionOfRealTypeFromDimension:f_tilde.dimensions[0] withDimensions:f_tilde.dimensions forEquation:equation];
	[k solve];
	
	// We have only positive frequencies, and therefore this portion of the sine function should be negative.
	GLFloat expected_realp[8] = {1.0, 0., 0., 0., 0., 0., 0., 0.};
	GLFloat expected_imagp[8] = {0., 0., -1.5, 0., 0., 0., 1.5, 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + i%f, found %f + i%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
	
    GLLinearTransform *matrix_inverse = [GLLinearTransform discreteTransformFromDimension: f_tilde.dimensions[0] toBasis: kGLDeltaBasis forEquation:equation];
	GLFunction *f_tilde_tilde = [matrix_inverse transform: f_tilde];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	GLFloat *output2 = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output2[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output2[i]);
		}
	}
}

- (void) test1DDiscreteCosineIIITransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLInteriorGrid nPoints: 8 domainMin: 0.0 length: 1.0];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[[x scalarMultiply: 2*2*M_PI] cos] scalarMultiply: 3.0] scalarAdd: 1.0];
    GLLinearTransform *matrix = [GLLinearTransform discreteTransformFromDimension: f.dimensions[0] toBasis: kGLCosineHalfShiftBasis forEquation:equation];
	GLFunction *f_tilde = [matrix transform: f];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {1.0, 0., 0., 0., 1.5, 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
    GLLinearTransform *matrix_inverse = [GLLinearTransform discreteTransformFromDimension: f_tilde.dimensions[0] toBasis: kGLDeltaBasis forEquation:equation];
	GLFunction *f_tilde_tilde = [matrix_inverse transform: f_tilde];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	output = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output[i]);
		}
	}
}

//- (void) test1DDiscreteSineIIITransform
//{
//	GLEquation *equation = [[GLEquation alloc] init];
//	
//	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLInteriorGrid nPoints: 8 domainMin: 0.0 length: 1.0];
//	
//	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
//	GLVariable *f = [[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0];
//	GLLinearTransform *matrix = [GLLinearTransform sineTransformMatrixFromDimension: f.dimensions[0] forEquation: equation];
//	GLVariable *f_tilde = [matrix transform: f];
//	[f_tilde solve];
//	GLFloat *output = f_tilde.pointerValue;
//	
//	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
//	GLFloat expected[8] = {0., 0., 0., 1.5, 0., 0., 0., 0.};
//	
//	for (int i=0; i<8; i++) {
//		if ( !fequal(output[i], expected[i]) ) {
//			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
//		}
//	}
//	
//	GLLinearTransform *matrix_inverse = [GLLinearTransform inverseSineTransformMatrixFromDimension: f_tilde.dimensions[0] forEquation:equation];
//	GLVariable *f_tilde_tilde = [matrix_inverse transform: f_tilde];
//	[f_tilde_tilde solve];
//	
//	// Check that the inverse transformation works as expected.
//	GLFloat *input = f.pointerValue;
//	output = f_tilde_tilde.pointerValue;
//	for (int i=0; i<8; i++) {
//		if ( !fequal(input[i], output[i]) ) {
//			XCTFail(@"Expected %f, found %f.", input[i], output[i]);
//		}
//	}
//}


- (void) test1DFastFourierTransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 8 domainMin: 0.0 length: 1.0];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLFunction *f_tilde = [f transformToBasis:@[@(kGLExponentialBasis)]];
	[f_tilde solve];
	GLSplitComplex output = f_tilde.splitComplex;
	
	GLFunction *k = [GLFunction functionOfRealTypeFromDimension:f_tilde.dimensions[0] withDimensions:f_tilde.dimensions forEquation:equation];
	[k solve];
	
	// We have only positive frequencies, and therefore this portion of the sine function should be negative.
	GLFloat expected_realp[5] = {1.0, 0., 0., 0., 0.};
	GLFloat expected_imagp[5] = {0., 0., -1.5, 0., 0.};
	
	for (int i=0; i<5; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + i%f, found %f + i%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
	
	GLFunction *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	GLFloat *output2 = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output2[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output2[i]);
		}
	}
}

- (void) test1DFastCosineITransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 8 domainMin: 0.0 length: 1.0];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[[x scalarMultiply: 2*2*M_PI] cos] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLFunction *f_tilde = [f transformToBasis:@[@(kGLCosineBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {1.0, 0., 0., 0., 1.5, 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLFunction *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	output = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output[i]);
		}
	}
}

- (void) test1DFastSineITransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	// The DST-I transform requires a non-standard interior grid.
	GLFloat sampleInterval = 1./(8.+1);
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 8 domainMin: sampleInterval length: 1.0-2*sampleInterval];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0];
	GLFunction *f_tilde = [f transformToBasis:@[@(kGLSineBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {0., 0., 0., 1.5, 0., 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLFunction *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	output = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output[i]);
		}
	}
}

- (void) test1DFastCosineIIITransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLInteriorGrid nPoints: 8 domainMin: 0.0 length: 1.0];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[[x scalarMultiply: 2*2*M_PI] cos] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLFunction *f_tilde = [f transformToBasis:@[@(kGLCosineHalfShiftBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {1.0, 0., 0., 0., 1.5, 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLFunction *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	output = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output[i]);
		}
	}
}

- (void) test1DFastSineIIITransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLInteriorGrid nPoints: 8 domainMin: 0.0 length: 1.0];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0];
	GLFunction *f_tilde = [f transformToBasis:@[@(kGLSineHalfShiftBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {0., 0., 0., 1.5, 0., 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLFunction *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	output = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output[i]);
		}
	}
}

- (void) test1DFastChebyshevTransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: 8 domainMin: -1.0 length: 2.0];
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[x times: x] times: @(2.0)] plus: @(-1.0)];
	GLFunction *f_tilde = [f transformToBasis:@[@(kGLChebyshevBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	GLFloat expected[8] = {0., 0., 0.5, 0., 0., 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLFunction *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
	[f_tilde_tilde solve];
	
	// Check that the inverse transformation works as expected.
	GLFloat *input = f.pointerValue;
	output = f_tilde_tilde.pointerValue;
	for (int i=0; i<8; i++) {
		if ( !fequal(input[i], output[i]) ) {
			XCTFail(@"Expected %f, found %f.", input[i], output[i]);
		}
	}
}

/************************************************/
/*		Linear Algebra							*/
/************************************************/

#pragma mark -
#pragma mark Linear Algebra
#pragma mark

-(void) testRealMatrixMultiplication
{
	GLDimension *xDim = [GLDimension dimensionXWithNPoints:2 length: 2.0];
	GLDimension *yDim = [GLDimension dimensionYWithNPoints:3 length: 3.0];
	
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLLinearTransform *A = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[yDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
	GLLinearTransform *B = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[yDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
	
	GLFloat *a = A.pointerValue;
	a[0] = 0.11; a[1] = 0.12; a[2] = 0.13;
	a[3] = 0.21; a[4] = 0.22; a[5] = 0.23;
	
	GLFloat *b = B.pointerValue;
	b[0] = 1011; b[1] = 1012;
	b[2] = 1021; b[3] = 1022;
	b[4] = 1031; b[5] = 1032;
	
	GLLinearTransform *C = [A times: B];
	GLFloat *output = C.pointerValue;
	[C solve];
	
	GLFloat expected[4] = { 367.76, 368.12, 674.06, 674.72 };
	   
	for (int i=0; i<4; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

-(void) testComplexMatrixMultiplication
{
	GLDimension *xDim = [GLDimension dimensionXWithNPoints:2 length: 2.0];
	GLDimension *yDim = [GLDimension dimensionYWithNPoints:3 length: 3.0];
	
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLLinearTransform *A = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: @[yDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
	GLLinearTransform *B = [GLLinearTransform transformOfType: kGLSplitComplexDataFormat withFromDimensions: @[xDim] toDimensions: @[yDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
	
	GLSplitComplex x = A.splitComplex;
	x.realp[0] = 1.0; x.realp[1] = 0.0; x.realp[2] = 1.0;
	x.realp[3] = 2.0; x.realp[4] = 0.0; x.realp[5] = 2.0;
	
	x.imagp[0] = 0.0; x.imagp[1] = 1.0; x.imagp[2] = 0.0;
	x.imagp[3] = 0.0; x.imagp[4] = 1.0; x.imagp[5] = 0.0;
	
	GLSplitComplex y = B.splitComplex;
	y.realp[0] = 0.0; y.realp[1] = 0.0;
	y.realp[2] = 1.0; y.realp[3] = 2.0;
	y.realp[4] = 1.0; y.realp[5] = 2.0;
	
	y.imagp[0] = 2.0; y.imagp[1] = 3.0;
	y.imagp[2] = 0.0; y.imagp[3] = 0.0;
	y.imagp[4] = 0.0; y.imagp[5] = 0.0;
	
	GLLinearTransform *C = [A times: B];
	GLSplitComplex output = C.splitComplex;
	[C solve];
	
	GLFloat expected_realp[4] = {1.0, 2.0, 2.0, 4.0};
	GLFloat expected_imagp[4] = {3.0, 5.0, 5.0, 8.0};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
}

- (void) testRealMatrixInversion
{
	GLDimension *xDim = [GLDimension dimensionXWithNPoints:3 length: 3.0];
	GLDimension *yDim = [GLDimension dimensionYWithNPoints:3 length: 3.0];
	
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLLinearTransform *A = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[yDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
	
	GLFloat *a = A.pointerValue;
	a[0] = 1.0; a[1] = 2.0; a[2] = 0.0;
	a[3] = 2.0; a[4] = 5.0; a[5] = -1.0;
	a[6] = 4.0; a[7] = 10.0; a[8] = -1.0;
	
	GLLinearTransform *invA = [A inverse];
	
	[invA solve];
	[equation waitUntilAllOperationsAreFinished];
	
	GLFloat *output = invA.pointerValue;
	
	GLFloat expected[9] = { 5., 2., -2., -2., -1., 1., 0., -2., 1. };
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

- (void) testRealMatrixEigendecomposition
{
	GLDimension *xDim = [GLDimension dimensionXWithNPoints:5 length: 5.0];
	GLEquation *equation = [[GLEquation alloc] init];
	GLLinearTransform *A = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[xDim] toDimensions: @[xDim] inFormat: @[@(kGLDenseMatrixFormat)] forEquation: equation matrix: nil];
	
	GLFloat a[5*5] =
		{	-1.01,   0.86,  -4.60,   3.31,  -4.81,
			3.98,   0.53,  -7.04,   5.29,   3.55,
			3.30,   8.26,  -3.89,   8.20,  -1.51,
			4.43,   4.96,  -7.66,  -7.33,   6.18,
			7.31,  -6.43,  -6.16,   2.47,   5.58};
	memcpy(A.pointerValue, a, A.nDataElements*sizeof(GLFloat));
	
	NSArray *system = [A eigensystem];
	GLFunction *eigenvalues = system[0];
	GLLinearTransform *S = system[1];
	
	GLFloat expected_eigval_realp[5] = { 2.858132878034353,
		2.858132878034353,
		-0.686674513305952,
		-0.686674513305952,
		-10.462916729456815};
	
	GLFloat expected_eigval_imagp[5] = {  10.762749830715675,
		-10.762749830715675,
		4.704261340628109,
		-4.704261340628109,
		0};
	
	GLFloat expected_eigvec_realp[5*5] =
    {   0.108064791301352,   0.108064791301352,   0.732233989783721,   0.732233989783721,   0.460646436627129,
		0.406312881322675,   0.406312881322675,  -0.026463011089023,  -0.026463011089023,   0.337703828285972,
		0.102357685061564,   0.102357685061564,   0.191648728080536,   0.191648728080536,   0.308743941854130,
		0.398631098084136,   0.398631098084136,  -0.079011062984309,  -0.079011062984309,  -0.743845837531074,
		0.539535056047412,   0.539535056047412,  -0.291604754325538,  -0.291604754325538,   0.158529281647889};
	
	GLFloat expected_eigvec_imagp[5*5] =
    {    0.168648343501007,  -0.168648343501007,   0,                   0,                 0,
		-0.259009768920532,   0.259009768920532,  -0.016946754378572,   0.016946754378572, 0,
		-0.508802314178709,   0.508802314178709,  -0.292565995475612,   0.292565995475612, 0,
		-0.091333452369541,   0.091333452369541,  -0.078075936426824,   0.078075936426824, 0,
		 0,                   0,                  -0.493102293052802,   0.493102293052802, 0};
	
	GLSplitComplex eigval = eigenvalues.splitComplex;
	for (int i=0; i<5; i++) {
		if ( !fequalprec(eigval.realp[i], expected_eigval_realp[i],1e-5) || !fequalprec(eigval.imagp[i], expected_eigval_imagp[i],1e-5) ) {
			XCTFail(@"Expected %f + I%f, found %f + I%f.", expected_eigval_realp[i], expected_eigval_imagp[i], eigval.realp[i], eigval.imagp[i]);
		}
	}
	
    GLFloat reltol = 1e-6;
	GLSplitComplex eigvec = S.splitComplex;
	for (NSUInteger i=0; i<5; i++) {
        BOOL didFlipSign = NO;
        for (NSUInteger j=0; j<5; j++) {
            BOOL isValid = fequalprec(eigvec.realp[j*5+i], expected_eigvec_realp[j*5+i],reltol) && fequalprec(eigvec.imagp[j*5+i], expected_eigvec_imagp[j*5+i],reltol);
            BOOL isValidFlipped = fequalprec(eigvec.realp[j*5+i], -expected_eigvec_realp[j*5+i],reltol) && fequalprec(eigvec.imagp[j*5+i], -expected_eigvec_imagp[j*5+i],reltol);
            if (j==0) {
                if (!isValid && isValidFlipped) {
                    didFlipSign = YES;
                }
            }
            
            if ( !((isValid && !didFlipSign) || (isValidFlipped && didFlipSign)) ) {
                XCTFail(@"(%lu,%lu), Expected %f + I%f, found %f + I%f.",i/5, i%5, expected_eigvec_realp[j*5+i], expected_eigvec_imagp[j*5+i], eigvec.realp[j*5+i], eigvec.imagp[j*5+i]);
            }
        }
	}
}

/************************************************/
/*		Differential							*/
/************************************************/

#pragma mark -
#pragma mark Differential
#pragma mark

- (void) test1DExponentialBasisDifferentiation
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 8 domainMin: 0.0 length: 1.0];
    xDim.name = @"x";
	
	GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLFunction *f = [[[[x times: @(2*2*M_PI)] sin] times: @(3.0)] plus: @(1.0)];
    
    GLFunction *fdiff = [[f x] spatialDomain];
    GLFunction *fdiff_expected = [[[x times: @(2*2*M_PI)] cos] times: @(3.0*2*2*M_PI)];
    
    [fdiff solve];
    [fdiff_expected solve];
    
    GLFloat *output = fdiff.pointerValue;
	GLFloat *expected = fdiff_expected.pointerValue;
    
	for (int i=0; i<4; i++) {
		if ( !fequalprec(output[i], expected[i], 1e-5) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

- (void) testFiniteDifferencing
{
	GLEquation *equation = [[GLEquation alloc] init];
	NSUInteger N=8;
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: N domainMin: 0.0 length: N-1];
    xDim.name = @"x";
	
	GLLinearTransform *diffX = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:xDim forEquation:equation];
	[diffX dumpToConsole];
	
	[[diffX densified] dumpToConsole];
	
//	GLFloat *a = diffX.pointerValue;
//	for (NSUInteger i=0; i<N; i++) {
//		for (NSUInteger iDiagonal=0; iDiagonal<3; iDiagonal++) {
//			printf("%6.2f\t", a[iDiagonal*N+i]);
//		}
//		printf("\n");
//	}
}

/************************************************/
/*		Optimizer                               */
/************************************************/

#pragma mark -
#pragma mark Optimizer
#pragma mark

- (void) testOperationOptimizer
{
    GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *dim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:4 domainMin:0.0 length:4.0];
	GLFunction *var = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: @[dim] forEquation:equation];
	GLFunction *result = [[var times: @(2.0)] plus: @(1.0)];
    
    Class newOperationClass = [GLVariableOperation variableOperationSubclassWithOperand: @[var] result: @[result]];
    GLVariableOperation *newOperation = [[newOperationClass alloc] initWithOperand: @[var]];
    
    GLFunction *newResult = newOperation.result[0];
	[newResult solve];
	GLFloat *output = newResult.pointerValue;
	
	GLFloat expected[4] = {1.0, 3.0, 5.0, 7.0};
	
	for (int i=0; i<4; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
    
//    GLOperationVisualizer *vizualizer = [[GLOperationVisualizer alloc] initWithTopVariables: @[var] bottomVariables:@[result]];
//    NSLog(@"%@", vizualizer.graphvisDescription);
}

- (void) testRungeKutta4thOrderIntegration
{
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 64 domainMin: 0 length: 20];
    xDim.name = @"x";
    GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation:equation];
    GLFloat x0 = 11;
    GLFunction *gaussian = [[[[x plus: @(-x0)] times: [x plus: @(-x0)]] negate] exponentiate];
    
    GLFloat cfl = 0.25;
    GLFloat timeStep = cfl * xDim.sampleInterval;
    
    GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: @[[gaussian frequencyDomain]] stepSize: timeStep fFromTY: ^(GLScalar *t, NSArray *ynew) {
        return @[[ynew[0] diff:@"x"]];
    }];
    
    // We stepped forward 20 steps in time
    gaussian = [[integrator stepForwardToTime: timeStep*20][0] spatialDomain];
    
    // Compute the analytical solution 20 steps forward
    GLFloat x10 = x0-integrator.currentTime;
    GLFunction *gaussian10 = [[[[x plus: @(-x10)] times: [x plus: @(-x10)]] negate] exponentiate];
    
    GLFloat *output = gaussian.pointerValue;
    GLFloat *expected = gaussian10.pointerValue;
    
    // We expected 4th order Runge-Kutta to give use relative accuracies of 10^(-4)
    for (int i=0; i<xDim.nPoints; i++) {
		if ( !fequalprec(output[i], expected[i], 1e-4) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

- (void) testRungeKutta23Integration
{
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 64 domainMin: 0 length: 20];
    xDim.name = @"x";
    GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation:equation];
    GLFloat x0 = 11;
    GLFunction *gaussian = [[[[x plus: @(-x0)] times: [x plus: @(-x0)]] negate] exponentiate];
    
    GLFloat cfl = 0.25;
    GLFloat timeStep = cfl * xDim.sampleInterval;
    
    GLRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[[gaussian frequencyDomain]] stepSize: timeStep fFromTY: ^(GLScalar *t, NSArray *ynew) {
        return @[[ynew[0] diff:@"x"]];
    }];
    
    // We stepped forward 20 steps in time
    gaussian = [[integrator stepForwardToTime: timeStep*20][0] spatialDomain];
    
    // Compute the analytical solution 20 steps forward
	// Do NOT extract the 'currentTime' from the integrator. It will have returned an interpolated value!
    GLFloat x10 = x0-timeStep*20;
    GLFunction *gaussian10 = [[[[x plus: @(-x10)] times: [x plus: @(-x10)]] negate] exponentiate];
    
    GLFloat *output = gaussian.pointerValue;
    GLFloat *expected = gaussian10.pointerValue;
    
    // We expected 4th order Runge-Kutta to give use relative accuracies of 10^(-4)
    for (int i=0; i<xDim.nPoints; i++) {
		if ( !fequalprec(output[i], expected[i], 1e-3) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
}

@end




