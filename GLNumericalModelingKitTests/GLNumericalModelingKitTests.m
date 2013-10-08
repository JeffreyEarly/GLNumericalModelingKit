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

#define fequal(a,b) (fabs((a) - (b)) < 10*FLT_EPSILON)

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

- (void)testVectorScalarAddition
{
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *xDim = [GLDimension dimensionXWithNPoints:4 length: 2.0];
	GLVariable *a = [GLVariable variableOfRealTypeWithDimensions: @[xDim] forEquation:equation];
	[a zero];
	
	GLVariable *b = [a scalarAdd: 1.0];
	[b solve];
	
	BOOL isTrue = 1;
	for (NSUInteger i=0; i<xDim.nPoints; i++) {
		isTrue &= (b.pointerValue[i] == 1.);
	}
	
    XCTAssertTrue( isTrue, @"No implementation for \"%s\"", __PRETTY_FUNCTION__);
}

- (void) test1DDiscreteFourierTransform
{
	GLEquation *equation = [[GLEquation alloc] init];
	
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 8 domainMin: 0.0 length: 1.0];
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLLinearTransform *matrix = [GLLinearTransform dftMatrixFromDimension: f.dimensions[0] forEquation: equation];
	GLVariable *f_tilde = [matrix transform: f];
	[f_tilde solve];
	GLSplitComplex output = f_tilde.splitComplex;
	
	GLVariable *k = [GLVariable variableOfRealTypeFromDimension:f_tilde.dimensions[0] withDimensions:f_tilde.dimensions forEquation:equation];
	[k solve];
	
	// We have only positive frequencies, and therefore this portion of the sine function should be negative.
	GLFloat expected_realp[8] = {1.0, 0., 0., 0., 0., 0., 0., 0.};
	GLFloat expected_imagp[8] = {0., 0., -1.5, 0., 0., 0., 1.5, 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + i%f, found %f + i%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
	
	GLLinearTransform *matrix_inverse = [GLLinearTransform idftMatrixFromDimension: f_tilde.dimensions[0] forEquation:equation];
	GLVariable *f_tilde_tilde = [matrix_inverse transform: f_tilde];
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
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[[x scalarMultiply: 2*2*M_PI] cos] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLLinearTransform *matrix = [GLLinearTransform cosineTransformMatrixFromDimension: f.dimensions[0] forEquation: equation];
	GLVariable *f_tilde = [matrix transform: f];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {1.0, 0., 0., 0., 1.5, 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLLinearTransform *matrix_inverse = [GLLinearTransform inverseCosineTransformMatrixFromDimension: f_tilde.dimensions[0] forEquation:equation];
	GLVariable *f_tilde_tilde = [matrix_inverse transform: f_tilde];
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
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLVariable *f_tilde = [f transformToBasis:@[@(kGLExponentialBasis)]];
	[f_tilde solve];
	GLSplitComplex output = f_tilde.splitComplex;
	
	GLVariable *k = [GLVariable variableOfRealTypeFromDimension:f_tilde.dimensions[0] withDimensions:f_tilde.dimensions forEquation:equation];
	[k solve];
	
	// We have only positive frequencies, and therefore this portion of the sine function should be negative.
	GLFloat expected_realp[5] = {1.0, 0., 0., 0., 0.};
	GLFloat expected_imagp[5] = {0., 0., -1.5, 0., 0.};
	
	for (int i=0; i<5; i++) {
		if ( !fequal(output.realp[i], expected_realp[i]) || !fequal(output.imagp[i], expected_imagp[i]) ) {
			XCTFail(@"Expected %f + i%f, found %f + i%f.", expected_realp[i], expected_imagp[i], output.realp[i], output.imagp[i]);
		}
	}
	
	GLVariable *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
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
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[[x scalarMultiply: 2*2*M_PI] cos] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLVariable *f_tilde = [f transformToBasis:@[@(kGLCosineBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {1.0, 0., 0., 0., 1.5, 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLVariable *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
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
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0];
	GLVariable *f_tilde = [f transformToBasis:@[@(kGLSineBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {0., 0., 0., 1.5, 0., 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLVariable *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
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
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[[x scalarMultiply: 2*2*M_PI] cos] scalarMultiply: 3.0] scalarAdd: 1.0];
	GLVariable *f_tilde = [f transformToBasis:@[@(kGLCosineHalfShiftBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {1.0, 0., 0., 0., 1.5, 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLVariable *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
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
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[x scalarMultiply: 2*2*M_PI] sin] scalarMultiply: 3.0];
	GLVariable *f_tilde = [f transformToBasis:@[@(kGLSineHalfShiftBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	// Note that we're expecting only *half* the power to appear (3.0/2) because we're assuming the other half is associated with the negative frequencies.
	GLFloat expected[8] = {0., 0., 0., 1.5, 0., 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLVariable *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
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
	
	GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: @[xDim] forEquation: equation];
	GLVariable *f = [[[x times: x] scalarMultiply: 2.0] scalarAdd: -1.0];
	GLVariable *f_tilde = [f transformToBasis:@[@(kGLChebyshevBasis)]];
	[f_tilde solve];
	GLFloat *output = f_tilde.pointerValue;
	
	GLFloat expected[8] = {0., 0., 0.5, 0., 0., 0., 0., 0.};
	
	for (int i=0; i<8; i++) {
		if ( !fequal(output[i], expected[i]) ) {
			XCTFail(@"Expected %f, found %f.", expected[i], output[i]);
		}
	}
	
	GLVariable *f_tilde_tilde = [f_tilde transformToBasis: @[@(kGLDeltaBasis)]];
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

@end
