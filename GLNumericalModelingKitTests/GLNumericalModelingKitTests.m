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

@end
