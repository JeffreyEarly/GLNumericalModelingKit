//
//  GLDifferentialOperator.m
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 4/27/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import "GLDifferentialOperator.h"

@implementation GLDifferentialOperator

- (GLBinaryOperation *) differentiationOperationFromVariable: (GLVariable *) operand
{
    return nil;
}

- (GLVariable *) differentiateVariable: (GLVariable *) operand
{
    GLBinaryOperation *diffOp = [self differentiationOperationFromVariable: operand];
    return diffOp.result;
}

@end