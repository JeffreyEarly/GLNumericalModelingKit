//
//  GLScalar.h
//  GLNumericalModelingKit
//
//  Created by Jeffrey J. Early on 10/11/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//


#import <GLNumericalModelingKit/GLVariable.h>

@interface GLScalar : GLVariable

+ (GLScalar *) scalarWithValue: (GLFloatComplex) aValue forEquation: (GLEquation *) anEquation;

- (GLScalar *) initWithType: (GLDataFormat) format forEquation: (GLEquation *) anEquation;
- (GLScalar *) initWithValue: (GLFloatComplex) aValue forEquation: (GLEquation *) anEquation;


- (id) dividedBy: (id) otherVariable;

@end
