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
+ (GLScalar *) scalarWithType: (GLDataFormat) format forEquation: (GLEquation *) anEquation;

- (GLScalar *) initWithType: (GLDataFormat) format forEquation: (GLEquation *) anEquation;
- (GLScalar *) initWithValue: (GLFloatComplex) aValue forEquation: (GLEquation *) anEquation;


- (id) dividedBy: (id) otherVariable;


/// C = -A
- (GLScalar *) negate;

/// C = abs(A)
- (GLScalar *) abs;

/// C = exp(A)
- (GLScalar *) exponentiate;

/// C = log(A)
- (GLScalar *) log;

/// C = sin(A)
- (GLScalar *) sin;

/// C = cos(A)
- (GLScalar *) cos;

/// C = atan(A)
- (GLScalar *) atan;

/// C = sqrt(A)
- (GLScalar *) sqrt;

@end
