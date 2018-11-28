classdef Lambda
    %UpperBoundary This class simply defines the possible valid
    %upper boundary conditions for the InternalModes classes.
    enumeration
        optimalIterated % minimize the expected mean-square error.
        optimalExpected % take a guess at minimizing the mean-square error based on the effective sample-size.
        fullTension     % 'full tension'. Essentially assume infinite effective sample size.
    end    
end

