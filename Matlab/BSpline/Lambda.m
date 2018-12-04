classdef Lambda
    %UpperBoundary This class simply defines the possible valid
    %upper boundary conditions for the InternalModes classes.
    enumeration
        optimalIterated         % minimize the expected mean-square error.
        optimalExpected         % take a guess at minimizing the mean-square error based on the effective sample-size.
        fullTensionExpected     % take a guess at the full tension solution assuming infinite effective sample size.
    end    
end

