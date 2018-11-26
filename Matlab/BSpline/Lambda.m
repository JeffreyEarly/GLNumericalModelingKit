classdef Lambda
    %UpperBoundary This class simply defines the possible valid
    %upper boundary conditions for the InternalModes classes.
    enumeration
        optimal % minimize the expected mean-square error.
        initial % take a guess at minimizing the mean-square error based on the effective sample-size.
        full    % 'full tension'. Essentially assume infinite effective sample size.
    end    
end

