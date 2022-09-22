classdef Lambda
    %UpperBoundary This class simply defines the possible valid
    %upper boundary conditions for the InternalModes classes.
    enumeration
        crossValidation         % minimize the expected mean-square error with cross-validation
        loglikelihood           % minimize using the log-likelihood
        optimalIterated         % minimize the expected mean-square error.
        
        optimalExpected         % take a guess at minimizing the mean-square error based on the effective sample-size.
        fullTensionIterated     % uses the anderson darling test on the interquartile region of the distribution
        fullTensionExpected     % take a guess at the full tension solution assuming infinite effective sample size.
        
        optimalRangedIterated       % minimize the expected mean-square error using points within the expected 99th percentile range
        optimalExpectedRobust       % Same as optimalExpected, but filters the signal to remove outliers.
        fullTensionExpectedRobust	% Same as fullTensionExpected, but filters the signal to remove outliers.
    end    
end

