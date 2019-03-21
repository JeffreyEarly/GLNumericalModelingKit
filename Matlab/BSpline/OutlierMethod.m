classdef OutlierMethod
    %OutlierMethod This class defines a few of the built-in outlier
    %detection strategies.
    enumeration
        doNothing                            % Don't do anything
        sigmaFullTensionMethod               % Prime the IRLS with a new sigma taken from the full tension solution.
        sigmaRejectionMethod                 % Prime the IRLS with a new sigma, tuned to outliers. Takes rejectionPDFRatio as argument.
        distributionMethod                   % Change the distribution to the added distribution
        distributionWeightAndSigmaMethod     % Change the distribution, sigma, and the weight function, tuned to outliers. Takes rejectionPDFRatio as argument.
        knotRemovalMethod                    % Remove knot points from suspected outliers. Takes rejectionPDFRatio as argument.
        weightingMethod                      % Weight the minimizers according the likelihood the point is an outlier
    end    
end

