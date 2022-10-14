classdef NormalDistribution < Distribution
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma
    end
    
    methods
        function self = NormalDistribution(sigma)
            self.sigma = sigma;
            self.pdf = @(z) exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi));
            self.cdf = @(z) 0.5*(1 + erf(z/(sigma*sqrt(2))));
            self.w = @(z) sigma*sigma*ones(size(z));
            self.variance = sigma*sigma;
            
            self.dPDFoverZ = @(z) -exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi))/(sigma*sigma);
            self.logPDF = @(z) -(z.*z)/(2*sigma*sigma)-length(z)*(log(sigma)+log(2*pi)/2);
        end
        
        function y = rand(self,sz)
           y = self.sigma*randn(sz); 
        end

    end
end

