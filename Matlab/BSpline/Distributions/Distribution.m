classdef (Abstract) Distribution < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pdf     % probability density function
        cdf     % cumulative distribution function
        rho     % autocorrelation function
        
        variance    % total variance

        w           % 'weight' function, 
        dPDFoverZ   % derivative of the pdf wrt z, divided by z
    end
    
    methods
        function z = LocationOfCDFPercentile(self, alpha)
            assert( alpha > 0 & alpha < 1,'alpha must be between 0 and 1');       
            z = fminsearch(@(z) abs(self.cdf(z)-alpha),0);
        end
        
        function var = VarianceInRange(self,zmin,zmax)
           var = integral( @(z) z.*z.*self.pdf(z),zmin,zmax);
        end
        
        function var = VarianceInPercentileRange(self,pctmin,pctmax)
            zmin = self.LocationOfCDFPercentile(pctmin);
            zmax = self.LocationOfCDFPercentile(pctmax);
            var = self.VarianceInRange(zmin,zmax);
        end
    end
end

