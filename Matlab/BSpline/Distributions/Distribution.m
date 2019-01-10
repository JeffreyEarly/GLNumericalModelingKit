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
        logPDF      % log of the pdf
    end
    
    methods
        function z = locationOfCDFPercentile(self, alpha)
            assert( alpha > 0 & alpha < 1,'alpha must be between 0 and 1');       
            z = fminsearch(@(z) abs(self.cdf(z)-alpha),0);
        end
        
        function var = varianceInRange(self,zmin,zmax)
           var = integral( @(z) z.*z.*self.pdf(z),zmin,zmax);
        end
        
        function var = varianceInPercentileRange(self,pctmin,pctmax)
            zmin = self.locationOfCDFPercentile(pctmin);
            zmax = self.locationOfCDFPercentile(pctmax);
            var = self.varianceInRange(zmin,zmax);
        end
        
        function totalError = kolmogorovSmirnovError(self,epsilon,zmin,zmax)
            if nargin == 4
                x = sort(epsilon( epsilon > zmin & epsilon < zmax ));
            else
                x = sort(epsilon);
            end
            
            n = length(x);
            y_data = (1:n)'/n;
            y = self.cdf(x);
            
            D = max(abs(y-y_data));
            
            totalError = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
        end
    end
end

