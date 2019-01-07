classdef AddedDistribution < Distribution
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        scalings
        distributions
    end
    
    methods
        function self = AddedDistribution(scalings,distribution1,distribution2,varargin)
            if length(scalings) == 1
                scalings(2) = 1-scalings;
            end
            if any(scalings<0) || sum(scalings) ~= 1.0
                error('The scalings must sum to 1.0, and all values must be between 0 and 1.');
            end
            self.scalings = scalings;
            
            nDistributions = 2 + length(varargin);
            if nDistributions ~= length(scalings)
                error('There must be a scaling value for all distrubtions.');
            end
            
            self.distributions = cell(nDistributions,1);
            self.distributions{1} = distribution1;
            self.distributions{2} = distribution2;
            for i = 1:length(varargin)
                self.distributions{i+2} = varargin{i};
            end
            
            self.pdf = @(z) self.SummedPDF(z);
            self.cdf = @(z) self.SummedCDF(z);
            self.dPDFoverZ = @(z) self.SummeddPDFoverZ(z);
            self.w = @(z) self.SummedWeight(z);
            self.logPDF = @(z) log(self.pdf(z));
            
            self.variance = 0;
            for i=1:length(self.distributions)
                self.variance = self.variance + self.scalings(i)*self.distributions{i}.variance;
            end
            
        end
        
        function pdf = SummedPDF(self,z)
            pdf = zeros(size(z));
            for i=1:length(self.distributions)
               pdf = pdf + self.scalings(i)*self.distributions{i}.pdf(z); 
            end
        end
                
        function cdf = SummedCDF(self,z)
            cdf = zeros(size(z));
            for i=1:length(self.distributions)
                cdf = cdf + self.scalings(i)*self.distributions{i}.cdf(z);
            end
        end
        
        function dPDFoverZ = SummeddPDFoverZ(self,z)
            dPDFoverZ = zeros(size(z));
            for i=1:length(self.distributions)
                dPDFoverZ = dPDFoverZ + self.scalings(i)*self.distributions{i}.dPDFoverZ(z);
            end
        end
        
        function w = SummedWeight(self,z)
            w = -self.pdf(z) ./ self.dPDFoverZ(z);
        end
        
    end
end