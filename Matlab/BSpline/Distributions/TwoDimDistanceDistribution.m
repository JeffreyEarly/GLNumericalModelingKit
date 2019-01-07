classdef TwoDimDistanceDistribution < Distribution
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        distribution1d
    end
    
    methods
        function self = TwoDimDistanceDistribution(distribution1d)
            self.distribution1d = distribution1d;
            
            [r_, pdf_, cdf_] = TwoDimDistanceDistribution.TwoDimDistributionFromOneDimDistribution(distribution1d);
            
            self.pdf = @(z) interp1(r_,pdf_,z,'linear');
            self.cdf = @(z) interp1(r_,cdf_,z,'linear');
%             self.w = @(z) sigma*sigma*ones(size(z));
%             self.variance = (4-pi)*sigma*sigma/2;
            
%             self.dPDFoverZ = @(z) -exp(-(z.*z)/(2*sigma*sigma))/(sigma*sqrt(2*pi))/(sigma*sigma);
        end
    end
    
    methods (Static)
        function [r, pdf, cdf] = TwoDimDistributionFromOneDimDistribution(distribution1d)
            % returns the pdf and cdf of the two-dimensional t-distribution
            % as a function of distance, r.
            
            % This function needs to fast, and the cdf needs to be
            % accurate. So we use a strategy of putting points where
            % necessary (logarthimically placed), and using the cdf as the
            % fundamental quantity (rather than the pdf) because we can
            % adjust for discretization more easily (b/c it's already
            % integrated).
            
            % scale our axis to be 20 times the standard deviation of the
            % 1 dimensional pdf.
            maxR = distribution1d.locationOfCDFPercentile(0.9999);
            
            % This appears to be big enough, although we can go bigger if
            % needed.
            N = 250;
            
            % create a logarithmic axis that includes zero
            r = [0;10.^(linspace(log10(maxR/1e5),log10(maxR),N))'];
            x = [-flip(r(2:end),1); r];
            
            % evaluate the cdf on that axis
            pdf_int = distribution1d.cdf(x);
            
            % use that to define the pdf...
            dx = diff(x);
            pdf = diff(pdf_int)./dx;
            %...so that sum(pdf.*dx)=1.00000
            
            % now create a 2d version
            pdf2d_norm = (pdf .* reshape(pdf,1,[])) .* (dx .* reshape(dx,1,[]));
            % again, where sum(pdf2d_norm(:)) = 1.0000
            
            % create a 2d distance metric, rho
            y = [-flip(r(2:end)); r(2:end)];
            rho = sqrt( y.*y + reshape(y.*y,1,[]));
            
            % we are going to bin using the diagonal of rho, call it s.
            s = diag(rho);
            s = [0; s((N+1):end)];
            
            % Now sum up the energy in each bin
            pdf1d = zeros( size(s) );
            midS = [0; s(1:(end-1))+diff(s)/2];
            for i = 1:(length(s)-1)
                pdf1d(i) = sum( pdf2d_norm(  rho >= midS(i) & rho <= midS(i+1) ) );
            end
            % it must be true that sum(pdf1d) = 1.0000
            
            % but we want that sum(pdf1d .* diff(s)) = 1.000
            pdf1d(2:end) = pdf1d(2:end)./diff(s);
            cdf1d = [0; cumsum(pdf1d(2:end).*diff(s))];
            
            % now take care of the the fact that we used the diagnoal of a
            % 2d rectangle, rather than the inscribed circle.
            
            pdf = pdf1d(s<maxR);
            cdf = cdf1d(s<maxR);
            r = s(s<maxR);
        end
        
    end
end

