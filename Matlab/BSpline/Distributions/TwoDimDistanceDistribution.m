classdef TwoDimDistanceDistribution < Distribution
    %UNTITLED2 Summary of this class goes here
    %   Take a distribution defined on z=[-Inf -Inf], and converts it into
    %   radial distribution of two independent distribution on x-y.
    %   The resulting cdf is far more accurate than the resulting pdf.
    
    properties
        distribution1d
    end
    
    methods
        function self = TwoDimDistanceDistribution(distribution1d)
            self.distribution1d = distribution1d;
            
            [r_, pdf_, cdf_] = TwoDimDistanceDistribution.TwoDimDistributionFromOneDimDistribution(distribution1d);
            
            self.pdf = @(z) interp1(r_,pdf_,z,'linear',0);
            self.cdf = @(z) interp1(r_,cdf_,z,'linear',0);
            self.zrange = [0 Inf];
            self.variance = self.varianceInRange(0,Inf);
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
            
            resolution = 1e-5;
            maxR = distribution1d.locationOfCDFPercentile(1-1e-10);
            N = 500;

            r = [0;10.^(linspace(log10(maxR*resolution),log10(maxR),N))'];

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
            pdf1d(2:end) = RunningFilter(pdf1d(2:end)./diff(s),3,'mean');
            cdf1d = [0; cumsum(pdf1d(2:end).*diff(s))];
            
            % now take care of the the fact that we used the diagnoal of a
            % 2d rectangle, rather than the inscribed circle.
            
            pdf = pdf1d(s<maxR);
            cdf = cdf1d(s<maxR);
            r = s(s<maxR);
            
            r(end+1) = 2*r(end);
            pdf(end+1) = 0;
            cdf(end+1) = 1;
        end
       
        function y = vecIntegral(f,alim,blim)
           y = zeros(size(blim));
           for i=1:length(y)
              y(i) = integral(f,alim,blim(i)); 
           end
        end
        
        function y = fInverseBisection(f, x, yMin,yMax, tol)
            %FINVERSEBISECTION(F, X)   Compute F^{-1}(X) using Bisection.
            % Taken from cumsum as part of chebfun.
            % chebfun/inv.m
            %
            % Copyright 2017 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.
            
            a = yMin*ones(size(x));
            b = yMax*ones(size(x));
            c = (a + b)/2;
            
            while ( norm(b - a, inf) >= tol )
                vals = feval(f, c);
                % Bisection:
                I1 = ((vals-x) <= -tol);
                I2 = ((vals-x) >= tol);
                I3 = ~I1 & ~I2;
                a = I1.*c + I2.*a + I3.*c;
                b = I1.*b + I2.*c + I3.*c;
                c = (a+b)/2;
            end
            
            y = c;
            
        end
        
    end
end

