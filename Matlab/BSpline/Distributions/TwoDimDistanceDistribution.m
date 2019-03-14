classdef TwoDimDistanceDistribution < Distribution
    %UNTITLED2 Summary of this class goes here
    %   Take a distribution defined on z=[-Inf -Inf], and converts it into
    %   radial distribution of two independent distribution on x-y.
    
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
            
            % scale our axis to be 20 times the standard deviation of the
            % 1 dimensional pdf.
%             resolution = 1e-10;
%             N_waypoints = 2;
%             N_points_between = 125;
%             N = N_waypoints*N_points_between+1;
%             
%             % Find out where the cdf is changing...
%             cdfMax = 1-resolution;
%             maxR = distribution1d.locationOfCDFPercentile(cdfMax);
%             minR = resolution*maxR;
%             cdfMin = distribution1d.cdf(minR);
%             waypointsR = [minR;zeros(N_waypoints-1,1);maxR];
%             for i=2:N_waypoints
%                 cdfValue = cdfMin + (i-1)*(cdfMax-cdfMin)/N_waypoints;
%                 waypointsR(i) = distribution1d.locationOfCDFPercentile(cdfValue);
%             end
%             
%             % ...and then place several points between those areas
%             r = zeros(N+1,1);
%             for i=1:N_waypoints
%                 indices = (((i-1)*N_points_between+1):(i*N_points_between))+1;
% %                 r_tmp = linspace(waypointsR(i),waypointsR(i+1),N_points_between+1)';
%                 r_tmp = 10.^(linspace(log10(waypointsR(i)),log10(waypointsR(i+1)),N_points_between+1))';
%                 r(indices) = r_tmp(1:end-1);
%             end
%             r(end) = maxR;           
            
resolution = 1e-10;
maxR = distribution1d.locationOfCDFPercentile(1-resolution);
N = 500;
            L = integral( @(z) z.*distribution1d.pdf(z),0,Inf );
            zmax = fminsearch( @(z) -z.*distribution1d.pdf(z),sqrt(distribution1d.variance) );
            f = @(z) TwoDimDistanceDistribution.vecIntegral( @(x) x.*distribution1d.pdf(x)/L,0,z );
            cdf_zmax = f(zmax);
            zverymax = TwoDimDistanceDistribution.fInverseBisection(f,1-resolution,0,maxR,1e-12);
            
            r1 = 10.^(linspace(log10(zmax*1e-3),log10(zmax),N/2))';
            dr = r1(end)-r1(end-1);
            r2 = [(r1(end):dr:zverymax)';zverymax];
            r = [0;r1;r2(2:end)];
            N = length(r)-1;
            
            % LEFT OFF HERE
            % The more the pdf is changing in R, the more points we need?
            
%             % to the pdf peak, we sample point logarithmically
%             r1 = 10.^(linspace(log10(zmax*1e-3),log10(zmax),floor(3*N/4)))';
%             dcdf = f(r1(end))-f(r1(end-1));
%             
%             % past the pdf peak, we space sample points evenly in cdf
%             % changes... BUT, we want dr to match at the transition. So
%             % what 
%             cdfsample2 = [(cdf_zmax:dcdf:(1-resolution))';(1-resolution)];
%             
%             r2 = TwoDimDistanceDistribution.fInverseBisection(f,cdfsample2,0,maxR,1e-12);
%            
%             r = [0;r1;r2(2:end)];
%             N = length(r)-1;
%             
%             cdfvalues = 10.^(linspace(log10(resolution),log10(1-resolution),N)');
%             cdfvalues = 10.^(linspace(log10(resolution),log10(1-resolution),N)');
%             r = [0;TwoDimDistanceDistribution.fInverseBisection(f,cdfvalues,0,maxR,1e-12)];
%             N = length(r)-1;
            
            %%%%%%%% LEFT OFF HERE
            % Need to compute the increments to grow with log, which will
            % set the smallest increment
            
%             r2 = TwoDimDistanceDistribution.fInverseBisection(f,linspace(0.25,1-resolution,N/2+1)',0,maxR,1e-12);
%             r1 = 10.^(linspace(log10(maxR*resolution),log10(r2(1)),N/2))';
%             r = [0;r1;r2(2:end)];
%             cdfsample1 = 10.^(linspace(log10(resolution),log10(cdf_zmax),N/2)'); 
%             
%             r = [0;TwoDimDistanceDistribution.fInverseBisection(f,[cdfsample1;cdfsample2(2:end)],0,maxR,1e-12)];
            
            %r = [0;TwoDimDistanceDistribution.fInverseBisection(distribution1d.cdf,linspace(0.5,1-resolution,N)',0,maxR,1e-12)];
            % create a logarithmic axis that includes zero
%             r = [0;10.^(linspace(log10(maxR*resolution),log10(maxR),N))'];
%             r = [0;10.^(linspace(log10(sqrt(maxR*resolution)),log10(sqrt(maxR)),N))'].^2;
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
%             
%             r(end+1) = Inf;
%             pdf(end+1) = 0;
%             cdf(end+1) = 1;
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

