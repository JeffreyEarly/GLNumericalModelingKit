classdef BSpline
    %BSPLINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = untitled2(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    
    methods (Static)

        function x = EvaluateSplineAtPoints(t, m, t_knot, K, j)
           %% Evaluate the spline at given points
           %
           % Returns the value x at points t for a spline of order K with
           % coefficient m. The knot points for the spline are t_knot. This
           % will return the j-th derivative (which can be zero).
           if j>K
               error('The derivative requested is higher than the spline order');
           end
           
           if issorted(t,'descend')
               t = flip(t);
               didFlip = 1;
           elseif issorted(t,'ascend')
               didFlip = 0;
           else
               error('The requested points are not sorted in ascending order');
           end
           
           % number of collocation points
           N = length(t);
           M = length(t_knot);
           N_splines = M - K;
           iSpline = discretize(t_knot(K:(N-K+1)));
           x = zeros(size(t));
           indices = iSpline(~isnan(iSpline));
           
           
           delta_r = zeros(N,K);
           delta_l = zeros(N,K);
           b = zeros(N,K); b(:,1) = 1;
           
           for j=1:(K-1) % loop through splines of increasing order
               delta_r(indices,j) = t_knot(indices+j) - t(indices);
               delta_l(indices,j) = t(indices) - t_knot(indices+1-j);
               
               saved = zeros(N,1);
               for r=1:j % loop through the nonzero splines
                   term = b(indices,r)/(delta_r(indices,r) + delta_l(indices,j+1-r));
                   b(indices,r) = saved + delta_r(indices,r)*term;
                   saved = delta_l(indices,j+1-r)*term;
               end
               b(indices,j+1) = saved;
               
               indices = max(1,i-j):i;
               XB(t_i,indices,j+1) = b(1:length(indices));
           end
           
           for m = (j+1):(K-1) % coefficient loop
               ilo = K-j;
               for jj=1:(K-j)
                  A(jj) = (A(jj+1)*delta_l(ilo) + A(jj)*delta_r(jj))/(delta_l(ilo)+delta_r(jj));
                  ilo = ilo - 1;
               end
           end
        end
        
        function [t_knot] = NaturalKnotsForSpline( t, K, DF )
            %% NaturalKnotsForSpline
            %
            % Returns the natural knot points for splines of order K. It defaults to 1
            % degree of freedom, but you can override this by setting DF to some
            % nonnegative integer value.
            %
            % This matches the definitions in interpolation with tension
            % splines paper
            if nargin < 3
                DF = 1;
            end
            
            if (DF < 1 || mod(DF,1) ~= 0)
                disp('DF must be a non-negative integer');
                return;
            end
            
            t = [t(1); t(1+DF:DF:end-DF); t(end)];
            
            if mod(K,2) == 1
                % Odd spline order, so knots go in between points.
                dt = diff(t);
                
                % This gives us N+1 knot points
                t_knot = [t(1); t(1:end-1)+dt/2; t(end)];
                
                % Now remove start and end knots
                for i=1:((K-1)/2)
                    t_knot(2) = [];
                    t_knot(end-1) = [];
                end
                
            else
                t_knot = t;
                
                % Now remove start and end knots
                for i=1:((K-2)/2)
                    t_knot(2) = [];
                    t_knot(end-1) = [];
                end
                
            end
            
            % Now we increase the multiplicity of the knot points at the beginning and
            % the end of the interval so that the splines do not extend past the end
            % points.
            t_knot = [repmat(t_knot(1),K-1,1); t_knot; repmat(t_knot(end),K-1,1)];
        end
        
        function B = BSpline( t, t_knot, K, Kout )
            %% BSpline
            %
            % Returns the basis splines of order K evaluated at point t,
            % given knot points t_knot. If you optionally provide Kout,
            % only the Kout derivatives will be returned.
            
            if any(diff(t_knot)<0)
                error('t_knot must be non-decreasing');
            end
            
            if isempty(Kout) || Kout > K
                Kout = K;
            end
            
            S = K-1;

            t_knot2 = t_knot + [diff(t_knot); 0];
            
            % numer of knots
            M = length(t_knot);
            
            % This is true assuming the original t_knot was strictly monotonically
            % increasing (no repeat knots) and we added repeat knots at the beginning
            % and end of the sequences.
            N_splines = M - K;
            
            % number of collocation points
            N = length(t);
            
            % Rows are the N collocation points
            % Columns are the M splines
            B = zeros(N,N_splines,K); % This will contain all splines and their derivatives
            XB = zeros(N,N_splines,Kout); % This will contain all splines through order K
            for t_i=1:N % loop through all N collocation points
                i = find( t_knot <= t(t_i) & t(t_i) < t_knot2, 1, 'last' );
                if isempty(i)
                    if t(t_i) < t_knot(1)
                        %             i = find( t_knot == t_knot(1), 1, 'last');
                        continue; %This continue means we don't need to set b(1) = 0; or check indices on the delta_r line
                    elseif t(t_i) == t_knot(end)
                        i = find( t_knot < t(t_i), 1, 'last');
                    else
                        %             i = find( t_knot < t_knot(end), 1, 'last');
                        continue; %b(1) = 0;
                    end
                end
                
                delta_r = zeros(K,1);
                delta_l = zeros(K,1);
                
                XB(t_i,i,1) = 1;
                
                b = zeros(K,1); b(1) = 1;
                for j=1:(K-1) % loop through splines of increasing order
                    delta_r(j) = t_knot(i+j) - t(t_i);
                    delta_l(j) = t(t_i) - t_knot(i+1-j);
                    
                    saved = 0;
                    for r=1:j % loop through the nonzero splines
                        term = b(r)/(delta_r(r) + delta_l(j+1-r));
                        b(r) = saved + delta_r(r)*term;
                        saved = delta_l(j+1-r)*term;
                    end
                    b(j+1) = saved;
                    
                    indices = max(1,i-j):i;
                    XB(t_i,indices,j+1) = b(1:length(indices));
                end
                
                indices = max(1,i-K+1):i;
                B(t_i,indices,1) = b(1:length(indices));
                
            end
            
            diff_coeff = @(a,r,m) (K-m)*(a(2)-a(1))/(t_knot(r+K-m) - t_knot(r));
            
            for r=1:N_splines
                % alpha mimics equation X.16 in deBoor's PGS, but localized to avoid
                % the zero elements.
                alpha = zeros(S+2,S+2); % row is the coefficient, column is the derivative (1=0 derivatives)
                alpha(2,1) = 1;
                for m=1:(Kout-1) % loop over derivatives
                    for i=1:(m+1) % loop over coefficients
                        a = alpha(:,m);
                        alpha(i+1,m+1) = diff_coeff(a(i:end),r+i-1,m);
                        if isinf(alpha(i+1,m+1)) || isnan(alpha(i+1,m+1))
                            alpha(i+1,m+1) = 0;
                        end
                        if r+i-1>N_splines
                            B0 = zeros(N,1);
                        else
                            B0 = XB(:,r+i-1,K-m);
                        end
                        B(:,r,m+1) = B(:,r,m+1) + alpha(i+1,m+1)*B0;
                    end
                end
            end
        end
        
    end
end

