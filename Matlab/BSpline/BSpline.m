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
           %
           % This is the function bvalue.f from PGS.
           % http://pages.cs.wisc.edu/~deboor/pgs/bvalue.f
           
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
           
           if didFlip == 1
               x = flip(x);
           end
        end
        
        function t_knot = KnotPointsForPoints( t, K, DF )
            %% KnotPointsForPoints
            %
            % Returns the natural knot points for splines of order K for
            % observations at time points t. These knot points are place so
            % that each spline has exactly enough support for a fully
            % constrained system.
            %
            % It defaults to 1 degree of freedom, but you can override this
            % by setting DF to some nonnegative integer value.
            %
            % This matches the definitions in interpolation with tension
            % splines paper.
            
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
        
        function B = Spline( t, t_knot, K, D )
            %% Spline
            %
            % Returns the basis splines of order K evaluated at point t,
            % given knot points t_knot. If you optionally provide D,
            % then D derivatives will be returned.
            %
            % size(B) = [length(t) M D] where M=length(t_knot)-K is the
            % number of splines
            
            if any(diff(t_knot)<0)
                error('t_knot must be non-decreasing.');
            end
            
            if nargin < 4
                D = 0;
            elseif D > K-1
                D = K-1;
            end
            
            % numer of knots
            M = length(t_knot);
            
            % This is true assuming the original t_knot was strictly monotonically
            % increasing (no repeat knots) and we added repeat knots at the beginning
            % and end of the sequences.
            N_splines = M - K;
            
            % number of collocation points
            N = length(t);
            
            % 1st index is the N collocation points
            % 2nd index is the the M splines
            B = zeros(N,N_splines,D+1); % This will contain all splines and their derivatives
            delta_r = zeros(N,K);
            delta_l = zeros(N,K);
            knot_indices = discretize(t,t_knot(1:(M-K+1)));
            
            % XB will contain all splines from (K-D) through order (K-1).
            % These are needed to compute the derivatives of the spline, if
            % requested by the user.
            %
            % The indexing is such that the spline of order m, is located
            % at index m-(K-D)+1
            if D > 0
                XB = zeros(N,N_splines,D);
                if D + 1 == K % if we go tho the max derivative, we need to manually create the 0th order spline.
                    to_indices = ((1:N)' + size(XB,1)*( (knot_indices-1) + size(XB,2)*(1-1)));
                    XB(to_indices) = 1;
                end
            end
            
            b = zeros(N,K);
            b(:,1) = 1;
            for j=1:(K-1) % loop through splines of increasing order: j+1
                delta_r(:,j) = t_knot(knot_indices+j) - t;
                delta_l(:,j) = t - t_knot(knot_indices+1-j);
                
                saved = zeros(N,1);
                for r=1:j % loop through the nonzero splines
                    term = b(:,r)./(delta_r(:,r)+delta_l(:,j+1-r));
                    b(:,r) = saved + delta_r(:,r).*term;
                    saved = delta_l(:,j+1-r).*term;
                end
                b(:,j+1) = saved;
                
                % Save this info for later use in computing the derivatives
                % have to loop through one index.
                if j+1 >= K-D && j+1 <= K-1 % if K-j == 1, we're at the end, j+1 is the spline order, which goes into slot j+1-(K-1-D). Thus, if j+1=K, and D=1, this goes in slot 2
                    for r = 1:(j+1)
                        % (i,j,k) = (1,knot_indices-j+(r-1),j+1) --- converted to linear indices
                        to_indices = ((1:N)' + size(XB,1)*( (knot_indices-j+(r-1)-1) + size(XB,2)*(j+1-(K-1-D)-1)));
                        from_indices = (1:N)' + size(b,1) * (r-1);
                        XB(to_indices) = b(from_indices);
                    end
                end
            end
            
            for r = 1:K
                to_indices = ((1:N)' + size(B,1) * ( (knot_indices-(K-1) + (r-1)-1) + size(B,2)*(1-1) ));
                from_indices = (1:N)' + size(b,1) * (r-1);
                B(to_indices) = b(from_indices);
            end
            
            diff_coeff = @(a,r,m) (K-m)*(a(2)-a(1))/(t_knot(r+K-m) - t_knot(r));
            
            if D > 0
                for r=1:N_splines
                    % alpha mimics equation X.16 in deBoor's PGS, but localized to avoid
                    % the zero elements.
                    alpha = zeros(K+1,K+1); % row is the coefficient, column is the derivative (1=0 derivatives)
                    alpha(2,1) = 1;
                    for m=1:D % loop over derivatives
                        for i=1:(m+1) % loop over coefficients
                            a = alpha(:,m);
                            alpha(i+1,m+1) = diff_coeff(a(i:end),r+i-1,m);
                            if isinf(alpha(i+1,m+1)) || isnan(alpha(i+1,m+1))
                                alpha(i+1,m+1) = 0;
                            end
                            if r+i-1>N_splines
                                B0 = zeros(N,1);
                            else
                                B0 = XB(:,r+i-1,D+1-m); % want the K-m order spline, in position D+1-m
                            end
                            B(:,r,m+1) = B(:,r,m+1) + alpha(i+1,m+1)*B0;
                        end
                    end
                end
            end
            
        end
        
    end
end

