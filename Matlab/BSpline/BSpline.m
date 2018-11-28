classdef BSpline < handle
    %BSPLINE Summary of this class goes here
    %   3 argument initialization
    %       f = BSpline(K,t_knot,m);
    %   where
    %       K           order of the spline
    %       t_knot      knot points
    %       m           coefficients for the splines
    
    
    properties (Access = public)
        K       % order of polynomial        
        m       % spline coefficients (MxD)
        t_knot  % spline knot points
        
        t_pp    % pp break points. size(t_pp) = length(t_knot) - 2*K + 1
        C       % piecewise polynomial coefficients. size(C) = [length(t_pp)-1, K]
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = BSpline(K,t_knot,m)
            self.K = K;
            
            self.t_knot = t_knot;
            self.m = m;            
        end
        
        function varargout = subsref(self, index)
            %% Subscript overload
            %
            % The forces subscript notation to behave as if it is
            % evaluating a function.
            idx = index(1).subs;
            switch index(1).type
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '()'
                    if length(idx) >= 1
                        t = idx{1};
                    end
                    
                    if length(idx) >= 2
                        NumDerivatives = idx{2};
                    else
                        NumDerivatives = 0;
                    end
                    
                    varargout{1} = BSpline.EvaluateFromPPCoefficients(t,self.C,self.t_pp,NumDerivatives);
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '.'       
                    [varargout{1:nargout}] = builtin('subsref',self,index);
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '{}'           
                    error('The BSpline class does not know what to do with {}.');
                otherwise
                    error('Unexpected syntax');
            end
            
        end
        
        
        function x_out = ValueAtPoints( self, t, NumDerivatives)
            if ~exist('NumDerivatives','var')
                NumDerivatives = 0;
            end
            x_out = BSpline.EvaluateFromPPCoefficients(t,self.C,self.t_pp,NumDerivatives);
        end
        
        function set.m(self,newM)
            self.m = newM;
            self.splineCoefficientsDidChange();
        end
        
        function self = splineCoefficientsDidChange(self)
            [self.C,self.t_pp] = BSpline.PPCoefficientsFromSplineCoefficients( self.m, self.t_knot, self.K );
        end
    end
    
    
    methods (Static)
        function [C,t_pp] = PPCoefficientsFromSplineCoefficients( m, t_knot, K )
            %% PPCoefficientsFromSplineCoefficients
            % Returns the piecewise polynomial coefficients in matrix C
            % from spline coefficients in vector m.
            %
            % size(t_pp) = length(t_knot) - 2*K + 1
            % size(C) = [length(t_pp)-1, K]
            
            Nk = length(t_knot);
            t_pp = t_knot(K:(Nk-K));
            B = BSpline.Spline( t_pp, t_knot, K, K-1 );
            
            % Build an array of coefficients for polyval, highest order first.
            C = zeros(length(t_pp),K);
            for i=1:K
                C(:,K-i+1) = B(:,:,i)*m;
            end
            
            t_pp = t_knot(K:(Nk-K+1));
        end
        
        function f = EvaluateFromPPCoefficients(t,C,t_pp, D)
            %% EvaluateFromPPCoefficients
            %
            % Returns the value of the function with derivative D
            % represented by PP coefficients C at locations t. t_pp
            % contains the intervals.
            %
            % The returned array f is the same size as t.
            %
            %
            if issorted(t,'ascend')
                didFlip = 0;
            elseif issorted(t,'descend')
                t = flip(t);
                didFlip = 1;
            else
                error('Not sorted')
                [t,I] = sort(t);
                d = 1:length(t);
                returnIndices = d(I);
                didFlip = 2;
            end
            
            K = size(C,2);
            f = zeros(size(t));
            
            if nargin < 4
                D = 0;
            elseif D > K-1
                % By construction the splines are zero for K or more derivs
                return;
            end
            
            scale = factorial((K-1-D):-1:0);
            indices = 1:(K-D);
            
            % startIndex and endIndex are indices into t/f
%             startIndex = 1;
%             for i=2:length(t_pp)
%                 endIndex = find(t <= t_pp(i),1,'last');
%                 f(startIndex:endIndex) = polyval(C(i-1,indices)./scale,t(startIndex:endIndex)-t_pp(i-1));
%                 startIndex = endIndex+1;
%             end
%             f(startIndex:end) = polyval(C(i-1,indices)./scale,t(startIndex:end)-t_pp(i-1));
            
% The above implementation is faster, but doesn't actually work in some
% cases, like evaluating a single point. The discretize function is slow.
            t_pp_bin = discretize(t,[-Inf; t_pp(2:end-1); Inf]);
            startIndex = 1;
            while startIndex <= length(t)
                iBin = t_pp_bin(startIndex);
                endIndex = find( t_pp_bin == iBin, 1, 'last');
                f(startIndex:endIndex) = polyval(C(iBin,indices)./scale,t(startIndex:endIndex)-t_pp(iBin));
                startIndex = endIndex+1;
            end
            
            % include an extrapolated points past the end.
            if startIndex <= length(t)
                f(startIndex:end) = polyval(C(iBin,indices)./scale,t(startIndex:end)-t_pp(iBin));
            end
            
            if didFlip == 0
                return;
            elseif didFlip == 1
                f = flip(f);
            else
                f = f(returnIndices);
            end
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


