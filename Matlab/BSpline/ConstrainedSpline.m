classdef ConstrainedSpline < BSpline
    %ConstrainedSpline Summary of this class goes here
    %   2 argument initialization
    %       f = ConstrainedSpline(t,x,K,t_knot,constraints);
    %   where
    %       t               array of values for the independent axis
    %       x               array of values for the dependent axis 
    %       K               spline order
    %       t_knot          array of knot points
    %       constraints     constraints = struct('t',[],'D',[]);
    %       f               cubic spline interpolant
    % 
    % constraints are such that f^(D)(t)=0
    
    
    properties (Access = public)
        distribution
        t
        x
        
        CmInv
        X
        W
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = ConstrainedSpline(t,x,K,t_knot,distribution,constraints)
           % Make sure (t,x) are the right shapes, and see how many
            % dimensions of x we are given.
            t = reshape(t,[],1);
            N = length(t);
            x = reshape(x,[],1);
            if length(x) ~= N
                error('x and t must have the same length.');
            end
            
            if any(isnan(x)) || any(isnan(t))
                error('Some of the data contain NaNs!')
            end
            
            if isempty(distribution)
                distribution = NormalDistribution(1);
            end
            
            % terminate the splines at the boundaries
            nl = find(t_knot <= t_knot(1),1,'last');
            nr = length(t_knot)- find(t_knot == t_knot(end),1,'first')+1;
            t_knot = [repmat(t_knot(1),K-nl,1); t_knot; repmat(t_knot(end),K-nr,1)];
            
            % Compute the splines at the observation points
            X = BSpline.Spline( t, t_knot, K, 0 );
            
            % Now handle the various contraints
            if isempty(constraints)
                F=[];
                S=[];
            else
                % Deal with *local* constraints
                if ~isfield(constraints,'t') || ~isfield(constraints,'D')
                    F=[];
                else
                    M = size(X,2); % number of splines
                    NC = length(constraints.t); % number of constraints
                    if length(constraints.D) ~= NC
                        error('t and D must have the same length in the constraints structure.');
                    end
                    F = zeros(NC,M);
                    Xc = BSpline.Spline( constraints.t, t_knot, K, K-1 );
                    for i=1:NC
                        F(i,:) = squeeze(Xc(i,:,constraints.D(i)+1));
                    end
                end
                
                % Deal with *global* constraints
                if ~isfield(constraints,'global')
                    S = [];
                else
                    M = size(X,2); % number of splines
                    switch constraints.global
                        case ShapeConstraint.none
                            S = [];
                        case ShapeConstraint.monotonicIncreasing
                            S = tril(ones(M)); % positive lower triangle
                        case ShapeConstraint.monotonicDecreasing
                            S = -tril(ones(M)); % negative lower triangle
                            S(:,1)=1; % except the first point
                        otherwise
                            error('Invalid global constraint.');
                    end
                end
            end
   
            cachedVars.X = X;
                
            if isa(distribution,'NormalDistribution')
                [m,CmInv,cachedVars] = ConstrainedSpline.ConstrainedSolution(t,x,K,t_knot,distribution,F,[],S,cachedVars);
            else
                [m,CmInv,cachedVars] = ConstrainedSpline.IteratedLeastSquaresTensionSolution(t,x,t_knot,K,distribution,F,S,cachedVars);
            end
          
            self@BSpline(K,t_knot,m);
            self.distribution = distribution;
            self.t = t;
            self.x = x;
            self.CmInv = CmInv;
            self.X = X;
            self.W = cachedVars.W;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Smoothing matrix and covariance matrix
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function S = smoothingMatrix(self)
            % The smoothing matrix S takes the observations and maps them
            % onto the estimated true values.
            S = (self.X*(self.CmInv\(self.X.')))*self.W;
        end
    end
    
    
    methods (Static)
        
        function tc = MinimumConstraintPoints(t_knot,K,T)
            %% MinimumConstraintPoints
            % Given some set of knot points for a spline of order K, if you
            % want to universally constrain that spline at degree T, this
            % function returns a set of constraint locations (in t) that
            % let you do that, assuming repeat knot points at the beginning
            % and end.
            
            t = unique(t_knot);
            D = K-1-T; % 0 if we're constraining at the same order
            if mod(D,2) == 0
                ts = t(1) + (t(2)-t(1))/(D/2 + 2)*(1:(D/2+1)).';
                te = t(end-1) + (t(end)-t(end-1))/(D/2 + 2)*(1:(D/2+1)).';
                ti = t(2:end-2) + diff(t(2:end-1))/2;
                tc = cat(1,ts,ti,te);
            else
                ts = t(1) + (t(2)-t(1))/((D-1)/2+1)*(0:((D-1)/2)).';
                te = t(end) - (t(end)-t(end-1))/((D-1)/2 + 1)*(0:((D-1)/2)).';
                ti = t(2:end-1);
                tc = cat(1,ts,ti,te);
            end
        end
        
        function cachedVars = PrecomputeSolutionMatrices(t,x,K,t_knot,distribution,F,W,S,cachedVars)
            %% PrecomputeSolutionMatrices
            % Computes several cachable variables.
            if ~exist('cachedVars','var') || isempty(cachedVars)
                cachedVars = struct('t',t,'x',x,'t_knot',t_knot,'K',K,'distribution',distribution,'X',[],'XS',[],'rho_t',[],'W',W,'XWX',[],'XWx',[],'FS',[]);
            end
            
            if ~isfield(cachedVars,'X') || isempty(cachedVars.X)
                % These are the splines at the points of observation
                cachedVars.X = BSpline.Spline( t, t_knot, K, 0 ); % NxM
            end
            
            if ~isfield(cachedVars,'XS') || isempty(cachedVars.XS)
                if ~isempty(S)
                    cachedVars.XS = cachedVars.X*S; % NxM
                else
                    cachedVars.XS = cachedVars.X;
                end
            end
            
            if ~isfield(cachedVars,'FS') || isempty(cachedVars.FS)
                if ~isempty(S) && ~isempty(F)
                    cachedVars.FS = F*S; % NxM
                else
                    cachedVars.FS = F;
                end
            end
            
            XS = cachedVars.XS;
            N = length(x);
                        
            if ~isempty(distribution.rho) && (~isfield(cachedVars,'rho_t') || isempty(cachedVar.rho_t))
                cachedVars.rho_t = distribution.rho(t - t.');
            end
            
            if ~isfield(cachedVars,'W') || isempty(cachedVars.W)
                % The (W)eight matrix.
                %
                % For a normal distribution the weight matrix is the
                % covariance matrix.
                %
                % For anything other than a normal distribution, it is
                % *not* the covariance matrix. During IRLS this matrix will
                % be changing as points are reweighted.
                if isempty(W)
                    if isempty(distribution.rho)
                        W = 1/(distribution.sigma0)^2;
                    else
                        rho_t = cachedVars.rho_t;
                        sigma = ones(size(x))*distribution.sigma0;
                        Sigma2 = (sigma * sigma.') .* rho_t;
                        W = inv(Sigma2);
                    end
                end
                cachedVars.W = W;
            end
            
            if ~isfield(cachedVars,'XWX') || isempty(cachedVars.XWX)
                if size(W,1) == N && size(W,2) == N
                    XWX = XS'*W*XS;
                elseif length(W) == 1
                    XWX = XS'*W*XS;
                elseif length(W) == N
                    XWX = XS'*diag(W)*XS; % (MxN * NxN * Nx1) = Mx1
                else
                    error('W must have the same length as x and t.');
                end
                cachedVars.XWX = XWX;
            end
            
            if ~isfield(cachedVars,'XWx') || isempty(cachedVars.XWx)
                if size(W,1) == N && size(W,2) == N
                    XWx = XS'*W*x;
                elseif length(W) == 1
                    XWx = XS'*W*x;
                elseif length(W) == N
                    XWx = XS'*diag(W)*x; % (MxN * NxN * Nx1) = Mx1
                else
                    error('W must have the same length as x and t.');
                end
                cachedVars.XWx = XWx;
            end
                        
        end
        
        function [m,CmInv,cachedVars] = ConstrainedSolution(t,x,K,t_knot,distribution,F,W,S,cachedVars)
            %% ConstrainedSolution
            %
            % Returns the spline fit solution with constraints, given by
            % the F matrix. It also supports IRLS with non-constant sigma.
            %
            %
            % N     # of observations
            % M     # of splines
            % NC    # of constraints
            %
            % inputs:
            % t             observation times (Nx1)
            % x             observations (Nx1)
            % K             spline order
            % t_knot        array of knot points
            % distribution  expected distribution of the errors
            % F             constraints (NCxM)
            % W             weight matrix for observations (NxN)
            % S             global constraint matrix (MxM)
            % cachedVars    precomputed matrices for computation
            %
            % output:
            % m         coefficients of the splines, Mx1
            % CmInv     Inverse of the covariance of coefficients, MxM
            
            cachedVars = ConstrainedSpline.PrecomputeSolutionMatrices(t,x,K,t_knot,distribution,F,W,S,cachedVars);
            
            NC = size(F,1);
            XWX = cachedVars.XWX;
            XWx = cachedVars.XWx;
     
            % set up inverse matrices
            E_x = XWX; % MxM
            F_x = XWx;
            
            if isempty(S)
                % No global constraints
                if NC > 0
                    E_x = cat(1,E_x,F); % (M+NC)xM
                    E_x = cat(2,E_x,cat(1,F',zeros(NC)));
                    F_x = cat(1,F_x,zeros(NC,1));
                    m_x = E_x\F_x;
                    m = m_x(1:size(X,2));
                else
                    m = E_x\F_x;
                end
            else
                M = size(XWX,2); % number of splines
                lb = zeros(M,1); lb(1) = -inf;
                ub = inf*ones(M,1);
                
                % These are the local constraints, exactly as above
                Aeq = cachedVars.FS;
                if isempty(Aeq)
                    beq = [];
                else
                    beq = zeros(NC,1);
                end
                
                % E_x should be symmetric, although sometimes it's not
                % exactly.
                H = (E_x+E_x')*0.5;
                
                % Note, switching to 'trust-region-reflective' algorithm
                % requires initial conditions. It seems to work to simply
                % use xi0=S\m0, where m0 is the unconstrained solution from
                % above.
                options = optimoptions('quadprog','Algorithm','interior-point-convex');
                x = quadprog(2*H,-2*F_x,[],[],Aeq,beq,lb,ub,[],options);
                m = S*x;
            end
            CmInv = E_x;
            
            % Here's the other way to solve the constraint problem
            %                 xi = optimvar('xi',M,'LowerBound',0,'UpperBound',Inf);
            %                 xi(1).LowerBound = -Inf;
            %                 xi(1).UpperBound = Inf;
            %                 B = X*S;
            %
            %                 f = sum(x.^2) - 2*sum(x.*(B*xi)) + sum((B*xi).^2);
            %
            %                 qprob = optimproblem('Objective',f);
            %                 opts = optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off');
            %                 [sol,qfval,qexitflag,qoutput] = solve(qprob,struct('xi',xi0),'options',opts);
            %
            %                 m = S*sol.xi;
        end
        
        function [m,CmInv,cachedVars] = IteratedLeastSquaresTensionSolution(t,x,t_knot,K,distribution,F,S,cachedVars)
            % Out first call is basically just a normal fit to a tension
            % spline. If W is not set, it will be set to either 1/w0^2 or
            % the correlated version
            cachedVars.W = []; cachedVars.XWX = []; cachedVars.XWx = [];
            [m,CmInv,cachedVars] = ConstrainedSpline.ConstrainedSolution(t,x,K,t_knot,distribution,F,[],S,cachedVars);
            
            X = cachedVars.X;
            sigma2_previous = (distribution.sigma0)^2;
            rel_error = 1.0;
            repeats = 1;
            while (rel_error > 0.01)
                sigma_w2 = distribution.w(X*m - x);
                
                if exist('cachedVars.rho_t','var')
                    Sigma2 = (sqrt(sigma_w2) * sqrt(sigma_w2).') .* cachedVars.rho_t;
                    W = inv(Sigma2);
                else
                    W = diag(1./sigma_w2);
                end
                
                % hose any cached variable associated with W...
                cachedVars.W = []; cachedVars.XWX = []; cachedVars.XWx = [];
                % ...and recompute the solution with this new weighting
                [m,CmInv,cachedVars] = ConstrainedSpline.ConstrainedSolution(t,x,K,t_knot,distribution,F,W,S,cachedVars);
                
                rel_error = max( abs((sigma_w2-sigma2_previous)./sigma_w2) );
                sigma2_previous=sigma_w2;
                repeats = repeats+1;
                
                if (repeats == 250)
                    disp('Failed to converge after 250 iterations.');
                    break;
                end
            end
            
        end
  
    end
end


