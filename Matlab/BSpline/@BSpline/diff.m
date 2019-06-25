function dspline = diff(spline,n)
%DIFF Differentiation of a BSpline

if nargin < 2
    n = 1;
end

if n == 0
    dspline = spline;
elseif n >= spline.K
    dspline = BSpline(spline.K,spline.t_knot,zeros(size(spline.m)));
    n = spline.K-1;
    dspline.K = spline.K-n;
    dspline.B = spline.B(:,:,1:(spline.K-n));
    dspline.C = zeros(size(spline.C(:,1:(spline.K-n))));
elseif ( ( n > 0 ) && ( round(n) == n ) )   % Positive integer
    %     dspline = BSpline(spline.K-1,spline.t_knot(2:end-1),spline.m(1:end-1));
    %     dspline.K = spline.K-n;
    %     dspline.B = spline.B(:,:,(n+1):(spline.K-n));
    %     dspline.C = spline.C(:,1:(spline.K-n));
    
    if 1 == 0
        ts = BSpline.PointsOfSupport(spline.t_knot,spline.K,0);
        t_knot = spline.t_knot(2:end-1);
        X = BSpline.Spline(ts,t_knot,spline.K-1);
        
        g = spline.ValueAtPoints(ts,1);
        m = X\g;
        dspline = BSpline(spline.K-1,t_knot,m);
    else
        D = n;
        m = spline.m;
        K = spline.K;
        t_knot = spline.t_knot;
        
        alpha = zeros(length(m),K+1);
        alpha(:,1) = m; % first column is the existing coefficients
        for r=1:length(m) % loop over splines
            for d=1:D % loop over derivatives
                alpha(r,1+d) = (K-d)*(alpha(r,d)-alpha(r,d))/(t_knot(r+K-d) - t_knot(r));
            end
        end
        
        % "a" are coefficients, "r" is the coefficient index, "d"
        % derivative
        diff_coeff = @(a,r,d) (K-d)*(a(2)-a(1))/(t_knot(r+K-d) - t_knot(r));
        
        for r=1:length(m)
            % alpha mimics equation X.16 in deBoor's PGS, but localized to avoid
            % the zero elements.
            alpha = zeros(K+1,K+1); % row is the coefficient, column is the derivative (1=0 derivatives)
            alpha(2,1) = 1;
            for d=1:D % loop over derivatives
                for i=1:(d+1) % loop over coefficients
                    a = alpha(:,d);
                    alpha(i+1,d+1) = diff_coeff(a(i:end),r+i-1,d);
                    if isinf(alpha(i+1,d+1)) || isnan(alpha(i+1,d+1))
                        alpha(i+1,d+1) = 0;
                    end
                    if r+i-1>N_splines
                        B0 = zeros(N,1);
                    else
                        B0 = XB(:,r+i-1,D+1-d); % want the K-m order spline, in position D+1-m
                    end
                    B(:,r,d+1) = B(:,r,d+1) + alpha(i+1,d+1)*B0;
                end
            end
        end
        
        
        
    end
else
    error('Can only differentiate with positive integers');
end

