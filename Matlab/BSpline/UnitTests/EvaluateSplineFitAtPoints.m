f = @(x) sin(2*pi*x/10);

K = 3; % order of spline
D = K-1; % number of derivates to return
t = (0:10)'; % observation points
t_knot = BSpline.KnotPointsForPoints(t,K);

B = BSpline.Spline( t, t_knot, K, D );

% Now find the coefficients of the spline fit
X = B(:,:,1);
m = X\f(t);

% Fill in the gaps with extra points
tq = linspace(min(t),max(t),10000)';
Bq = BSpline.Spline( tq, t_knot, K, D );
Xq = Bq(:,:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this first test, we just use the matrix form of the splines.
figure
subplot(K,1,1)
scatter(t,f(t)), hold on
plot(tq,Xq*m)
title('Evaluated with Splines')
for iPlot = 1:D
    subplot(K,1,iPlot+1)
    plot(tq,Bq(:,:,iPlot+1)*m)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this second test, we use the PP coefficients to evaluate the fit, and
% also extend the end points to test extrapolation.

[C,t_pp] = BSpline.PPCoefficientsFromSplineCoefficients( m, t_knot, K );

tq = linspace(min(t)-1,max(t)+1,10000)';

figure
subplot(K,1,1)
scatter(t,f(t)), hold on
plot(tq,BSpline.EvaluateFromPPCoefficients(tq,C,t_pp))
title('Evaluated with PPs')
for iPlot = 1:D
    subplot(K,1,iPlot+1)
    plot(tq,BSpline.EvaluateFromPPCoefficients(tq,C,t_pp,iPlot))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% In this third test, we the class/object based interface.

f_spline = BSpline(t,f(t),K);

figure
subplot(K,1,1)
scatter(t,f(t)), hold on
plot(tq,f_spline(tq) )
title('Evaluated with PPs via OO interface')
for iPlot = 1:D
    subplot(K,1,iPlot+1)
    plot(tq,f_spline(tq,iPlot))
end

profile on
for i=1:10000
    f = BSpline.EvaluateFromPPCoefficients(tq,C,t_pp);
end
profile viewer