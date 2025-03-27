RealToComplexTransform.makeMexFiles();

%%
N=32; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/N;

% Create a simple function, x = 2 + 3*sin(4*pi*t);
% The lowest resolved frequency is one cycle per unit time
x = 1*ones(size(t)) + 3*sin(2*(2*pi)*t);

dft = RealToComplexTransform(size(x),dims=1,nCores=8,planner="measure");

%%

xbar = dft.transformForward(x);
xback = dft.transformBack(xbar);
xback2 = double(zeros(dft.realSize));
xback2 = dft.transformBackIntoArray(xbar,xback2);
