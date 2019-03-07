%%%% additional properties
% splines
N_max = 750 % maximum time series length
N_fits
t_boundary % boundaries defining which spline to use

%%%% from init

minOverlapTime = 1200;
minOverlapPoints = 2*self.K;
minOverlap = max(round(minOverlapTime/median(diff(self.t))),minOverlapPoints);

if 2*minOverlap>self.N_max
    self.N_max = 2*minOverlap;
    fprintf('Data is dense.\n');
end

self.N_fits = floor((N+minOverlap)/(self.N_max+minOverlap))+1;
N_overlaps = self.N_fits-1;
N_fit_length = floor((N - N_overlaps*minOverlap)/self.N_fits)+minOverlap;

fit_starts = [1; 1+(1:(self.N_fits-1))'*N_fit_length-minOverlap] ;
fit_ends = fit_starts+N_fit_length;
fit_ends(end) = N;

if self.N_fits > 1
    fprintf('Splitting into %d overlapping fits.\n',self.N_fits);
end

self.t_boundary = zeros(length(fit_starts)+1,1);
self.t_boundary(1) = self.t(1);
self.t_boundary(2:end-1) = self.t(fit_starts(2:end)) + (self.t(fit_ends(1:end-1))-self.t(fit_starts(2:end)))/2;
self.t_boundary(end) = self.t(end);

self.spline_x = cell(self.N_fits,1);
self.spline_y = cell(self.N_fits,1);
for iFit=1:self.N_fits
    indices = (fit_starts(iFit):fit_ends(iFit))';
    if self.shouldUseRobustFit == 1
        self.spline_x{iFit} = RobustTensionSpline(self.t(indices),self.q(indices),self.distribution,'K',self.K,'T',self.T);
        self.spline_y{iFit} = RobustTensionSpline(self.t(indices),self.r(indices),self.distribution,'K',self.K,'T',self.T);
    else
        self.spline_x{iFit} = TensionSpline(self.t(indices),self.q(indices),self.distribution,'K',self.K,'T',self.T);
        self.spline_y{iFit} = TensionSpline(self.t(indices),self.r(indices),self.distribution,'K',self.K,'T',self.T);
    end
end

function [x,y] = xyAtTime(self,time)
x=zeros(size(time));
y=zeros(size(time));
spline_id = discretize(time,self.t_boundary);
unique_spline_id = unique(spline_id);
for iSpline=1:length(unique_spline_id)
    indices = (spline_id == unique_spline_id(iSpline));
    x(indices) = self.spline_x{unique_spline_id(iSpline)}(time(indices));
    y(indices) = self.spline_y{unique_spline_id(iSpline)}(time(indices));
end
end