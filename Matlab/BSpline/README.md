Interpolating splines, smothing splines and b-splines 
==============

A series of classes for interpolating and smoothing data using b-splines.

The `InterpolatingSpline` class is useful for interpolating between points when the data is not noisy, while the `SmoothingSpline` class is useful for smoothing noisy data. Both classes are subclasses of `BSpline`, which can be used to generate b-splines from any set of knot points. 

If you use these classes, please cite the following paper,
- Early, J. J., & Sykulski, A. M. (2020). [Smoothing and Interpolating Noisy GPS Data with Smoothing Splines](https://journals.ametsoc.org/view/journals/atot/37/3/JTECH-D-19-0087.1.xml), *Journal of Atmospheric and Oceanic Technology*, 37(3), 449-465.

There is a [five minute overview of the smoothing and interpolating noisy GPS data](https://jeffreyearly.com/smoothing-and-interpolating-noisy-gps-data/) paper.

### Table of contents
1. [Quick start](#quick-start)
2. [Basis spline](#basis-spline)
2. [Interpolating spline](#interpolating-spline)
2. [Smoothing spline](#smoothing-spline)
2. [GPS smoothing spline](#gps-smoothing-spline)

------------------------

Quick start
------------

### Interpolating spline

The `InterpolatingSpline` class is useful for interpolating smooth functions. Initialize with data (x,y),
```matlab
x = linspace(0,20,10)';
y = sin(2*pi*x/10);
spline = InterpolatingSpline(t,x)
```
and now you can evaluate the interpolated value,
```matlab
x_dense = linspace(0,20,100)';
y_dense = spline(x_dense);
```
By default the class uses cubic spline, but you can initialize with whatever order of spline you want, e.g.,
```matlab
spline = InterpolatingSpline(t,x,'K',5)
```

### Smoothing spline

If your data is noisy, you'll want to use the `SmoothingSpline` class instead. In this example we sample a smooth function (sine) and contaminate it Gaussian noise,
```matlab
N = 30;
L = 10;
sigma = 0.1;

x = linspace(0,L,N)';
y = sin(2*pi*x/L) + sigma*randn(N,1);
```
and then initialize the `SmoothingSpline` class with the data and the standard deviation of the noise,
```matlab
spline = SmoothingSpline(x,y,NormalDistribution(sigma))
```
That's it! The `spline` object can be evaluated at any point in the domain, just as with the interpolating spline class.

### GPS smoothing spline

The `GPSSmoothingSpline` class offers some GPS specific additions to the `SmoothingSpline` class including t-distributed errors and outlier detection.

If your data does not contain any obvious outliers, just initialize the  `GPSSmoothingSpline` with the time, latitude, and longitude,
```matlab
spline = GPSSmoothingSpline(t,lat,lon);
```
but if you suspect outliers, then you should add the additional flag to detect outliers,
```matlab
spline = GPSSmoothingSpline(t,lat,lon,'shouldUseRobustFit',1);
```
The smoothed and interpolated data can then be returned in either projected coordinates (x,y), or geographic coordinates (lat,lon),
```matlab
[x,y] = spline.xyAtTime(t);
[lat,lon] = spline.latLonAtTime(t);
```
The default projection used is a transverse Mercator projector with central meridian chosen to be the central longitude of the dataset. This value, as well as the false northing and false easting, can be overriden by setting flags at initialization.

Overview
------------

The `BSpline` class is a primitive class that creates b-splines given some set of knot points, and then evaluates the splines given some set of coefficients for the splines. This class is just for generating b-splines and doesn't do anything with data.

The following three classes all directly inherit from `BSpline`  class,

- An `InterpolatingSpline` uses local b-splines to interpolate between data points. The (K-1)th derivative of a K order spline is piecewise continuous. This is a generalization of Matlab's cubic spline interpolation function.
- The `ConstrainedSpline` class does a least-squares fit to given data points with a chosen set of b-splines. Constraints can be added at any point in the domain and the class can also accommodate non-gaussian errors.
- A `SmoothingSpline` can be used to smooth noisy data and attempt to recover the "true" underlying function.

The `BivariateSmoothingSpline` essentially takes  (t,x,y) as input data and creates smoothing splines for both x, y *after* removing a mean  `ConstrainedSpline` fit from both directions. The idea is to make the data stationary and therefore treat tension parameter optimization isotropically. 

The `GPSSmoothingSpline`  inherits from `BivariateSmoothingSpline` but takes latitude and longitude as arguments, and assume the noise follows a Student's t-distribution.

All of these classes include extensive documentation within the code itself.

Basis spline
------------

The `BSpline` class is a primitive class that creates b-splines given some set of knot points, and then evaluates the splines given some set of coefficients for the splines.

The class would be initialized as follows,
```matlab
K = 3; % order of spline
D = K-1; % number of derivates to return
t = (0:10)'; % observation points

% increase the multiplicity of the end knots for higher order splines
t_knot = [repmat(t(1),K-1,1); t; repmat(t(end),K-1,1)];

nSplines = length(t)+K-2;

% coefficients for the bsplines---set all of them to zero for now.
m = zeros(nSplines,1);

% initialize the BSpline class
B = BSpline(K,t_knot,m);
```

If you set the coefficient of one of the splines to 1, e.g.,
```matlab
m = zeros(nSplines,1);
m(3) = 1;
B.m = m;
```
Then you can plot that particular spline,
```matlab
tq = linspace(min(t),max(t),1000)';
figure, plot(tq,B(tq),'LineWidth',2)
```
Here's an image of all the splines and their first derivatives plotted,
<p align="center"><img src="figures/bspline.png" width="400" /></p>

**Note the usage here that calling `B(t)` evaluates the splines at points `t`.** This same notation also works for derivatives where  `B(t,n)` will return the `n`-th derivative the spline at points `t`.

This class serves as a useful building block for other classes.

Interpolating spline
------------

An interpolating spline uses local b-splines to interpolate across points. The (K-1)th derivative of a K order spline is piecewise continuous.

### Example

Let's start by defining a  function,
```matlab
L = 10;
f = @(x) sin(2*pi*x/L);
```
create a grid of observation points,
```matlab
N = 10;
x = linspace(0,2*L,N)';
```
and then sampling the function at the observation points
```matlab
y = f(x);
```
Let's create a dense grid in x just to better visualize the true function, and the sample points.
```matlab
x_dense = linspace(0,2*L,10*N)';

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
legend('true function', 'observations')
```
Finally, let's use an interpolating spline,
```matlab
spline = InterpolatingSpline(x,y);

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
plot(x_dense,spline(x_dense))
legend('true function', 'observations', 'spline fit')
```
This spline fit is identical to Matlab's built in griddedInterpolant with the 'spline' option,
```matlab
ginterp = griddedInterpolant(x,y,'spline');

residual = ginterp(x_dense)-spline(x_dense);
relative_error = max(abs(residual))/max(abs(y)) % returns O(1e-16)
```
However, unlike Matlab's implementation, this class lets you specify the order of the spline to something other than cubic (K=4). For example, we can do a K=2 fit (which is just linear interpolation), but also a (K=5) fit.
```matlab
spline_2 = InterpolatingSpline(x,y,'K',2);
spline_5 = InterpolatingSpline(x,y,'K',5);

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
plot(x_dense,spline_2(x_dense))
plot(x_dense,spline_5(x_dense))
legend('true function', 'observations', 'spline fit, K=2', 'spline fit, K=5')
```

<p align="center"><img src="figures/interpolatingspline.png" width="400" /></p>

### Options

The `InterpolatingSpline` class takes name/value pairs at initialization to set the spline order (or degree).

- `'K'` spline order, default is 4.
- `'S'` spline degree (order-1), default is 3.

Smoothing spline
------------

A smoothing spline can be used to smooth noisy data and attempt to recover the "true" underlying function.

### Example

Let's start by defining a new function,
```matlab
L = 10;
f = @(x) sin(2*pi*x/L);
```
create a grid of observation points,
```matlab
N = 30;
x = linspace(0,L,N)';
```
and then sampling the function at the observation points while adding Gaussian noise,
```matlab
sigma = 0.1;
y = f(x) + sigma*randn(N,1);
```
Let's create a dense grid in x just to better visualize the true function, and then plot the function and noisy observations,
```matlab
x_dense = linspace(0,L,10*N)';

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
legend('true function', 'noisy data')
```

<p align="center"><img src="figures/noisydata.png" width="400" /></p>

Finally, let's use a smoothing spline to try to smooth the data and plot the results,
```matlab
spline = SmoothingSpline(x,y,NormalDistribution(sigma));

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
plot(x_dense,spline(x_dense),'LineWidth',2)
legend('true function', 'noisy data', 'smoothing spline fit')
```
The `NormalDistribution` class indicates that the noise is expected to be normally distributed. In general, the arrgument must be a  [`Distribution`](../Distributions) subclass.

<p align="center"><img src="figures/noisydatawithtensionspline.png" width="400" /></p>

If you don't know the error distribution, you can use cross-validation (CV) to estimated the expected mean-square error and minimize that instead.
```matlab
spline.minimize( @(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromCV );
```
The result will still be senstive to which distribution you choose, e.g., a `SmoothingSpline` intialized with a t-distribution will produce different results than a `SmoothingSpline` initialized with a normal distribution, even when using cross-validation.

### Options

The `SmoothingSpline` class takes name/value pairs at initialization to set the spline order (or degree).

- `'K'` spline order, default is 4.
- `'S'` spline degree (order-1), default is 3.
- `'T'` tension degree, default is to use the same degree as the spline.
- `'lambda'` smoothing parameter, either pass a numeric value, or the `Lambda` enumeration. Default is `Lambda.optimalIterated`.
- `'mu'` mean tension.
- `'knot_dof'` number of degrees of freedom to be used in placing knot points. Either specify an integer, or `'auto'` to have the algorithm attempt to choose an appropriate number. Default is 1.
- `'t_knot'` array of manually placed knot points.

The `Lambda` enumeration has the following values,

- `optimalIterated`  which minimize the expected mean-square error, but may take a while for lots of data.
- `optimalExpected`  which takes a guess at minimizing the mean-square error based on the effective sample-size.
- `fullTensionExpected`  which takes a guess at the full tension solution assuming infinite effective sample size.

### Methods for minimization

A smoothing spline varies the smoothing parameter `smoothing' to minimize some penalty function, such as the mean square error or expected mean square error. The primary method for this is,
```matlab
[lambda,fval] = minimize(self,penaltyFunction)
```
which takes a `penaltyFunction` functional handle. This function handle *must* take a smoothing spline object as its argument (so it can adjust the smoothing parameter) and return a scalar value (presumably your measure of error). A simple example is given above with,
```matlab
spline.minimize( @(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromCV );
```

The class also includes several convenience functions that make use of the `minimization` function.
```matlab
[lambda,mse] = minimizeMeanSquareError(self,t_true,x_true) % minimize against true values
[lambda,mse] = minimizeExpectedMeanSquareError(self) % minimize using the expected MSE from the distribution.
[lambda,mse] = minimizeExpectedMeanSquareErrorInPercentileRange(self,pctmin,pctmax); % minimize the expected MSE, but only using points within given the CDF percentile range.
```
The method `minimizeExpectedMeanSquareErrorInPercentileRange` is remarkably good at excluding outliers.

Bivariate smoothing spline
------------

The `BivariateSmoothingSpline` assume that there is one indepent variable (time) and two isotropic dependent variables (x,y). The class encapsulates two properties, `spline_x` and `spline_y` which are just independent `SmoothingSpline` objects for the x and y data. The smoothing parameter `lambda` is always set isotropically and the noise is assumed to be isotropic as well.

The class is initialized with,
```matlab
spline = BivariateSmoothingSpline(t,x,y,distribution);
```
where the argument `distribution` must be a  [`Distribution`](../Distributions) subclass describing the expected distribution of errors.

Note that upon initialization, this class uses the `ConstrainedSpline` to remove a mean flow from (x,y), with the intention of isotropizing the velocity. It is not clear if this is always appropriate.

### Options

The `BivariateSmoothingSpline` class takes name/value pairs at initialization to set the spline order (or degree).

- `'K'` spline order, default is 4.
- `'S'` spline degree (order-1), default is 3.
- `'T'` tension degree , default is to use the same degree as the spline.
- `'gpsNoiseDistribution'` an instance of the `Distribution` class characterizing the GPS noise.
- `'shouldUseRobustFit'` either 1 or 0. Default is 0.

### Methods for outliers

These have no consequence on the fit, but can be useful for diagnostics.

- `outlierThreshold` [get/set] the distance at which something is considered an outlier. By default this is set using the CDF to find the value at 1-1/10000.
- `outlierIndices` [get] contains the indices of the outliers that were detected.
- `nonOutlierIndices` [get] contains the indices of points not considered outliers

This methods are especially useful for plotting and visualization, where you may want to highlight outliers in some fashion.



GPS smoothing spline
------------

The `GPSSmoothingSpline` class is useful for smoothing noisy gps data and removing outliers. The class is initialized with,
```matlab
spline = GPSSmoothingSpline(t,lat,lon);
```
where `t` is time, and `lat,lon` are the latitude and longitude. Internally, the latitude and longitude are *projected* using a transverse Mercator projection to positions in meters `x,y`. By default, the class will,

1. identify outliers if you pass `'shouldUseRobustFit',1` and,
2. smooth the data to the appropriate value.

The results can then be evaluated at any time,
```matlab
[x_smooth,y_smooth] = spline.xyAtTime(t);
```
or
```matlab
[lat_smooth,lon_smooth] = spline.latLonAtTime(t);
```

### Properties

The  `GPSSmoothingSpline` it inherits from `BivariateSmoothingSpline` and primarily adds built-in projection between geographic (lat,lon) and projected coordinates (x,y). It also assumes the GPS noise follows a t-distribution, `StudentTDistribution(8.5,4.5);`, although this can be overriden at initialization.

The `distanceError` property is root mean square of the spline errors under full tension. This property is only populated when looking for outliers.


### Options

The `GPSSmoothingSpline` class takes the same name/value pairs at initialization as `BivariateSmoothingSpline` above, as well as the following:

- `'lon0'` central meridian used for the transverse Mercator projection.
- `'x0'` false easting.
- `'y0'` false northing.
- `'gpsNoiseDistribution'` an instance of the `Distribution` class characterizing the GPS noise.

Note that if you are processing multiple GPS tracks with the intention of using projected coordinates (x,y), then you *must* set `lon0`, `x0`, and `y0` to be the same for all the tracks in order for them to share the same coordinate system. For example in my code to process multiple drifters, I grab these values from the first drifter, then apply them to the remaining drifters,
```matlab
if iDrifter == 1
        spline = GPSSmoothingSpline(t_drifter,drifters.lat{iDrifter},drifters.lon{iDrifter},'shouldUseRobustFit',1);
        [x,y] = spline.xyAtTime(t_interp);
        
        lon0 = spline.lon0;
        x0 = min(x);
        y0 = min(y);
        x = x-x0;
        y = y-y0;
    else
        spline = GPSSmoothingSpline(t_drifter,drifters.lat{iDrifter},drifters.lon{iDrifter},'shouldUseRobustFit',1,'lon0',lon0,'x0',x0,'y0',y0);
        [x,y] = spline.xyAtTime(t_interp(nonExtrapIndices));
    end
```
If you only ever use geographic coordinates (lat,lon), then this doesn't matter.
