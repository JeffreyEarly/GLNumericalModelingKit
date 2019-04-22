Distributions 
==============

A series of classes for modeling distributions and stochastic processes.

If you use these classes, please cite the following paper,
- J. Early and A. Sykulski. Smoothing and interpolating noisy GPS data with smoothing splines. IEEE Transactions on Signal Processing. In prep.

### Table of contents
1. [Quick start](#quick-start)
2. [Distribution class](#distribution-class)


------------------------

Quick start
------------

The `Distribution` class is an abstract class, so in practice one only uses its subclasses such as `NormalDistribution` directly. For example, to create a Gaussian distribution,
```matlab
sigma = 1;
distribution = NormalDistribution(sigma);
```
and now you can create some random values,
```matlab
n = 1000;
z = distribution.rand([n 1]);
```
and plot those values,
```matlab
figure
histogram(z,'Normalization','pdf')
```
which can be compared against the expected pdf,
```matlab
zdist = linspace(min(z),max(z),100)';
hold on
plot(zdist,distribution.pdf(zdist))
```
<p align="center"><img src="figures/normaldistribution.png" width="400" /></p>

All distribution subclasses include the `cdf` as a function, as well as the total variance of the process.

### Correlated values

It is also sometimes useful to correlate the random values in order to create a stochastic process. Let's define an autocorrelation sequence,
```matlab
tau = 10;
distribution.rho = @(z) exp(-(z/tau).^2);
```
and now generate a signal,
```matlab
t = (0:500).';
z = distribution.noise(t);
```
<p align="center"><img src="figures/correlatednoise.png" width="400" /></p>

The signal now looks smooth on time scales less than O(tau).


Overview
------------

The `Distribution` class is an abstract class


Distribution class
------------



