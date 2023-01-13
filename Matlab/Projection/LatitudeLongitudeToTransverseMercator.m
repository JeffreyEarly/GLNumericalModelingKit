function [x,y] = LatitudeLongitudeToTransverseMercator(lat, lon, options)
arguments
    lat (:,1) double {mustBeNumeric,mustBeReal}
    lon (:,1) double {mustBeNumeric,mustBeReal}
    options.lon0 (1,1) double {mustBeNumeric,mustBeReal}
    options.k0 (1,1) double {mustBeNumeric,mustBeReal} = 0.9996;
end
k0 = options.k0;
if ~isfield(options.lon0)
    lon0 = min(lon) + (max(lon)-min(lon))/2;
else
    lon0 = options.lon0;
end

% These are the *defined* values for WGS84
WGS84a=6378137;
WGS84invf=298.257223563;

% Convert from degrees to radians
phi = lat*pi/180;

% Compute a few trig functions that we'll be using a lot
s = sin(phi);
c = cos(phi);
t = tan(phi);
s2 = s.*s;
c2 = c.*c;
t2 = t.*t;

% Compute v and e2 -- pieces of this could be precomputed and #define
f = 1 / WGS84invf;
e2 = f*(2 - f);
v = WGS84a ./ sqrt( 1 - e2*s2);
e2 = e2 / (1 - e2); % From this point forward e2 will actually be (e^prime)^2
e2c2 = e2*c2;

deltaLambda = (lon - lon0)*pi/180;
d2c2 = deltaLambda.*deltaLambda.*c2;

% Terms to compute x.
T7 = 1 - t2 + e2c2;
T8 = 5 +t2.*(t2- 18) + e2c2.*(14 - 58*t2); % + 13.*e4*c4 + 4.*e6*c6 - 64.*t2*e4*c4 - 24.*t2*e6*c6;
T9 = 61 - t2.*(479 - t2.*(179 - t2));

x = k0 .* v .* c .* deltaLambda .* (1 + (d2c2/6).*( T7 + (d2c2/20).*( T8 + (d2c2/42).*T9 )));

% Terms to compute y.
T3 = 5 - t2 + e2c2.*(9 + 4*e2c2);
T4 = 61 - t2.*(58 - t2) + 270*e2c2 - 330*t2.*e2c2; % + 445.*e4*c4 + 324.*e6*c6 - 680.*t2*e4*c4 + 88.*e8*c8 - 600.*t2*e6*c6 - 192.*t2*e8*c8;
T5 = 1385 - t2.*(3111 - t2.*(543 - t2));

y = k0 * MeridionalArcPROJ4(phi) + (k0 * v .* s .* c / 2) .* deltaLambda.*deltaLambda .* (1 + (d2c2/12).*( T3 + (d2c2/30).*( T4 + (d2c2/56.).*T5)));
end