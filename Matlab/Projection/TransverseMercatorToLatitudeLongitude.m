function [lat,lon] = TransverseMercatorToLatitudeLongitude(x, y, lon0, k0)
% These are the *defined* values for WGS84
WGS84a=6378137;
WGS84invf=298.257223563;

phi = InverseMeridionalArcPROJ4(y/k0);

if ( abs( phi ) >= 2*pi )
    if y<0
        lat = -90;
    else
        lat = 90;
    end
    lon = 0;
    return;
end

% Compute a few trig functions that we'll be using a lot
s = sin(phi); s2 = s.*s;
c = cos(phi); c2 = c.*c;
t = tan(phi); t2 = t.*t;

% Compute v and e2 -- pieces of this could be precomputed and #define
f = 1 / WGS84invf;
e2 = f*(2 - f);
v = WGS84a ./ sqrt( 1 - e2*s2);
e2 = e2 / (1 - e2); % From this point forward e2 will actually be (e^prime)^2
e2c2 = e2*c2;

T11 = 5 + 3*t2 + e2c2.*(1 - 4*e2c2 - 9*t2);
T12 = 61 + t2.*(90 - 25.*e2c2 + 45*t2) + 46*e2c2;
T13 = 1385 + t2.*(3633. + t2.*(4095 + t2.*1575.));

d = x./(v*k0);
d2 = d.*d;

phi = phi - 0.5*(1 + e2c2).*t.*d2.*(1 - (d2/12).*( T11 - (d2/30).*( T12 - (d2/56).*T13)));
lat = phi*180./pi;

T15 = 1 + 2*t2 + e2c2;
T16 = 5 + t2.*(28 + 24*t2 + 8*e2c2) + 6*e2c2;
T17 = 61 + t2 .* (662 + t2 .* (1320 + 720 * t2));

lon = lon0 + 180/pi*(d./c).*(1 - (d2/6).*( T15 - (d2/20).*( T16 - (d2/42).*T17 )));