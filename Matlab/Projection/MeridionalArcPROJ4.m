function y = MeridionalArcPROJ4(phi)
% These are the *defined* values for WGS84
WGS84a=6378137;
WGS84invf=298.257223563;

% Constants
C00=1;
C02=.25;
C04=.046875;
C06=.01953125;
C08=.01068115234375;
C22=.75;
C44=.46875;
C46=.01302083333333333333;
C48=.00712076822916666666;
C66=.36458333333333333333;
C68=.00569661458333333333;
C88=.3076171875;

f = 1 / WGS84invf;
es = f*(2 - f);
t = es * es;
aPrime = C00 - es * (C02 + es * (C04 + es * (C06 + es * C08)));
bPrime = es * (C22 - es * (C04 + es * (C06 + es * C08)));
cPrime = t * (C44 - es * (C46 + es * C48));
t = t*es;
dPrime = t * (C66 - es * C68);
ePrime = t * es * C88;

sphi = sin(phi);
cphi = cos(phi);
cphi = cphi.*sphi;
sphi = sphi.*sphi;
y =  WGS84a*(aPrime * phi - cphi .* (bPrime + sphi.*(cPrime + sphi.*(dPrime + sphi.*ePrime))));
end