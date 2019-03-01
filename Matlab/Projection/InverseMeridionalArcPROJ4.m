function phi = InverseMeridionalArcPROJ4(y)
% These are the *defined* values for WGS84
WGS84a=6378137;
WGS84invf=298.257223563;

f = 1 / WGS84invf;
es = f*(2 - f);

k = 1/(1-es);
phi = y/WGS84a;

i = 0;
MAX_ITER = 10;
while (i < MAX_ITER)
    s = sin(phi);
    t = 1 - es * s .* s;
    t = (MeridionalArcPROJ4(phi) - y) .* (t .* sqrt(t)) * k / WGS84a;
    phi = phi - t;
    if max(abs(t)) < 1e-11
        break;
    end
    i = i+1;
end
