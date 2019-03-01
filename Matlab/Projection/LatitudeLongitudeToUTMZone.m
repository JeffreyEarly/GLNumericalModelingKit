function [x,y] = LatitudeLongitudeToUTMZone( lat, lon, zone)
%% LatitudeLongitudeToUTMZone
% You pick the zone and the standard UTM parameters will be chosen for your
% given latitude and longitude.
[x,y] = LatitudeLongitudeToTransverseMercator( lat, lon, (-183 + 6 * zone),0.9996 );
x = x + 500000; % False easting
end