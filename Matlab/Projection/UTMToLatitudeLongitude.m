function [lat,lon] = UTMToLatitudeLongitude( x, y, zone)
x = x-500000; % False easting
[lat, lon] = TransverseMercatorToLatitudeLongitude( x, y, 6 * zone - 183, 0.9996);
