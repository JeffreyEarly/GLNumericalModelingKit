function [x,y,zone] = LatitudeLongitudeToUTM( lat, lon)
%% LatitudeLongitudeToUTM
% Given a latitude and longitude, this returns the zone, easting and
% northing for the appropriate UTM projection.
zone = floor(lon / 6.0) + 31;
[x,y] = LatitudeLongitudeToUTMZone( lat, lon, zone );
end