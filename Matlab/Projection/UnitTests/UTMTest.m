lat_expected = 45;
lon_expected = 45;

x_expected = 500000;
y_expected = 4982950.4;

[x,y,zone] = LatitudeLongitudeToUTM(lat_expected,lon_expected);

(x - x_expected)/x_expected
(y - y_expected)/y_expected

[lat,lon] = UTMToLatitudeLongitude(zone,x,y)