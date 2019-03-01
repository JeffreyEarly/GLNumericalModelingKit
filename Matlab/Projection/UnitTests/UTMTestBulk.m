clear
tracks = open('raw_rho1_drifters.mat');

iDrifter = 1;
lat = tracks.lat{iDrifter};
lon = tracks.lon{iDrifter};

[x,y,zone] = LatitudeLongitudeToUTM(lat,lon);

[lat_back,lon_back] = UTMToLatitudeLongitude(x,y,zone);

max((lat-lat_back)./lat)
max((lon-lon_back)./lon)