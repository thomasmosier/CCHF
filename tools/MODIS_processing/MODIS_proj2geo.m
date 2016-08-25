function [lon, lat] = MODIS_proj2geo(x, y)

%https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
%modis-land.gsfc.nasa.gov/MODLAND_grid.html

rEarth = 6371007.181; %units = meters
m2deg = 360/(2*pi*rEarth); %units = meters


lat = y*m2deg;
lon = x*m2deg./cosd(lat);


