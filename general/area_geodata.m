function area = area_geodata(lon,lat,sTyp)
% Copyright 2013, 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the GlobalClimateData Downscaling Package.
% 
% The GlobalClimateData Downscaling Package is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The Downscaling Package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%Calculate cartesian area of grid over the Earth.

%Inputs: 
%'lat' and 'lon' vectors of latitude and longitude.
%Optional argument indicates whether the input vectors indicate edges ('e') 
%or centroids ('c') of the grid.  Default is centroid.

%Output: 
%'area' = area of each lat/lon cell (units = 'r's units ^2), 

%Ensure that Lat is row vector
[rLat, cLat] = size(lat);
if rLat < cLat
    lat = lat';
end

if nargin == 2 || isempty(sTyp) || ~isempty(regexpi(sTyp,'c'))
    %Calculate half the difference 
    if numel(lat) > 1
        dLat = 0.5*diff(lat,1,1);
    elseif numel(lat) == 1
        dLat = 90;
    end
    if numel(lon) > 1
        dLon = 0.5*diff(lon,1,2);
    elseif numel(lon) == 1
        dLon = 180;
    end

    %Calculate bounding latitudes:
    if numel(lat) == 1
        latBS = lat - dLat;
        latBN = lat + dLat;
    elseif lat(1) > lat(2)
        latBS = lat - [dLat; dLat(end,:)];
        latBN = lat + [dLat(1,:); dLat];
    elseif lat(1) < lat(2)
        latBS = lat - [dLat(1,:); dLat];
        latBN = lat + [dLat; dLat(end,:)];
    else
       error('area_Geodata:lat','The latitudes in the georeferenced dataset are repeated.') 
    end

    %Ensure that Lon is column vector
    [rLon, cLon] = size(lon);
    if rLon > cLon && (rLon == 1 || cLon == 1)
        lon = lon';
    end
    
    %Calculate bounding longitudes:
    if numel(lon) == 1
        lonBW = lon - dLon;
        lonBE = lon + dLon;
    elseif lon(1) < lon(2)
        lonBW = lon - [dLon(:,1), dLon];
        lonBE = lon + [dLon, dLon(:,end)];
    elseif lon(1) > lon(2)
        lonBW = lon - [dLon, dLon(:,end)];
        lonBE = lon + [dLon(:,1), dLon];
    else
       error('area_Geodata:lon','The longitudes in the geo-referenced dataset are repeated.') 
    end
elseif ~isempty(regexpi(sTyp,'e'))
    latBN = lat(1:end-1);
    latBS = lat(2:end);
    lonBW = lon(1:end-1);
    lonBE = lon(2:end);
end

if any(size(lat) == 1)
    [xSW,ySW] = meshgrid(lonBW,latBS);
    [xNE,yNE] = meshgrid(lonBE,latBN);

    %Calculate fractional area over sphere:
    %Specify the reference sphere as Earth with units of 'meters' (output in
    %units of meters^2)
    area = single(areaquad(double(ySW),double(xSW),double(yNE),double(xNE),referenceSphere('earth','meters')));
else
    area = single(areaquad(double(latBS),double(lonBW),double(latBN),double(lonBE),referenceSphere('earth','meters')));
end
% area = areaquad(latBS,lonBW,latBN,lonBE,referenceSphere('earth','meters'));
