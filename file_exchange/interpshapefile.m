function value = interpshapefile(S, lat, lon, attribute)
%INTERPSHAPEFILE Interpolates data based on a shapefile data structure
%
% value = interpshapefile(S, lat, lon, attribute)
%
% This program determines the attribute value at the given points, based on
% the data in the geographic data structure S (created via
% shaperead).  
%
% If the specified attribute holds a numeric scalar, value will be a
% numeric array; points located outside the polygons will receive a NaN.
% Otherwise, value will be a cell array; points outside the polygons will
% hold an empty array.
%
% NOTE: This program assumes that the geographic region covered by two
% adjacent vertices is small enough that the lines connecting the polygon
% vertices can be approximated using straight lines in a cartesian plane,
% and that the region does not cross the poles.
%
% Input variables:
%   
%   S:              geographic data structure returned by a shaperead.  At
%                   minimum, the structure must have the Lat, Lon, and
%                   specified attribute fields.  If the Geometry field is
%                   present, the function uses only those elements
%                   designated Polygons; if the field is not present (e.g.
%                   if it was deleted in previous processing), it assumes
%                   all elements are polygons.
%
%   lat:            latitude coordinates of points to interpolate
%
%   lon:            longitude coordinates of points to interpolate
%
%   attribute:      name of attribute field to be interpolated
%
% Output variables:
%
%   value:          array of interpolated values corresponding to each
%                   lat/lon input.  If the specified attribute holds a
%                   numeric scalar, value will be a numeric array; points
%                   located outside the polygons will receive a NaN.
%                   Otherwise, value will be a cell array; points outside
%                   the polygons will hold an empty array.    
%
% Example:
% 
%   States = shaperead('usastatelo', 'UseGeoCoords', true);
%   latcorners = [max([States.Lat]) min([States.Lat])];
%   loncorners = [min([States.Lon]) max([States.Lon])];
%   [lat, lon] = track(latcorners, loncorners);
%   statename = interpshapefile(States, lat, lon, 'Name');
%   
%   axesm('mercator');
%   geoshow('usastatelo.shp', 'FaceColor', [.5 .5 1]);
%   plotm(lat, lon, 'r.');
%   textm(lat, lon, statename);
%
% See also:
%
%   shaperead

% Copyright 2005 Kelly Kearney

%------------------------------
% Check input
%------------------------------


%Check if shape structure has 'lon'/'lat' or 'X'/'Y' fields
strX = '';
strY = '';
if all(isfield(S, {'lon','lat'}))
    strX = 'lon';
    strY = 'lat';
elseif all(isfield(S, {'X','Y'}))
    strX = 'X';
    strY = 'Y';
elseif all(isfield(S, {'x','y'}))
    strX = 'x';
    strY = 'y';
else
    error('interpshapefile:unknownCrd','The coordinate fields could not be found.');
end


if ~isstruct(S) || ~isfield(S, attribute) || isempty(strX)
    error('The specified attribute field could not be found.');
end

if ~all(size(lat) == size(lon)) && any(size(lat) == 1) && any(size(lon) == 1)
    [lon, lat] = meshgrid(lon,lat);
elseif ~all(size(lat) == size(lon))
    error('interpshapefile:latLon','lat and lon must be matrices of the same size');
end

%------------------------------
% Remove any non-polygon 
% objects
%------------------------------


% If Geometry field is not present, all elements are assumed to be
% polygons.

if isfield(S, 'Geometry')
    ispolygon = strcmp('Polygon', {S.Geometry});
    S(~ispolygon) = [];
end

%------------------------------
% Find attribute value at each
% lat/lon point
%------------------------------


    

nclockwise = cellfun(@(x,y) length(find(ispolycw(x,y))), {S.(strX)}, {S.(strY)});
repVals = cellfun(@(nrep, val) repmat({val},1,nrep), num2cell(nclockwise), {S.(attribute)}, 'UniformOutput', false);
repVals = [repVals{:}];

[in, indexCell] = inpolygons(lon, lat, [S.(strX)], [S.(strY)]);

try
    index = cell2mat(indexCell);
catch
    warning('interpshapefile:overlap','Overlapping polygons; using first match for each point');
    index = zeros(size(indexCell));
    for icell = 1:numel(index)
        index(icell) = indexCell{icell}(1);
    end
end

isGood = (index ~= 0);
if isnumeric(S(1).(attribute)) && isscalar(S(1).(attribute))
    value = NaN(size(lon));
    value(isGood) = cell2mat(repVals(index(isGood)));
else
    value = cell(size(lon));
    value(isGood) = repVals(index(isGood));
end

