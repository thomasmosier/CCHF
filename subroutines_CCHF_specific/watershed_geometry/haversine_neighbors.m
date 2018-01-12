% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function d = haversine_neighbors(lon,lat,r,varargin)


if any(size(lon) == 1) && any(size(lon) == 1)
    [lon, lat] = meshgrid(lon,lat);
elseif any(size(lon) == 1) && ~any(size(lon) == 1) || ~any(size(lon) == 1) && any(size(lon) == 1)
    error('haversine:matrixAndVec','One input appears to be a vector and the other a matrix.  Both inputs must have same form.')
end

[ic, id] = ixneighbors(lon);

d = 2*r*asin(sqrt(sind(0.5*(lat(ic)-lat(id))).^2 + cosd(lat(ic)).*cosd(lat(id)).*sind(0.5*(lon(id)-lon(ic))).^2));

if ~isempty(varargin(:)) 
    if regexpbl(varargin{1},'sparse')
        d = sparse(ic, id, d, numel(lon), numel(lon));
    elseif regexpbl(varargin{1},'full')
        d = full(sparse(ic, id, d, numel(lon), numel(lon)));
    else
        error('haversine:outputType',['Output format ' varargin{1} ' is not recognized.'])
    end
end