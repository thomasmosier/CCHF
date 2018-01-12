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

function d = haversine_reference(lon0,lat0,lonVec,latVec,r)

szGrid = size(lonVec);
    
if ~any(size(lonVec) == 1) && ~any(size(latVec) == 1) && isequal(size(latVec),size(lonVec))
    lonVec = reshape(lonVec,[],1);
    latVec = reshape(latVec,[],1);
elseif any(size(lonVec) == 1) && ~any(size(latVec) == 1) || ~any(size(lonVec) == 1) && any(size(latVec) == 1)
    error('haversine:matrixAndVec','One input appears to be a vector and the other a matrix.  Both inputs must have same form.')
end

lon0 = repmat(lon0,numel(lonVec),1);
lat0 = repmat(lat0,numel(latVec),1);

d = 2*r*asin(sqrt(sind(0.5*(lat0-latVec)).^2 + cosd(lat0).*cosd(latVec).*sind(0.5*(lonVec-lon0)).^2));

if all(szGrid ~= 1)
   d = reshape(d,szGrid); 
end
