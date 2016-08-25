function latInd = NC_lat_ind(ds)
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


latIndTemp = strcmpi(ds.variables,'latitude');

if sum(latIndTemp) == 1
    latInd = find(latIndTemp == 1); 
elseif sum(latIndTemp) == 0 %If no variables found or multiple matches, change search
    latIndTemp = strcmpi(ds.variables,'lat');
    
    if sum(latIndTemp) == 1
        latInd = find(latIndTemp == 1); 
    end
end
if sum(latIndTemp) > 1 || sum(latIndTemp) == 0
    %Allows user to choose which variable is correct
    latInd = NC_man_sel(latInd, ds, 'latitude');
end