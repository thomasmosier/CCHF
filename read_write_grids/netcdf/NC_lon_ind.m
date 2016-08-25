function lonInd = NC_lon_ind(ds)
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


lonInd = strcmpi(ds.variables,'lon');
if sum(lonInd) == 0
    lonInd = strcmpi(ds.variables,'longitude');
end
if sum(lonInd) == 0
   error('gcm_load:latFind',['The script was unable to detect which '...
       'variable of the GCM time-series file represents longitude.']); 
end