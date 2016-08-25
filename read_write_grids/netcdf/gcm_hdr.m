function hdr = gcm_hdr(gcmLat, gcmLon)
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


%Designed to be used in conjunction with 'gcm_load' to create an ESRI
%format header. Note that it assumes uniform lat and lon spacing.

szLon = size(gcmLon);
if szLon(1) < szLon(2)
    gcmLon = gcmLon';
end

gcmStep = mean([abs(gcmLat(1:end-1) - gcmLat(2:end));abs(gcmLon(1:end-1) - gcmLon(2:end))]);

if gcmLat(1) > gcmLat(end)
    yll = gcmLat(end)-0.5*gcmStep;
else
    yll = gcmLat(1)-0.5*gcmStep;
end

hdr = [ length(gcmLon); ...
        length(gcmLat); ...
        gcmLon(1)-0.5*gcmStep; ...
        yll; ...
        gcmStep; ...
        -9999 ];
end