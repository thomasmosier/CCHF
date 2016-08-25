function [longW, longE, latS, latN] = ESRI_hdr_bounds(ESRIhdr)
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


%Reads ESRI formatted header (First 6 rows of ESRI ASCII file) and returns
%a vector of length 4 that includes the geographic bounding box in decimal
%degrees at the center of the exterior-most rows/columns within the grid.

%ESRII Header Format:
% NCOLS xxx
% NROWS xxx
% XLLCORNER xxx
% YLLCORNER xxx
% CELLSIZE xxx
% NODATA_VALUE xxx

longW = ESRIhdr(3) + 0.5*ESRIhdr(5);
longE = longW + ESRIhdr(5)*(ESRIhdr(1)-1);
latS = ESRIhdr(4) + 0.5*ESRIhdr(5);
latN = latS + ESRIhdr(5)*(ESRIhdr(2)-1);
end