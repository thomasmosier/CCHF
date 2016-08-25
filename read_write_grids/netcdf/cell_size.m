function [cellSize, latGrid] = cell_size(hdr, r, pathWrite)
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


[~, latGrid] = geo_mesh(hdr);

%Check if grid exists and matches current header
if ~isempty(pathWrite) && exist(pathWrite, 'file')
   [cellSize, hdrLoad] = read_ESRI(pathWrite);
   if sum(hdrLoad == hdr) == 6
       return;
   end
end

cellSize = r^2*(pi/180*hdr(5))*(cosd(90-latGrid-0.5*hdr(5))-cosd(90-latGrid+0.5*hdr(5)));
    %`hdrDem(5)' is in degrees.  Must be in radians to calculate arc-length. 
    
if ~isempty(pathWrite)
	write_ESRI_v2(cellSize, hdr, pathWrite);
end
     
end