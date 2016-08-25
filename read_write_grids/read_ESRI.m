function [dataESRI,hdrESRI,hdrMeta] = read_ESRI(pathESRI)
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


%Requires standard ESRI header and .asc or .txt file.
%ESRI Hdr format: 
    %line 1 = number of columns
    %line 2 = number of rows
    %line 3 = longitude, lower left corner (center)
    %line 4 = latitude, lower left corner (center)
    %line 5 = cellsize
    %line 6 = no data value

if isunix
    pathESRI = [filesep pathESRI];
end

fidESRI = fopen( pathESRI, 'r');
    hdrESRI   = textscan(fidESRI,'%s %f',6);
    hdrMeta = hdrESRI{1};
    hdrESRI = hdrESRI{2};
   
    dataESRIRaw = fscanf(fidESRI,'%f');
    fclose(fidESRI);
dataESRI = single(reshape(dataESRIRaw, hdrESRI(1), hdrESRI(2))'); 
    %When using 'reshape', swap row and column indices then transpose.
    dataESRI(dataESRI == hdrESRI(6)) = NaN;
    clear('dataESRIRaw');
end