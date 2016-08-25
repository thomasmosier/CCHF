function [hdrEsri, metaEsri] = read_wc_bin_hdr(pathWCfile)
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

fidWcHdr = fopen(pathWCfile,'r'); 
dimWC = textscan(fidWcHdr,'%s %f',2,'HeaderLines',2);   %[rows, colums]
boundsWC = textscan(fidWcHdr,'%s %f',5,'HeaderLines',6);    %[NoData, upper-left longitude (center-o-box), upper-left latitutde (center-o-box), dx, dy]
%Write ESRI format header:
%NCOLS      xxx
%NROWS      xxx
%XLLCORNER  xxx
%YLLCORNER  xxx
%CELLSIZE   xxx
%NODATA_VALUE

metaEsri = {'NCOLS';'NROWS';'XLLCORNER';'YLLCORNER';'CELLSIZE';'NODATA_VALUE'};

hdrEsri = [ dimWC{2}(2), ...     %ncols
            dimWC{2}(1), ...     %nrows
            boundsWC{2}(2) - boundsWC{2}(4)/2, ...  %xllcorner
            boundsWC{2}(3) - (dimWC{2}(1)-0.5)*boundsWC{2}(5), ...    %yllcorner
            boundsWC{2}(5), ...  %cellsize
            boundsWC{2}(1)]';    %NODATA_value
    clear dimWC boundsWC
fclose(fidWcHdr);
end