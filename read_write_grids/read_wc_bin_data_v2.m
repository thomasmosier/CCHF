function [dataMat, hdrClip, bndsDdWc, bndsIndWc] = read_wc_bin_data_v2(pathWcData, hdrRaw, gridDd, metVar, cropSide)    
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

%Loads WorldClim binary file, gridded to the entire world, and clips to
%region specified by 'gridBoundsDd', which is a vector of length four with
%the format [longW, longE, latS, latN].


fidWcData = fopen( pathWcData ); %Open binary WC file. 

wcVec = fread(fidWcData,'*int16');  %Read WC binary file.
dataMat = single(reshape(wcVec, hdrRaw(1), hdrRaw(2) )');    %After transpose, m (1st arg) is the number of columns and n (2nd arg) the number of rows.
    
if ~isempty(gridDd)
    [bndsDdWc, bndsIndWc] = adj_bounds(hdrRaw, gridDd, cropSide);    %Find indices to crop WC data to.
    dataMat = dataMat(bndsIndWc(4):bndsIndWc(3), bndsIndWc(1):bndsIndWc(2));  %Crop WC to desired region.
else
    bndsDdWc = [hdrRaw(3) + hdrRaw(5)/2, hdrRaw(3) + hdrRaw(5)*(0.5 + numel(dataMat(1,:))), hdrRaw(4) - hdrRaw(5)/2, hdrRaw(4) - (numel(dataMat(:,1)) + 0.5)*hdrRaw(5)];
    bndsIndWc = [1, numel(dataMat(1,:)), 1, numel(dataMat(:,1))];
end
    
dataMat( dataMat == hdrRaw(6) ) = NaN;
 
hdrClip(1) = length(dataMat(1,:));
hdrClip(2) = length(dataMat(:,1));
if ~isempty(gridDd)
    hdrClip(3) = bndsDdWc(1) - hdrRaw(5)/2;
    hdrClip(4) = bndsDdWc(3) - hdrRaw(5)/2;
else
    hdrClip(3) = hdrRaw(3);
    hdrClip(4) = hdrRaw(4);
end
hdrClip(5) = hdrRaw(5);
hdrClip(6) = hdrRaw(6);

if regexpbl(metVar,{'tmn','tmean','tmp','tmx'}) 
    dataMat = dataMat / 10;
elseif regexpbl(metVar,{'pr'}) || regexpbl(metVar,{'alt','elev','DEM'})
    %Do nothing
else
    warning('WorldClim_read:metVar',['An unknown meteorological '...
        'variable,' metVar ', has been used.  This may need to be '...
        'scaled; however, currently is not.  Add code to accomodate this']);
end
