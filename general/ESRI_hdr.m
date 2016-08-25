function hdr = ESRI_hdr(vecX, vecY, strOffset)
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


%ESRI Hdr format: 
    %line 1 = number of columns
    %line 2 = number of rows
    %line 3 = longitude, lower left corner (center or corner)
    %line 4 = latitude, lower left corner (center or corner)
    %line 5 = cellsize
    %line 6 = no data value

%If 'vecX' or 'vecY' are actually matrices, convert them to vectors:
if sum(size(vecX) == 1) == 0
    vecX = vecX(1,:);
end

if sum(size(vecY) == 1) == 0
    vecY = vecY(:,1);
end
    
%Check that ordering of vectors is correct:
if vecX(1) > vecX(2) 
    vecX = fliplr(vecX);
end

if vecY(2) > vecY(1) 
    vecY = flipud(vecY);
end

%Ensure that lon vector is column and lat vector is row vector:
vecX = reshape(vecX,1,[]);
vecY = reshape(vecY,[],1);

%Calculate the average step size:
step = mean(abs([diff(reshape(vecX,1,[])), diff(reshape(vecY,1,[]))]));
if round(100*(step - abs(vecY(2) - vecY(1)))/step) ~= 0 || round(100*(step - abs(vecX(1) - vecX(2)))/step) ~= 0
    warning('ESRI_hdr:avgStep',['An average stepsize is being used in '...
        'defining the header; however, the stepsize is non-uniform.'])
end

%Create Header:
hdr = NaN(6,1);
hdr(1) = length(vecX(:));
hdr(2) = length(vecY(:));
if ~isempty(regexpi(strOffset,'cor'))
    hdr(3) = vecX(1) - 0.5*step;
    hdr(4) = vecY(end) - 0.5*step;
elseif ~isempty(regexpi(strOffset,'cen'))
    hdr(3) = vecX(1);
    hdr(4) = vecY(end);  
else
    error('ESRI_hdr:offset',[strOffset ' is an unknown ESRI header type.']);
end
hdr(5) = step;
hdr(6) = -9999;