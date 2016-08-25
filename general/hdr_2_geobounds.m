function [geoExtent] = hdr_2_geobounds(hdr, align, bound)
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


%I'm not currently sure if bound = 'in' is useful but have added it
%nonetheless.

if regexpi(align, 'corner')
   %Do nothing
elseif regexpi(align, 'center')
    %Convert hdr to corner:
    hdr(3) = hdr(3) - 0.5*hdr(5);
    hdr(4) = hdr(4) - 0.5*hdr(5);
else
    error('Unknown ESRI header alignment.');
end
    
if regexpi(bound, 'out')
    lonW = hdr(3);
    lonE = hdr(3) + hdr(5)*(hdr(1) + 1);
    latS = hdr(4);
    latN = hdr(4) + hdr(5)*(hdr(2) + 1);
elseif regexpi(bound, 'in')
    lonW = hdr(3) + hdr(5);
    lonE = hdr(3) + hdr(5)*hdr(1);
    latS = hdr(4) + hdr(5);
    latN = hdr(4) + hdr(5)*hdr(2);  
    
elseif regexpi(bound, 'center')
    lonW = hdr(3) + 0.5*hdr(5);
    lonE = hdr(3) + hdr(5)*(hdr(1) + 0.5);
    latS = hdr(4) + 0.5*hdr(5);
    latN = hdr(4) + hdr(5)*(hdr(2) + 0.5);    
else
    error('Unknown bound option selected');
end
    
geoExtent = [lonW, lonE, latS, latN];
    
end