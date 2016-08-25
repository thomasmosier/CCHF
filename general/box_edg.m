function edge = box_edg(center)
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


%Orient data as column vector:
flip = 0;
if length(center(:,1)) > length(center(1,:))
    flip = 1;
    center = center';
end

if center(1) < center(2)
    dCenter = abs(diff(center));
    edge = [center(1) - 0.5*dCenter(1), ...
            center(1:end-1) + 0.5*dCenter, ...
            center(end) + 0.5*dCenter(end)];
elseif center(2) < center(1)
    dCenter = abs(diff(center));
    edge = [center(1) + 0.5*dCenter(1), ...
        center(1:end-1) - 0.5*dCenter, ...
        center(end) - 0.5*dCenter(end)];
else
    error('box_edg:repeatInd',['The centroid dataset has repeated ' ...
        'indexes.  This case has not been coded for.'])
    
end

%if data transposed, transpose back:
if flip == 1
   edge = edge';
end