% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function slope = slope_2pts(deltaH,pt1,pt2)
%deltaH = Elevation at pt 1 minus at pt 2
%pt 1 = [lon,lat] at pt 1
%pt 2 = [lon,lat] at pt 2
%slopeO = slope across

%If vector inputs given, ensure individual points are seperate rows
if length(deltaH(1,:)) > length(deltaH(:,1))
    deltaH = deltaH';
end
if length(pt1(:,1)) == 2 && length(pt1(1,:)) ~= 2
    pt1 = pt1';
end
if length(pt2(:,1)) == 2 && length(pt2(1,:)) ~= 2
    pt2 = pt2';
end

rEarth = 6371000; %Radius of Earth in meters.

%User Haversine Formula to calculate distance between two points along
%Earth's surface:
dHaver = 2*rEarth*asind(sqrt(sind(0.5*(pt1(2)-pt2(2)))^2 + cosd(pt1(2))*cosd(pt2(2))*(sind(0.5*(pt1(1)-pt2(1))).^2)));

slope = deltaH ./ dHaver;

% slope = deltaH ./ distance(pt1(:,2),pt1(:,1),pt2(:,2),pt2(:,1),referenceSphere('earth','meters'));