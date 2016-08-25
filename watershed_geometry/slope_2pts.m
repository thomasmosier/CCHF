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