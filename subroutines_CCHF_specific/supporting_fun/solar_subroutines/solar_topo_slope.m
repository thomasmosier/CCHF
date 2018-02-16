function [latSlope, lonDiff] = solar_topo_slope(lat, aspect, slopeAng, type)
%Find "equivalent" corrected latitude and longitude for sloped surfaces

%Slope angle should be as degrees inclination. Slope angle calculated by
%function 'watershed' is descent. Therefore, take absolute value
slopeAng = abs(slopeAng);

%Correct latitude to account for topography:
if regexpbl(type, 'flat')
    latSlope = lat;
    lonDiff = zeros(size(aspect));
elseif regexpbl(type, 'slope')
    %Ensure lat is grid rather than vector
    if  ~isequal(size(lat), size(aspect)) 
        if any(size(lat) == 1) && ~any(size(aspect) == 1)
            [~, lat] = meshgrid((1:numel(aspect(1,:))), lat);
        else
            error('solorTopoCorrect:latVector','The lat input is a vector, but grid is needed.')
        end
    end
    
    indFlat = find(aspect == -1); %-1 is value assigned to flat surfaces.
    %Correct latitude grid to account for slope and aspect:
        %Use negative of slope angle because equation calls for angle of
        %incline (see example on pg. 399) 
    latSlope = asind(sind(slopeAng).*cosd(aspect).*cosd(lat) + cosd(slopeAng).*sind(lat));
    indLatP = find(isnan(latSlope));
    %Set cells that are nan and cells that are flat to the original latitude
    latSlope(indLatP) = lat(indLatP);
    latSlope(indFlat) = lat(indFlat);
    
    
    %Difference in longitude between sloped surface and equivalent horizontal
    lonDiff = atand(...
        (sind(aspect).*sind(slopeAng)) ...
        ./ (cosd(slopeAng).*cosd(lat) - cosd(aspect).*sind(slopeAng).*sind(lat)) ...
        );
    lonDiff(isnan(lonDiff)) = 0;
else
    error('raadSlopeCorr:unknownTopoOpt', [Type ' is an unknown topographic option.']);
end

% %For testing:
% %Constant slope angle:
% lat = 30*ones([10, 1]);
% aspect = linspace(0, 360, 10)';
% slopeAng = ones([10, 1]);
% 
% %Constant aspect:
% lat = 55*ones([10, 1]);
% aspect = 135*ones([10, 1]);
% slopeAng = linspace(0, 30, 10)';
% 
% [lat, latSlope, lat - latSlope, lonDiff]