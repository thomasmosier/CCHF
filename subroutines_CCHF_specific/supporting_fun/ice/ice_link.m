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

function ice_link(sHydro)
%Translates changes from main grid to ice grid, which can be different.
%Ensures conservation of mass between the grids.

global sCryo

%Translate change in ice water equivalent modelled in mass and energy 
%function (using main grid) to ice grid:
if isequal(sHydro.lat, sCryo.iceLat) && isequal(sHydro.lon, sCryo.iceLon) %Ice and main grids are same
    sCryo.icdwe = sCryo.icdwe;
    sCryo.icwe = sCryo.icwe + sparse(double(sCryo.icdwe));
else %Ice and main grids are different
    sCryo.icgrddwe = zeros(size(sCryo.icgrddem));
    
%     %Calculate areas of ice grid cells
%     if ~isfield(sCryo,'iceArea')
%         sCryo.iceArea = area_geodata(sCryo.iceLon, sCryo.iceLat);
%     end
    
    warning('ice_link:needToEdit','This needs to be edited');
    %Loop over all grid cells in the main grid to distribute changes
    for ii = 1 : numel(sHydro.dem)
        indIce = find(sCryo.icgrdwe(sCryo.main2ice{ii}) > 0);
                
%         %Calculate area ratio between sHydro.area.*sCryo.icx and
%         %sum(sCryo.iceArea)
%         %Use this to distribute water equivalent change from main grid to
%         %fine grid.
%         %**This assumes that all ice water change (positive and negative)
%         %goes to cells that currently have ice.
%         areaRatio = sHydro.area(ii).*sCryo.icx(ii) ...
%             / sum(sCryo.iceArea(sCryo.main2ice{ii}(indIce)));
        
        %Distribute ice change from main grid to ice grid:
        sCryo.icgrddwe(indIce) = sCryo.icgrddwe(ii); %This grid is depth 
        %equivalent (regardless of fractional glacier coverage)
        sCryo.icgrdwe(indIce) = sCryo.icgrdwe(indIce) + sCryo.icgrddwe(indIce);
    end
end
    
    
