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

function mb_main2ice_grid(sHydro, sMeta)
%Translates changes from main grid to ice grid, which can be different.
%Ensures conservation of mass between the grids.

global sCryo

%Define variables used further down:
varLon = 'longitude';
varLat = 'latitude';
varIgrdLon = 'igrdlon';
varIgrdLat = 'igrdlat';


if ~isfield(sCryo, 'icmb')
    error('iceLink:noMb','The cryosphere structure array has no mass balance field (at main hydrologic grid scale). This is required for computing glacier dynamics.')
end

%Translate change in ice water equivalent modelled in using main grid to ice grid:
if isequal(sHydro.(varLat), sCryo.(varIgrdLat)) && isequal(sHydro.(varLon), sCryo.(varIgrdLon)) %Ice and main grids are same
    %mass balance of igrd same as mb of main grid
    sCryo.igrdmb = sCryo.icmb;
    %Update water equivalent of igrd
    sCryo.igrdwe  = sCryo.igrdwe + sCryo.igrdmb;
else %Ice and main grids are different
    sCryo.igrdmb = zeros(size(sCryo.igrddem), 'single');
    
    %Loop over all grid cells in the main grid to distribute changes
    indIc = find(sCryo.icmb ~= 0);
    for ii = 1 : numel(indIc)
        indIgrd = find(sCryo.igrdwe(sCryo.main2igrd{indIc(ii)}) > 0);
                
        %Distribute ice change from main grid to ice grid
        %Everything is units of depth, so area differences do not matter here
        %Mass balance:
        sCryo.igrdmb(indIgrd) = sCryo.icmb(indIc(ii)); %This grid is depth 
        %Water equivalent of igrd:
        sCryo.igrdwe(indIgrd) = full(sCryo.igrdwe(indIgrd)) + sCryo.igrdmb(indIgrd);
    end
    clear ii
end


%%Ensure ice water equivalent cannot have less than zero depth of ice:
sCryo.igrdwe(sCryo.igrdwe < 0) = 0;


%%Update fractional coverage of glaciers:
%re-initialize grid:
sCryo.icx = zeros(size(sCryo.icx), 'single');

%Test if same number of igrd cells in each main grid cell:
varGridSame = 'nigrdsame';
if ~isfield(sCryo, varGridSame)
    testIgrd = countmember(single((1:numel(sHydro.dem))), sCryo.igrd2main);
    
    if all(testIgrd(:) == testIgrd(1))
        sCryo.(varGridSame) = 1;
    else
        sCryo.(varGridSame) = 0;
    end
end

%Find indices of main grid with glacier:
indMnX = sCryo.igrd2main(sCryo.igrdwe > 0);
%Unique main grid cells with glaciers
[indMnXUnq] = unique(indMnX);
%Number of glacier grid cells contained in each main grid cell with
%glaciers:
mnXfreq = countmember(indMnXUnq,indMnX);
    
%Calculate fraction:
if sCryo.(varGridSame) == 1
    %Assign to ice fraction in main grid:
    sCryo.icx(indMnXUnq) = mnXfreq/100;
else
    %Check number of glacier grid cells in each main grid cell:
    nIgrdInMain = countmember(single((1:numel(sCryo.icx))'),sCryo.igrd2main);
    %Assign to ice fraction in main grid:
    sCryo.icx(indMnXUnq) = mnXfreq./nIgrdInMain(indMnXUnq);
end


%%Update icegrid surface elevation:
indUpdate = find(sCryo.igrdmb ~= 0);
if ~isempty(indUpdate)
    rhoI = find_att(sMeta.global,'density_ice'); 
    rhoW = find_att(sMeta.global,'density_water');

    sCryo.igrddem(indUpdate) = sCryo.igrdbsdem(indUpdate) + (rhoW/rhoI)*full(sCryo.igrdwe(indUpdate));
end
    
