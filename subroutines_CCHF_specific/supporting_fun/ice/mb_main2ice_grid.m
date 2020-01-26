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

%Main hydrologic grid (also has some ice areas)
varMnLatMov = 'iclateral';
varMnIceWe = 'icwe';
varMnIceDwe = 'icdwe'; %Change in water equivalent (main grid)
varMnIcx = 'icx';
varMnMb = 'icmb'; %Glacier mass balance (annual basis)

%Ice grid (possibly finer than main hydrologic grid)
varIgrdArea = 'igrdarea';
varIgrdDl = 'igrdfdrdl';
varIgrdFdr = 'igrdfdr';
varIgrdDem = 'igrddem';
varIgrdFdrM = 'igrdfdrmult';
varIgrdWgtF = 'igrd2mainFlowWgt';

varIgrdWe = 'igrdwe';
varIgrdMb = 'igrdmb'; %Change in water equivalent (ice grid)


if ~isfield(sCryo, 'icmb')
    error('iceLink:noMb','The cryosphere structure array has no mass balance field (at main hydrologic grid scale). This is required for computing glacier dynamics.')
end


if isequal(sHydro.(varLat), sCryo.(varIgrdLat)) && isequal(sHydro.(varLon), sCryo.(varIgrdLon)) 
    blGrdResSame = 1;
else
    blGrdResSame = 0;
end

%Translate change in ice water equivalent modelled in using main grid to ice grid:
if blGrdResSame == 1 %Ice and main grids are same
    %mass balance of igrd same as mb of main grid
    sCryo.(varIgrdMb) = sCryo.(varMnMb);
    %Update water equivalent of igrd
    sCryo.(varIgrdWe)  = full(sCryo.(varIgrdWe)) + sCryo.(varIgrdMb);
else %Ice and main grids are different
    sCryo.(varIgrdMb) = zeros(size(sCryo.igrddem), 'single');
    
    %Loop over all grid cells in the main grid to distribute changes
    %Note: There are several ice grid cells for each main grid cell
    indIc = find(sCryo.icmb ~= 0);
    for ii = 1 : numel(indIc)
        indIgrd = find(sCryo.(varIgrdWe)(sCryo.main2igrd{indIc(ii)}) > 0);
                
        %Distribute ice change from main grid to ice grid
        %Everything is units of depth, so area differences do not matter here
        %Mass balance:
        sCryo.(varIgrdMb)(indIgrd) = sCryo.icmb(indIc(ii)); %This grid is depth 
        %Water equivalent of igrd:
        sCryo.(varIgrdWe)(indIgrd) = full(sCryo.(varIgrdWe)(indIgrd)) + sCryo.(varIgrdMb)(indIgrd);
    end
    clear ii
end


%%Ensure ice water equivalent cannot have less than zero depth of ice:
sCryo.( varIgrdWe)(sCryo.( varIgrdWe) < 0) = 0;
sCryo.(varMnIceWe)(sCryo.(varMnIceWe) < 0) = 0;


%%Update fractional coverage of glaciers:
if blGrdResSame
    sCryo.(varMnIcx)(sCryo.(varMnIceWe) <= 0 | isnan(sCryo.(varMnIceWe))) = 0;
    sCryo.(varMnIcx)(sCryo.(varMnIcx) < 0) = 0;
else
    %re-initialize grid:
    sCryo.(varMnIcx) = zeros(size(sCryo.(varMnIcx)), 'single');

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
    indMnX = sCryo.igrd2main(sCryo.(varIgrdWe) > 0);
    %Unique main grid cells with glaciers
    [indMnXUnq] = unique(indMnX);
    %Number of glacier grid cells contained in each main grid cell with
    %glaciers:
    mnXfreq = countmember(indMnXUnq,indMnX);

    %Calculate fraction:
    if sCryo.(varGridSame) == 1
        %Assign to ice fraction in main grid:
        sCryo.(varMnIcx)(indMnXUnq) = mnXfreq/100;
    else
        %Check number of glacier grid cells in each main grid cell:
        nIgrdInMain = countmember(single((1:numel(sCryo.(varMnIcx)))'),sCryo.igrd2main);
        %Assign to ice fraction in main grid:
        sCryo.(varMnIcx)(indMnXUnq) = mnXfreq./nIgrdInMain(indMnXUnq);
    end
end


%Enforce physical constraints
sCryo.(varMnIcx)( sCryo.(varMnIcx)  < 0) = 0;


%%Update icegrid surface elevation:
indUpdate = find(sCryo.(varIgrdMb) ~= 0);
if ~isempty(indUpdate)
    rhoI = find_att(sMeta.global,'density_ice'); 
    rhoW = find_att(sMeta.global,'density_water');

    sCryo.igrddem(indUpdate) = sCryo.igrdbsdem(indUpdate) + (rhoW/rhoI)*full(sCryo.(varIgrdWe)(indUpdate));
end
    
