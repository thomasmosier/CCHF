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

function varargout = glaciermove_slide(sHydro, varargin)

global sCryo

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    return
else
    sMeta = varargin{1};
end


varLon = 'longitude';
varLat = 'latitude';
varIgrdLon = 'igrdlon';
varIgrdLat = 'igrdlat';

%Main hydrologic grid (also has some ice areas)
varMnLatMov = 'iclateral';
varMnIceWe = 'icwe';
varMnIceDwe = 'icdwe';
varMnIcx = 'icx';

%Ice grid (possibly finer than main hydrologic grid)
varIgrdArea = 'igrdarea';
varIgrdDl = 'igrdfdrdl';
varIgrdFdr = 'igrdfdr';
varIgrdDem = 'igrddem';
varIgrdFdrM = 'igrdfdrmult';
varIgrdWgtF = 'igrd2mainFlowWgt';
varIgrdVel = 'igrdvel';
varIgrdWe = 'igrdwe';
varIgrdMb = 'igrdmb';


szIce = size(sCryo.(varIgrdDem));
dt = time2sec(1,sMeta.dt,sMeta.dateCurr); 


%Set velocity minimum threshold (if less than this, don't bother with calculation):
velThresh = 0.01/(3.154*10^7); %(1 cm per year; in units of m/s)


%Check for velocity field:
if ~isfield(sCryo, varIgrdVel)
   error('glaciermoveSlide:noVel', ['Glacier movement can only be calculated if ' ...
       'glacier velocity is also calculated. This is done in the ' char(39) glaciervel char(39) ' process representation.']) 
end


%Calculate distance between grid cells along flow path:
if ~isfield(sCryo, varIgrdDl)
    rEarth = find_att(sMeta.global,'radius_Earth');
    
    dl = haversine_neighbors(sCryo.(varIgrdLon), sCryo.(varIgrdLat), rEarth,'sparse');

    %Create field that is distance along flowpath:
    %Notes: Ice flow direction used for estimating flow distance includes 
    %single flow path; For actually calculating displacement, multiple flow 
    %paths are used. This is therefore an approximation. Eeach row of both 
    %flow direction grids sums to 1
    %(1) Element-wise multiplication ensures dl is matched with each flow
    %path. 
    %(2) Sum over dim 2 retrieves dl corresponding to each originating grid
    %cell along the flow path
    %(3) Reshape formulates the dl vector into an array with shape of igrd
    sCryo.(varIgrdDl) = reshape(sum(sCryo.(varIgrdFdr).*dl, 2), size(sCryo.(varIgrdDem)));
%     sCryo.(varDl)(sCryo.(varDl) == 0) = nan; 
end

%Calculate multiple flow directions 
%Note: values in output are fractions that sum to 1. The values are based on magnitude of descent. 
if ~isfield(sCryo, varIgrdFdrM)
    %Create storage path:
    dirStore = sMeta.foldstorage;
    if ~exist(dirStore, 'dir')
       mkdir(dirStore); 
    end
    pathStore = fullfile(dirStore, ...
        [sMeta.region{sMeta.siteCurr}, '_glaciermove_slide_weight_', sMeta.iceGrid '_' ...
        num2str(numel(sHydro.(varLat))) 'x' num2str(numel(sHydro.(varLon))) '_' num2str(szIce(1)) 'x' num2str(szIce(2))  '.mat']);

    if ~exist(pathStore, 'file') %File does not exist, so create
        %display message that this will take a long time
        if ~regexpbl(sMeta.mode, 'calib')
            tic
            if numel(sCryo.igrdslope(:)) > 10000
                disp('Calculating gracier flow grids. This may take some time. Grids will be stored for subsequent model runs.');
            end
        end
        
        %Flow and weighting are calculated only on first iteration. This is
        %an approximation, because it may change if glacier changes
        %dramatically, but likely not over ~100 year periods since it 
        %should depend on basal elevation.
        [~, iceFdrM, ~] = ...
                    wflowacc(sCryo.(varIgrdLon), sCryo.(varIgrdLat), sCryo.(varIgrdDem), ...
                    'type','multi', 'edges','open', 'coord','geographic', 'routeflats', 'yes'); 
        %Find indices where mass balance crosses main grid cells
        [indIgrdF, indIgrdT] = find(iceFdrM ~= 0);
        %Find weighting:
        [fracMnFTemp, indMnF, indMnT] = ice2main_wgt(sCryo.(varIgrdLon), sCryo.(varIgrdLat), ...
            sHydro.(varLon), sHydro.(varLat), sCryo.igrd2main, indIgrdF, indIgrdT);
    
        %Display time that has been taken
        if ~regexpbl(sMeta.mode, 'calib')
            tAv = round2(toc/60,2);
            if tAv > 0.5
                disp(['It took ' num2str(tAv) ' minutes to calculate the avalanche grids.']);
            end
        end
        
        fracMnF = sparse(double(indMnF),double(indMnT),double(fracMnFTemp), numel(sHydro.dem), numel(sHydro.dem));

        save(pathStore, 'iceFdrM', 'fracMnF');
    else %File exists, so load
       load(pathStore) 
    end
    
    sCryo.(varIgrdFdrM) = iceFdrM; %Multiple flow path of glaciers at igrd (Number of ice grid cells by number of ice grid cells)
    sCryo.(varIgrdWgtF) = fracMnF; %Fractional weighting of flow at main (Number of main grid cells by number of main grid cells)
end


%Determine if ice grid is same res as main hydrologic grid 
if isequal(sHydro.(varLat), sCryo.(varIgrdLat)) && isequal(sHydro.(varLon), sCryo.(varIgrdLon)) 
    blGrdResSame = 1;
    %In this case, consider fractional glacier coverage within grid cell
else
    blGrdResSame = 0;
end


%Find grid cells where velocity exceeds threshold
indVel = find(sCryo.(varIgrdVel) > velThresh);
%Calculate lateral transport at both glacier and main grids if any moving
%ice
if ~isempty(indVel)
    %%CALCULATE LATERAL MOVEMENT AT ICE GRID
    %Calculate displacement as fraction of grid cell flow distance
    fracMove = dt*(sCryo.(varIgrdVel)(indVel)./sCryo.(varIgrdDl)(indVel)); %fractional displacement (unitless)
    %Set nan values to zero
    fracMove(isnan(fracMove)) = 0;
    fracMove(isinf(fracMove)) = 0;
    
    
    %Initialize displaced ice:
    dispIce = zeros(szIce(1), szIce(2));
    %Calculate displaced ice (displaced volume = frac move * elevation * area * fraction cover):
    if blGrdResSame == 1
        %Delta Ice = fractional displacement * Area grid cell * fractional
        %glacier coverage * Depth of icea
        dispIce(indVel) = double(full(transpose(fracMove))*(sCryo.(varIgrdWe)(indVel).*sHydro.area(indVel).*sCryo.(varMnIcx)(indVel)));
    else
        %Delta Ice volume = fractional displacement * Area grid cell * Depth of icea
        
        %Calculate grid cell area:
        if ~isfield(sCryo, varIgrdArea)
            sCryo.(varIgrdArea) = area_geodata(sCryo.(varIgrdLon), sCryo.(varIgrdLat));
        end
        
        %Calculate delta ice volume:
        dispIce(indVel) = double(full(transpose(fracMove))*(sCryo.(varIgrdWe)(indVel).*sCryo.(varIgrdArea)(indVel)));  
    end
    
    %Calculate ice that is added to downslope grid cells (assumes same outline; additional ice adds to depth)
    if blGrdResSame == 1
        %New ice depth = (delta ice / (area grid cell * fractional coverage))
        matNewIce = reshape(dispIce, [], 1) ./ (sHydro.area(:).*sCryo.(varMnIcx)(:));
        %Use flow direction to move downslope:
        matNewIceDownSlp = (full(transpose(sCryo.(varIgrdFdrM))*sparse(matNewIce)));
    else
        %New ice depth = (delta ice / area grid cell)
        matNewIce = reshape(dispIce, [], 1) ./ sCryo.(varIgrdArea)(:);
        %Use flow direction to move downslope:
        matNewIceDownSlp = (full(transpose(sCryo.(varIgrdFdrM))*sparse(matNewIce)));
    end
    
    %Add displaced ice to ice water equivalent
    sCryo.(varIgrdWe)(:) = sCryo.(varIgrdWe)(:) + matNewIceDownSlp; %Units of meters cubed

    %Add ice to delta icea water equivalent
    sCryo.(varIgrdMb)(:) = sCryo.(varIgrdMb)(:) + matNewIceDownSlp;

    %Remove displaced ice from uphill grid cell (assumes same outline; removes from depth)
    sCryo.(varIgrdWe)(:) = sCryo.(varIgrdWe)(:) - matNewIce; %Units of meters cubed

    %Ensure no negative ice depths:
    sCryo.(varIgrdWe)(sCryo.(varIgrdWe) < 0) = 0;
    
    
    %%TRANSLATE MASS BALANCE FROM ICE GRID TO MAIN GRID
    %Initialize movement grid
    sCryo.(varMnLatMov) = zeros(size(sHydro.dem), 'single');
    if blGrdResSame == 1
        %Ice moved downhill minus ice removed from uphill
        sCryo.(varMnLatMov)(:) = matNewIceDownSlp - matNewIce;
    else
        warning('glaciermoveSlide:translationIgrdToMn', ...
            ['Check the translation between ice movement at the ice '...
            'grid to main grid. For exmaple, ensure conservation of mass.']);
        sCryo.(varMnLatMov)(:) = ...
            transpose(sCryo.(varIgrdWgtF))*double(sCryo.(matNewIce)(:)) ...
            - double(full(sum(sCryo.(varIgrdWgtF), 2))).*sCryo.(matNewIce)(:);
    end
    
    %Add lateral ice movement to ice arrays (main grid):
    sCryo.(varMnIceWe)  = sCryo.(varMnIceWe) + sCryo.(varMnLatMov);
    sCryo.(varMnIceDwe) = sCryo.(varMnIceDwe) + sCryo.(varMnLatMov);
    
    %Set any negative values to zero:
    sCryo.(varMnIceWe)(sCryo.(varMnIceWe) < 0) = 0;
    
%     %For testing:
%     figure; imagesc(sCryo.(varMnLatMov)); colorbar;
%     figure; imagesc(sCryo.(varMnIcx).*sCryo.(varMnIceWe)); colorbar;
%     figure; histogram(sCryo.(varMnLatMov));
end