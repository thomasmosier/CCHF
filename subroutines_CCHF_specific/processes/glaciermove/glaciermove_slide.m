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
varIceLon = 'igrdlon';
varIceLat = 'igrdlat';
varDl = 'igrdfdrdl';
varIceFdr = 'igrdfdr';
varIceDem = 'igrddem';
varIceFdrM = 'igrdfdrmult';
varIceWgtF = 'igrd2mainFlowWgt';
varVel = 'igrdvel';
varIceWE = 'igrdwe';
varIceMb = 'igrdmb';
varMnMb = 'icmb';
varMnLat = 'iclateral';
varMnIceWe = 'icwe';
varMnIcx = 'icx';

szIce = size(sCryo.(varIceDem));
dt = time2sec(1,sMeta.dt,sMeta.dateCurr); 


%Set velocity minimum threshold (if less than this, don't bother with calculation):
velThresh = 0.01/(3.154*10^7); %(1 cm per year; in units of m/s)

%Calculate distance between grid cells along flow path:
if ~isfield(sCryo, varDl)
    rEarth = find_att(sMeta.global,'radius_Earth');
    
    dl = haversine_neighbors(sCryo.(varIceLon), sCryo.(varIceLat), rEarth,'sparse');

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
    sCryo.(varDl) = reshape(sum(sCryo.(varIceFdr).*dl, 2), size(sCryo.(varIceDem)));
%     sCryo.(varDl)(sCryo.(varDl) == 0) = nan; 
end

%Calculate multiple flow directions 
%Note: values in output are fractions that sum to 1. The values are based on magnitude of descent. 
if ~isfield(sCryo, varIceFdrM)
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
                    wflowacc(sCryo.(varIceLon), sCryo.(varIceLat), sCryo.(varIceDem), ...
                    'type','multi', 'edges','open', 'coord','geographic', 'routeflats', 'yes'); 
        %Find indices where mass balance crosses main grid cells
        [indIgrdF, indIgrdT] = find(iceFdrM ~= 0);
        %Find weighting:
        [fracMnFTemp, indMnF, indMnT] = ice2main_wgt(sCryo.(varIceLon), sCryo.(varIceLat), ...
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
    
    sCryo.(varIceFdrM) = iceFdrM; %Multiple flow path of glaciers at igrd (Number of ice grid cells by number of ice grid cells)
    sCryo.(varIceWgtF) = fracMnF; %Fractional weighting of flow at main (Number of main grid cells by number of main grid cells)
end


%Find grid cells where velocity exceeds threshold
indVel = find(sCryo.(varVel) > velThresh);
%Calculate lateral transport at both glacier and main grids if any moving
%ice
if ~isempty(indVel)
    %%CALCULATE LATERAL MOVEMENT AT ICE GRID
    %Calculate displacement as fraction of grid cell flow distance
    fracMove = dt*(sCryo.(varVel)(indVel)./sCryo.(varDl)(indVel)); %fractional displacement (unitless)
    %Set nan values to zero
    fracMove(isnan(fracMove)) = 0;
    fracMove(isinf(fracMove)) = 0;
    
    %Calculate displaced ice:
    dispIce = sparse(szIce(1), szIce(2));
    dispIce(indVel) = transpose(fracMove)*sCryo.(varIceWE)(indVel); 
    %Add ice to downslope grid cells:
    sCryo.(varIceWE)(:) = sCryo.(varIceWE)(:) + transpose(sCryo.(varIceFdrM))*reshape(dispIce, [], 1); %Units of meters
    %Remove ice from uphill grid cell:
    sCryo.(varIceWE)(:) = sCryo.(varIceWE)(:) - dispIce(:); %Units of meters

    %Add displaced ice to ice mass balance field (glacier grid):
    sCryo.(varIceMb)(:) = sCryo.(varIceMb)(:) + full(transpose(sCryo.(varIceFdrM))*reshape(dispIce, [],1) - dispIce(:));
    
    %%TRANSLATE MASS BALANCE FROM ICE GRID TO MAIN GRID
    %Initialize movement grid
    sCryo.(varMnLat) = zeros(size(sHydro.dem), 'single');
    %Calculate movement at main grid:
    sCryo.(varMnLat)(:) = ...
        transpose(sCryo.(varIceWgtF))*double(sCryo.(varMnIceWe)(:)) ...
        - double(full(sum(sCryo.(varIceWgtF), 2))).*sCryo.(varMnIceWe)(:);
%     %For testing:
%     figure; imagesc(sCryo.(varMnLat)); colorbar;
%     figure; imagesc(sCryo.(varMnIcx).*sCryo.(varMnIceWe)); colorbar;
%     figure; histogram(sCryo.(varMnLat));
    
    %Add lateral ice movement to ice array:
    sCryo.(varMnIceWe) = sCryo.(varMnIceWe) + sCryo.(varMnLat);
    %Set any negative values to zero:
    sCryo.(varMnIceWe)(sCryo.(varMnIceWe) < 0) = 0;
end