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

function varargout = avalanche_angle(sHydro, varargin)
%Model avalanching of snow when slope is greater than critical angle
%This uses the glacier grid (instead of the main grid)


global sCryo

varLon = 'longitude';
varLat = 'latitude';
varIgrdLon = 'igrdlon';
varIgrdLat = 'igrdlat';

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'angle_critical', 1, 30, 10, 'avalanche','cryo'});
    
    return
else
    sMeta = varargin{1};
    angCrit = find_att(varargin{1}.coef, 'angle_critical');  
    
    %Stored as degrees for readability. Convert to rise over run:
    angCrit = tand(angCrit);
end



szMn = size(sHydro.dem);
if isfield(sCryo, 'igrddem')
    szIce = size(sCryo.igrddem);
else
    error('avalancheAngle:noIceDem', ['There is no ice DEM. This is ' ...
        'needed for the avalanche model to work. Set Avalanche to none ' ...
        'or load ice DEM.']);
end

%Create or load avalnche flow direction array
if ~isfield(sCryo,'snavfdr')
    %Create storage path:
    dirStore = sMeta.foldstorage;
    if ~exist(dirStore, 'dir')
       mkdir(dirStore); 
    end
    pathStore = fullfile(dirStore, ...
        [sMeta.region{sMeta.siteCurr}, '_avalanche_', sMeta.iceGrid '_' ...
        num2str(szMn(1)) 'x' num2str(szMn(2)) '_' num2str(szIce(1)) 'x' num2str(szIce(2))  '.mat']);
    
    if ~exist(pathStore, 'file') %File does not exist, so create
        %display message that this will take a long time
        if ~regexpbl(sMeta.mode, {'calib','valid'})
            tic
            if numel(sCryo.igrdslope(:)) > 10000
                disp('Calculating avalance grids. This may take some time. Grids will be stored for subsequent model runs.');
            end
        end
    
        %In flow direction grid, row is grid cell that flow is originating from
        %and column is cell the flow is going to
        [indIgrdF, indIgrdT] = find(sCryo.igrdfdr ~= 0);  
%         indFdr = sub2ind(size(sCryo.igrdslopefdr), indIgrdF, indIgrdT);
        indAval = find(sCryo.igrdslopefdr(indIgrdF) <= -angCrit | sCryo.igrdslope(indIgrdF) >= angCrit); 
        indIgrdF = indIgrdF(indAval);
        indIgrdT = indIgrdT(indAval);

        %Avalanching only matters when cross "main" hydrologic grid 
        indCross = find(sCryo.igrd2main(indIgrdF) ~= sCryo.igrd2main(indIgrdT));
        indIgrdF = indIgrdF(indCross);
        indIgrdT = indIgrdT(indCross);

        %Find unique pairs of main grid indices (multiple ice grid cells may flow to same set of main grid cells):
        %NOTE: the same main grid cell may flow to multiple main grid cells
        %due to differences in flow direction at fine grid scale.
        [temp, ~, ~] = unique([full(sCryo.igrd2main(indIgrdF)), full(sCryo.igrd2main(indIgrdT))],'rows');
        indMnF = temp(:,1);
        indMnT = temp(:,2);

        
        %%Calculate area-weighting factor for each main grid cell where avalanching crosses border. 
        %These factors will be used during each time step to redistribute
        %snow
        
        %Need areas for the cells that snow is avalanching from
        indMnArea  = unique(indMnF);
        indIceArea = unique(indIgrdF);
        
        %Initialize area arrays (same size as whole domain)
        areaMain = sparse(szMn(1), szMn(2));
        areaIce  = sparse(szIce(1), szIce(2));
        
        %Calculate areas
        %It could be the case this this function is not correct
        areaMain(indMnArea) = area_grid_pt( indMnArea,    sHydro.(varLon),    sHydro.(varLat));
        areaIce(indIceArea) = area_grid_pt(indIceArea, sCryo.(varIgrdLon), sCryo.(varIgrdLat));

        %%Calculate fractional contributions for cells that snow is avalanching from:
        %Note: fractional contributions come from sum of areas in ice grid
        %cells that avalance
        fracMnF = zeros(numel(indMnF), 1, 'single');
        
        indIgrd2MainAvF = sCryo.igrd2main(indIgrdF);
        indIgrd2MainAvT = sCryo.igrd2main(indIgrdT);
        
        %Loop over all main grid cells that are being looped over
        for ii = 1 : numel(indMnF)
            %Find ice grid indices contributing to current grid cell and
            %calculate fraction
            fracMnF(ii) = sum(full(areaIce(indIgrdF(indIgrd2MainAvF == indMnF(ii) & indIgrd2MainAvT == indMnT(ii)))))...
                /areaMain(indMnF(ii));
        end
        
        %Display time that has been taken
        if ~regexpbl(sMeta.mode, {'calib','valid'})
            tAv = round2(toc/60,2);
            if tAv > 0.5
                disp(['It took ' num2str(tAv) ' minutes to calculate the avalanche grids.']);
            end
        end
    
        save(pathStore, 'indMnF', 'indMnT', 'fracMnF');
    else %File exists, so load
       load(pathStore) 
    end

    
    %Create snow avalance "flow direction (and fraction)" grid
    sCryo.snavfdr = sparse(double(indMnF),double(indMnT),double(fracMnF), prod(szMn), prod(szMn));
end


% %Test resulting grids:
% unqCellsFrom = unique(indMnF);
% 
% for ii = 1 :numel(unqCellsFrom)
%     fracCurr = sum(fracMnF(unqCellsFrom(ii) == indMnF));
%     if fracCurr > 1
%         disp([num2str(unqCellsFrom(ii)) ' = ' num2str(fracCurr)])
%     end
% end


%%REDISTRIBUTE AVALANCHED SNOW AT MAIN GRID RESOLUTION:
%Reminder of nomenclature
% strIndMnF = indices in main grid that are losing snow
% strIndMnT = indices in main grid that are gaining snow
% strFrcMnF = fraction of snow in main grid that is being lost

%Calculate snow avalanching:
sCryo.snav = zeros(szMn, 'single');
%First term is snow gains, second term is snow losses 
sCryo.snav(:) = transpose(sCryo.snavfdr)*double(reshape(sCryo.snw,[],1)) - double(full(sum(sCryo.snavfdr, 2))).*sCryo.snw(:);

% %For testing:
% figure; imagesc(sCryo.snav); colorbar;
% figure; histogram(sCryo.snav);

%Add snow avalanching to snow array:
sCryo.snw = sCryo.snw + sCryo.snav;
%Set any negative values to zero:
sCryo.snw(sCryo.snw < 0) = 0;