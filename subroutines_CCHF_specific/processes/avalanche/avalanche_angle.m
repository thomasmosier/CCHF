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


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'angle_critical', 30, 70, 50, 'avalanche','cryo'});
    
    return
else
    sMeta = varargin{1};
    angCrit = find_att(varargin{1}.coef, 'angle_critical');  
    angCritDisp = round2(angCrit, 3);
    %Stored as degrees for readability. Convert to rise over run:
    angCrit = tand(angCrit);
end

%Define variables used here:
varLon = 'longitude';
varLat = 'latitude';
varIgrdLon = 'igrdlon';
varIgrdLat = 'igrdlat';
varSnavFdr = 'snavfdr';
varSnav = 'snav';



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
    angCritDispFlr = floor(angCritDisp);
    dirStore = fullfile(sMeta.foldstorage, 'avalanche_angle', [num2str(angCritDispFlr) 'to' num2str(angCritDispFlr+1)]);
    if ~exist(dirStore, 'dir')
       mkdir(dirStore); 
    end
    pathStore = fullfile(dirStore, ...
        [sMeta.region{sMeta.siteCurr}, '_avalanche_', sMeta.iceGrid '_' ...
        num2str(szMn(1)) 'x' num2str(szMn(2)) '_' num2str(szIce(1)) 'x' ...
        num2str(szIce(2))  '_critical_angle_' num2str(angCritDisp) '.mat']);

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
        indAval = find(sCryo.igrdslopefdr(indIgrdF) <= -angCrit | sCryo.igrdslopefdr(indIgrdF) >= angCrit); 
        indIgrdF = indIgrdF(indAval);
        indIgrdT = indIgrdT(indAval);

        [fracMnF, indMnF, indMnT] = ice2main_wgt(sCryo.(varIgrdLon), sCryo.(varIgrdLat), ...
            sHydro.(varLon), sHydro.(varLat), sCryo.igrd2main, indIgrdF, indIgrdT);
        
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
    sCryo.(varSnavFdr) = sparse(double(indMnF),double(indMnT),double(fracMnF), prod(szMn), prod(szMn));
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
sCryo.(varSnav) = zeros(szMn, 'single');
%First term is snow gains, second term is snow losses 
sCryo.(varSnav)(:) = transpose(sCryo.(varSnavFdr))*double(reshape(sCryo.snw,[],1)) - double(full(sum(sCryo.(varSnavFdr), 2))).*sCryo.snw(:);

% %For testing:
% figure; imagesc(sCryo.(varSnav)); colorbar;
% figure; histogram(sCryo.(varSnav));

%Add snow avalanching to snow array:
sCryo.snw = sCryo.snw + sCryo.(varSnav);
%Set any negative values to zero:
sCryo.snw(sCryo.snw < 0) = 0;