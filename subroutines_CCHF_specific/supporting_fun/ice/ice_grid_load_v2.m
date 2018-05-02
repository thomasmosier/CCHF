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

function [sIceInit] = ice_grid_load_v2(sHydro, sPath, sMeta)

%Ice field naming convention
%".ic" is at same resolution at main hydrologic grid.
%"igrd" is either at same resolution at main hydrologic grid or finer
%resolution. It may include some sparse arrays
szMain = size(sHydro.dem); 
% icBlMain = zeros(szMain, 'single');

varLon = 'longitude';
varLat = 'latitude';

%Load glacier coverage data and convert to depth:
%If glaciers or boolean clipping being used, create grids:
if isfield(sPath, 'ice') %&& ~regexpbl(sMeta.mode, 'parameter')

    %Create path to save ice structure to:
    [foldIceOutline, fileIceOutline, ~] = fileparts(sPath.ice);
%     [~, fileDem, ~] = fileparts(sPath.dem);
    if regexpbl(sMeta.iceGrid, 'same')
        pathIceStruct = fullfile(foldIceOutline, ...
            [fileIceOutline, '_', sMeta.region{sMeta.siteCurr}, '_', sMeta.iceGrid '.mat']);
    elseif regexpbl(sMeta.iceGrid, 'fine')
        if isfield(sPath, 'iceDem')
            [~, fileIceDem, ~] = fileparts(sPath.iceDem);
            
            pathIceStruct = fullfile(foldIceOutline, ...
                [fileIceOutline, '_', sMeta.region{sMeta.siteCurr}, '_', sMeta.iceGrid '_' fileIceDem '.mat']);
        else
            error('iceGridLoad:noIceDEM',['A fine resolution ice grid '...
                'has been selected, but no corresponding fine resolution ice DEM is present.']);
        end
    end

    
    if exist(pathIceStruct, 'file')
       load(pathIceStruct);
    else
        %Ice grid can either be 'same' as the main model grid or can be
        %'fine'. The advtange of 'fine' is that it better captures
        %local glacier dynamics
        
        disp(['Loading and processing ice grid information. '...
            'This may take quite some time (a few minutes for grids '...
            'with 200 indices, but a few hours for very large domains) ' ...
            'depending on the size of the domain and the options chosen.' ...
            char(10) 'Output will automatically be saved to skip this '...
            'process in subsequent runs for this model setup.']);

        edgeLatMain = box_edg(sHydro.(varLat));
        edgeLonMain = box_edg(sHydro.(varLon));


        %Use edg to crop ice dem:
        if regexpbl(sMeta.iceGrid, 'fine')
            %Crop fine grid (in case larger than main)
            lonCrp = [min(edgeLonMain), max(edgeLonMain)];
            latCrp = [min(edgeLatMain), max(edgeLatMain)];
            
            sIceDem = read_geodata_v2(sPath.iceDem, 'dem', lonCrp, latCrp, nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
                sIceDem.dem = squeeze(sIceDem.dem);
            
            lonIceIn = sIceDem.(varLon);
            latIceIn = sIceDem.(varLat);
            
            %Find scale factor between two grids:
            nRefin = round(nanmean([abs(nanmean(diff(sHydro.(varLon))))/abs(nanmean(diff(sIceDem.(varLon)))), ...
                abs(nanmean(diff(sHydro.(varLat))))/abs(nanmean(diff(sIceDem.(varLat))))]));
            
            %Create output grid such that edges of fine resolution are
            %consistent with edges of main grid
            [sIceDem.(varLon), sIceDem.(varLat)] = grid_refine(edgeLonMain, edgeLatMain, nRefin, 'edge');

            %Interpolate fine resolution ice grid to ensure complete
            %alignment with main hydrologic grid
            [lonInMesh,   latInMesh] = meshgrid(lonIceIn, latIceIn);
            [lonOutMesh, latOutMesh] = meshgrid(sIceDem.(varLon), sIceDem.(varLat));

            sIceDem.dem = interp2(lonInMesh, latInMesh, sIceDem.dem, ...
                lonOutMesh, latOutMesh, 'linear');
        elseif regexpbl(sMeta.iceGrid, 'same')
            sIceDem.dem = squeeze(sHydro.dem);
            sIceDem.(varLon) = sHydro.(varLon);
            sIceDem.(varLat) = sHydro.(varLat);
        else
            error('CCHF_vackbone:unknownIceGridType', ['The ice '...
                'grid method is ' iceGrid ', which is not known.']);
        end


        %Load ice extents (either shape or dem) and get indices containing ice:
        if regexpbl(sPath.ice,'.shp')
%             disp('Reading ice shapefile.');
            shpIce = shaperead(sPath.ice);

            if isfield(shpIce,'Area')
                strIceLd = 'Area';
            else
                warning('iceGridLd:noArea','The input shapefile does not have a glacier area field. Therefore, this attribute will not be recorded.');
                if isfield(shpIce,'Slope')
                    strIceLd = 'Slope';
                elseif isfield(shpIce,'Aspect')
                    strIceLd = 'Aspect';
                else
                    error('CCHF_backbone:unknownShpField',['A suitable '...
                        'field for querrying the ice shapefile was not found.'])
                end
            end
                warning('off','interpshapefile:overlap');
                %Use original interpshapefile, which is faster than
                %"interpshapefile_geo" (used below)
%             disp('Interpolating shapefile to grid.');
            [iceExist, iceArea] = shapefile_geo_v2(shpIce, strIceLd, sIceDem.(varLat), sIceDem.(varLon));
%             iceExist = interpshapefile_geo_v2(shpIce, sIceDem.(varLat), sIceDem.(varLon), strIceLd, 1);
            indIce = find(~isnan(iceExist));
%             disp('finished first ice shp interp')
        else
            sIceInput = read_geodata_v2(sPath.ice, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
                sIceInput.data = squeeze(sIceInput.data);
            %Check that coordinates are same as iceDEM:
            lonRes = log10(round(1/nanmean(abs(diff(sIceInput.(varLon))))))+1;
            latRes = log10(round(1/nanmean(abs(diff(sIceInput.(varLat))))))+1;
            if numel(sIceInput.(varLon)) ~= numel(sIceDem.(varLon)) || numel(sIceInput.(varLat)) ~= numel(sIceDem.(varLat)) 
                error('CCHF_backbone:diffIceres','The ice DEM and ice presence grid have different sizes but the option has been state such that the grids are expected to have the same size.');
            elseif sum(ismember(round2(sIceInput.(varLon), lonRes), round2(sIceDem.(varLon), lonRes))) < 0.9*numel(sIceDem.(varLon)) || sum(ismember(round2(sIceInput.(varLat), latRes), round2(sIceDem.(varLat), latRes))) < 0.9*numel(sIceDem.(varLat))
               error('CCHF_backbone:misalignedIceDem','The ice DEM and ice presence grid do not appear to align. Maybe rounding error?');
            end

            %Grids basicaly align, so ensure they match completely.
            [indLonIcePre, indLonIceDem] = ismember(round2(sIceInput.(varLon), lonRes), round2(sIceDem.(varLon), lonRes));
            [indLatIcePre, indLatIceDem] = ismember(round2(sIceInput.(varLat), latRes), round2(sIceDem.(varLat), latRes));

            if all(indLonIcePre == 1 | indLonIcePre == 0)
                indLonIcePre = find(indLonIcePre == 1);
                indLatIcePre = find(indLatIcePre == 1);
            end
            if all(indLonIceDem == 1 | indLonIceDem == 0)
                indLonIceDem = find(indLonIceDem == 1);
                indLatIceDem = find(indLatIceDem == 1);
            end

            sIceInput.data = sIceInput.data(indLatIcePre, indLonIcePre);
            sIceInput.(varLon)  = sIceInput.(varLon)(indLonIcePre);
            sIceInput.(varLat)  = sIceInput.(varLat)(indLatIcePre);

            sIceDem.dem = sIceDem.dem(indLatIceDem, indLonIceDem);
            sIceDem.(varLon)  = sIceDem.(varLon)(indLonIceDem);
            sIceDem.(varLat)  = sIceDem.(varLat)(indLatIceDem);

            indIce = find(~isnan(sIceInput.data) & sIceInput.data ~= 0);
        end

        %Populate ice presence grid (at main resolution):
        szIce = size(sIceDem.dem);
        
        %Initialize ice structure array output:
        sIceInit = struct;
        %Initialize water equivalent as sparse array:
        sIceInit.igrdwe = sparse(szIce(1), szIce(2));
        %Initalize to 1 where there is ice
        sIceInit.igrdwe(indIce) = 1;
        
%         %Record ice area (if loaded)
%         if regexpbl(sPath.ice,'.shp') && regexpbl(strIceLd, 'Area')
%             sIceInit.igrdarea = iceArea;
%         end


        %Estimate fractional glacier extent within each grid cell
        nScl = 5;
        
        [lonSub, latSub] = grid_refine(edgeLonMain, edgeLatMain, nScl, 'edge');

        %This can be taken out:
%         nScl = 5;
%         lonSub = nan(nScl*numel(sHydro.(varLon)),1);
%         lonSubInd = cell(numel(sHydro.(varLon)),1);
%         for ll = 1 : numel(sHydro.(varLon))
%             cntrSub = (ll-1)*(nScl) + 1;
%             lonSubInd{ll} = (cntrSub : cntrSub + nScl-1);
% 
%             lonSubCurr = linspace(edgeLon(ll), edgeLon(ll+1), nScl+1);
%             lonSubCurr = lonSubCurr(1:nScl) + 0.5*(mean(abs(diff(lonSubCurr(1:nScl)))));
%             lonSub(lonSubInd{ll}) = lonSubCurr;
%         end
%         latSub = nan(nScl*numel(sHydro.(varLat)),1);
%         latSubInd = cell(numel(sHydro.(varLat)),1);
%         for ll = 1 : numel(sHydro.(varLat))
%             cntrSub = (ll-1)*(nScl) + 1;
%             latSubInd{ll} = (cntrSub : cntrSub + nScl-1);
% 
%             latSubCurr = linspace(edgeLat(ll), edgeLat(ll+1), nScl+1);
%             latSubCurr = latSubCurr(1:nScl) - 0.5*(mean(abs(diff(latSubCurr(1:nScl)))));
%             latSub(latSubInd{ll}) = latSubCurr;
%         end

%         disp('Interpolating ice shapefile to finer mesh in order to estimate fractional ice area in main grid.');
        if regexpbl(sPath.ice,'.shp')
                warning('off','interpshapefile:overlap');
            %This part eats a lot of memory. Perform in chunks
            iceSubExist = shapefile_geo_bl(shpIce, latSub, lonSub);
            %iceSubExist = interpshapefile_geo_v2(shpIce, latSub, lonSub, strIceLd, 1);
%                 iceSubExist(~isnan(iceSubExist)) = 1;
%                 iceSubExist(isnan(iceSubExist)) = 0;
        else
            warning('CCHF_backbone:iceResample',['The input ice grid '...
                'is being resampled for the purpose of estimating ' ...
                'fractional ice coverage.' char(10) 'This uncertainty can be '...
                'avoided through using a shapefile instead.' char(10) 'It '...
                'should not be a problem if the glacier grid is '...
                'significantly finer resolution than the main grid.']);

            tempIceBl = zeros(size(sIceInit.igrdwe),'single');
            tempIceBl(full(sIceInit.igrdwe ~= 0) & full(~isnan(sIceInit.igrdwe))) = 1;
            [lonMesh,latMesh] = meshgrid(sIceDem.(varLon), sIceDem.(varLat));
            [lonSubMesh,latSubMesh] = meshgrid(lonSub, latSub);
            iceSubExist = interp2(lonMesh,latMesh, tempIceBl, lonSubMesh,latSubMesh,'nearest');
        end
    
        %
        %Checking edits here
        %
        %Estimate fractional ice cover (main grid resolution):
%         disp('Estimating fractional ice area in main grid.')
        nSclSq = nScl^2;
        
        sIceInit.icx = sparse(szMain(1), szMain(2));

        for ll = 1 : numel(sHydro.(varLat))
            indLatN = find(latSub <= edgeLatMain(ll), 1, 'first');
            indLatS = find(edgeLatMain(ll+1) <= latSub, 1, 'last');
    
            for mm = 1 : numel(sHydro.(varLon))
                indLonW = find(edgeLonMain(mm) <= lonSub, 1, 'first');
                indLonE = find(lonSub <= edgeLonMain(mm+1), 1, 'last');
    
                sIceInit.icx(ll,mm) = sum2d(iceSubExist(indLatN:indLatS,indLonW:indLonE))/nSclSq;
            end
        end
%         disp('finished fractional ice shp interp')

%         %Note all main grid cells that contain any ice
%         sIceInit.icbl = sparse(szMain(1), szMain(2));
%         sIceInit.icbl(sIceInit.icx > 0) = 1;

        %Create sparse ice grid (including frame) and delete previous structures:                
        %Find indices of frame around existing ice (allows glacier morphology)

%         disp('Calculating surface flow properties for ice grid.');
        
%         %Set # of indices in frame around existing glaciers to crop ice DEM (want approximately 1
%         %km)
%             rEarth = find_att(sMeta.global,'radius_Earth');
%         frIce = ceil(1000/haversine_reference(sIceDem.(varLon)(1),sIceDem.(varLat)(1),sIceDem.(varLon)(2),sIceDem.(varLat)(2), rEarth));
%             frIce = max(frIce,2);
%         %In for loop, set all indices within frame to 0. After for loop, find indices of all 0 elements.    
%         [rIce, cIce] = ind2sub(szIce, indIce);
%         indDimIce = reshape((1:(szIce(1)*szIce(2))), szIce);
% 
%         for ll = 1 : numel(indIce)
%             rIceCurr = (rIce(ll)-frIce:rIce(ll)+frIce);
%                 rIceCurr(rIceCurr < 1) = [];
%                 rIceCurr(rIceCurr > szIce(1)) = [];
%             cIceCurr = (cIce(ll)-frIce:cIce(ll)+frIce);
%                 cIceCurr(cIceCurr < 1) = [];
%                 cIceCurr(cIceCurr > szIce(2)) = [];
% 
%             indDimIce(rIceCurr, cIceCurr) = 0;
%         end
% 
%         indIceFr = find(indDimIce == 0);
% 
%         sIceInit.igrddem =  sparse(szIce(1), szIce(2));
%         sIceInit.igrddem(indIceFr) = sIceDem.dem(indIceFr);
%         sIceInit.igrdlat = sIceDem.(varLat);
%         sIceInit.igrdlon = sIceDem.(varLon);
%             clear('sIceInput','sGlaciers');
%             clear('sIceDem');
            
        sIceInit.igrddem = sIceDem.dem;
        sIceInit.igrdlat = sIceDem.(varLat);
        sIceInit.igrdlon = sIceDem.(varLon);
            clear('sIceInput','sGlaciers');
            clear('sIceDem');


        %Calculate ice surface slope and flow direction (based on surface elevation)
        [~, sIceInit.igrdfdr, iceSurfSlope] = ...
            wflowacc(sIceInit.igrdlon, sIceInit.igrdlat, sIceInit.igrddem,'type','single','edges','open','coord','geographic');

        %Calculate slope in the flow direction:
        sIceInit.igrdslopefdr = -sparse(reshape(sum(sIceInit.igrdfdr.*iceSurfSlope, 2), szIce));
%             sIceInit.igrdslopefdr(sIceInit.igrddem == 0) = 0;
        %Calculate slope in all directions:
        sIceInit.igrdslope = sparse(iceSurfSlope);
%             sIceInit.igrdslope(sIceInit.igrddem == 0, :) = 0;
%             sIceInit.igrdslope(:, sIceInit.igrddem == 0) = 0;


        %%Find correspondance between glacier grid and main grid
        %Check which of main grid cell the centroid of the glacier
        %grid cell is contained within

        %Initialize two gris:
        %sparse grid = tracks which main grid cells contain 
            %glaciers 
        %cell array = Dim. of main grid; each non-empty cell-array
            %index contains the indices of the glacier-grid
            %centroids (including frame) contained with that main 
            %grid cell

        sIceInit.ic2igrd = cell(szMain);
        sIceInit.igrd2ic = nan(szIce(1),szIce(2), 'single');

        for ll = 1 : numel(sHydro.dem)
            [rC, cC] = ind2sub(szMain, ll);

            rMatch = find(sIceInit.igrdlat < edgeLatMain(rC) & sIceInit.igrdlat > edgeLatMain(rC+1));
            cMatch = find(sIceInit.igrdlon > edgeLonMain(cC) & sIceInit.igrdlon < edgeLonMain(cC+1));

            %Only proceed is glacier grid cells possibly exist
            %within current main grid cell
            if ~isempty(rMatch) && ~isempty(cMatch)
               %Ensure column orientation (matters for 'combvec'):
               rMatch = rMatch(:)';
               cMatch = cMatch(:)';
               subCurr = combvec(rMatch,cMatch);
               indCurr = sub2ind(szIce, subCurr(1,:), subCurr(2,:));

               %glacier dem is sparse array. Therefore, find all
               %non-zero grid cells (assumes no ice at sea-level)
               indCurr = indCurr(sIceInit.igrddem(indCurr) ~= 0);

               %If glaciers (or frame around glaciers) exist in
               %current grid cell, record indices
               if ~isempty(indCurr)
                   %Sparse array in which the value indicates which main
                   %grid cell the centroid of the ice grid cell is within
                   sIceInit.igrd2ic(indCurr) = ll;

                   %Cell array (inherently sparse) that has same
                   %resolution as main grid and contains a list of any
                   %glacier cells within that grid:
                   sIceInit.ic2igrd{ll} = indCurr;
%                    icBlMain(ll) = 1;
               end
            end
        end

        %Initialize arrays to track changes in ice depth
        sIceInit.icdwe = sparse(szMain(1),szMain(2));
        sIceInit.icdwe = sparse(szIce(1),szIce(2));
        
        save(pathIceStruct, 'sIceInit');
    end
else %In this case, ice grid same as main grid
    sIceInit.igrdwe  = sparse(szMain(1),szMain(2));
    sIceInit.icdwe = sparse(szMain(1),szMain(2));
    sIceInit.icx = zeros(szMain,'single');
end

szIce = size(sIceInit.igrddem);
szMain = size(sHydro.dem);

%Load debris thickness grid (if present):
if isfield(sPath, 'icdbr')
    [foldDbrSv, fildDbrSv, ~] = fileparts(sPath.icdbr);
    
    pathDbrSv = fullfile(foldDbrSv, ...
        [fildDbrSv, '_', sMeta.region{sMeta.siteCurr}, '_', num2str(szMain(1)) , 'x', num2str(szMain(2)), '.mat']);
    
    if exist(pathDbrSv, 'file')
        load(pathDbrSv);
    else
        sDebris = read_geodata_v2(sPath.icdbr, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
            sDebris.data = squeeze(sDebris.data);
        if ~isequal(size(sDebris.data), size(sIceInit.igrddem))
            resDebris = [mean(abs(diff(sDebris.(varLat)))), mean(abs(diff(sDebris.(varLon))))];
            resMain = [mean(abs(diff(sHydro.(varLat)))), mean(abs(diff(sHydro.(varLon))))];
            resIceGrd = [mean(abs(diff(sIceInit.igrdlat))), mean(abs(diff(sIceInit.igrdlon)))];

            prec = -order(min(resDebris))+1;
            if isequal(round2(resDebris, prec), round2(resIceGrd, prec)) || isequal(round2(resDebris, prec), round2(resMain, prec)) || all(resDebris < resMain)
                sDebris.data = geodata_area_wgt(sDebris.(varLon), sDebris.(varLat), sDebris.data, sHydro.(varLon), sHydro.(varLat), 'nansum');
            else
                error('iceGridLoaded:mismatchIceDebris', ...
                    'The ice and debris cover grids are not the same size.');
            end
        end
        
        save(pathDbrSv, 'sDebris');
    end

    sIceInit.icdbr = sDebris.data;
%     sIceInit.icdbr = sparse(szIce(1),szIce(2));
%     sIceInit.icdbr(indIce) = sDebris.data(indIce);
    

    disp(['Debris cover thickness grid loaded from ' sPath.icdbr '. This grid uses the main hydrologic grid.']);
%     warning('ice_grid_init:debrisCover',['Debris cover thickness grid loaded. Ensure this '...
%         'is used in the ice processs representations '...
%         '(including heat calculation).']);
end

%Load glacier lake fraction grid (if present)
if isfield(sPath, 'icpndx')
    [foldPndSv, fildPndSv, ~] = fileparts(sPath.icpndx);
    
    pathPndSv = fullfile(foldPndSv, ...
        [fildPndSv, '_', sMeta.region{sMeta.siteCurr}, '_', num2str(szMain(1)) , 'x', num2str(szMain(2)), '.mat']);
    
    if exist(pathPndSv, 'file')
        load(pathPndSv);
    else  
        sPond = read_geodata_v2(sPath.icpndx, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
            sPond.data = squeeze(sPond.data);
        if ~isequal(size(sPond.data), size(sIceInit.igrddem))
            resPndx = [mean(abs(diff(sPond.(varLat)))), mean(abs(diff(sPond.(varLon))))];
            resMain = [mean(abs(diff(sHydro.(varLat)))), mean(abs(diff(sHydro.(varLon))))];
            resIceGrd = [mean(abs(diff(sIceInit.igrdlat))), mean(abs(diff(sIceInit.igrdlon)))];

            prec = -order(min(resPndx))+1;
            if isequal(round2(resPndx, prec), round2(resIceGrd, prec)) || isequal(round2(resPndx, prec), round2(resMain, prec)) || all(resPndx < resMain)
                sPond.data = geodata_area_wgt(sPond.(varLon), sPond.(varLat), sPond.data, sHydro.(varLon), sHydro.(varLat), 'nansum');
            else
                error('CCHF_backbone:mismatchIcePond', ...
                    'The ice and pond fraction grids are not the same size.');
            end
        end
        
        save(pathPndSv, 'sPond');
    end

    sIceInit.icpndx = sPond.data;
    sIceInit.icpndx(isnan(sIceInit.icpndx)) = 0;
%     sIceInit.icdbr = sparse(szIce(1),szIce(2));
%     sIceInit.icdbr(indIce) = sDebris.data(indIce);

    disp(['Ice pond fraction loaded from ' sPath.icpndx '. This grid uses the main hydrologic grid.']);
%     warning('ice_grid_init:icPond',['Ice pond fraction loaded. Ensure this '...
%         'is actually used in the ice processs representations '...
%         '(including heat calculation).'])
end

%Load ice thickness (water equivalent) grid (if present)
if isfield(sPath, 'icwe')
    if isfield(sPath, 'icwe')
        [foldIcweSv, fildIcweSv, ~] = fileparts(sPath.icwe);
    elseif isfield(sPath, 'igrdwe')
        [foldIcweSv, fildIcweSv, ~] = fileparts(sPath.igrdwe);
    end
    
        
    pathIcweSv = fullfile(foldIcweSv, ...
        [fildIcweSv, '_', sMeta.region{sMeta.siteCurr}, '_', num2str(szIce(1)) , 'x', num2str(szIce(2)), '.mat']);
    
    if exist(pathIcweSv, 'file')
        load(pathIcweSv);
    else 
        if isfield(sPath, 'icwe')
            sIceWe = read_geodata_v2(sPath.icwe, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
        elseif isfield(sPath, 'igrdwe')
            sIceWe = read_geodata_v2(sPath.igrdwe, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
        end
        sIceWe.data = squeeze(sIceWe.data);
        
        if ~isequal(size(sIceWe.data), size(sIceInit.igrddem))
            resIceWe  = [mean(abs(diff( sIceWe.(varLat)))), mean(abs(diff( sIceWe.(varLon))))];
            resMain   = [mean(abs(diff( sHydro.(varLat)))), mean(abs(diff( sHydro.(varLon))))];
            resIceGrd = [mean(abs(diff(sIceInit.igrdlat))), mean(abs(diff(sIceInit.igrdlon)))];

            prec = -order(min(resIceWe))+1;
            if isequal(round2(resIceWe, prec), round2(resIceGrd, prec)) || isequal(round2(resIceWe, prec), round2(resMain, prec)) || all(resIceWe < resMain)
                sIceWe.data = geodata_area_wgt(sIceWe.(varLon), sIceWe.(varLat), sIceWe.data, sIceInit.igrdlon, sIceInit.igrdlat, 'nansum');
            else
                error('CCHF_backbone:mismatchIceWe', ...
                    'The ice and ice WE grids are not the same size.');
            end
        end
        
        sIceWe.data(isnan(sIceWe.data)) = 0;
        
        save(pathIcweSv, 'sIceWe');
    end

    %Populate sparse matrix:
    sIceInit.igrdwe = sparse(double(sIceWe.data));
    
    disp(['Ice water equivalent thickness loaded from ' ...
        fullfile(foldIcweSv, fildIcweSv) '. This grid uses the potentially finer resolution igrd.']);
%     warning('ice_grid_init:icWe',['Ice water equivalent loaded. Ensure this '...
%         'is actually used in the ice processs representations '...
%         '(including heat calculation).'])
end

disp('Finished ice grid initialization');