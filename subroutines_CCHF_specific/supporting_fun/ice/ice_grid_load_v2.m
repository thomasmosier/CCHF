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

function [sIceInit, icBlMain] = ice_grid_load_v2(sHydro, sPath, sMeta)

szDem = size(sHydro.dem); 
icBlMain = zeros(szDem, 'single');

%Load glacier coverage data and convert to depth:
%If glaciers or boolean clipping being used, create grids:
if isfield(sPath, 'ice') %&& ~regexpbl(sMeta.mode, 'parameter')

    %Create path to save ice structure to:
    [foldIce, fileIce, ~] = fileparts(sPath.ice);
%     [~, fileDem, ~] = fileparts(sPath.dem);
    pathIceStruct = fullfile(foldIce, [fileIce, '_', sMeta.region{sMeta.siteCurr}, '.mat']);
    
    if exist(pathIceStruct, 'file')
       load(pathIceStruct);
    else
        %Ice grid can either be 'same' as the main model grid or can be
        %'fine'. The advtange of 'fine' is that it better captures
        %local glacier dynamics
        
        disp(['Loading ice grid. This may take quite some time (a few minutes for grids with 200 indices, but a few hours for very large domains) ' ...
            'depending on the size of the domain and the options chosen.' ...
            char(10) 'Output will automatically be saved to skip this '...
            'process in subsequent runs for this model setup.']);

        edgeLat = box_edg(sHydro.lat);
        edgeLon = box_edg(sHydro.lon);


        %Use edg to crop ice dem:
        if regexpbl(sMeta.iceGrid, 'fine')
            sIceDem = read_geodata_v2(sPath.iceDem, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
                sIceDem.data = squeeze(sIceDem.data);

            %FORCE EDGES OF FINE GRID TO ALIGN WITH EDGES OF MAIN GRID:
            %Check lon/lat alignment (display warning if fine grid does not
            %align with main grid)
            edgLonIceIn = box_edg(sIceDem.lon);
            deltaLonAlign = nan(numel(edgeLon),1);
            edgLatIceIn = box_edg(sIceDem.lat);
            deltaLatAlign = nan(numel(edgeLat),1);

            for ii = 1 : numel(deltaLonAlign)
                deltaLonAlign(ii) = min(abs(edgLonIceIn - edgeLon(ii)));
            end
            for ii = 1 : numel(deltaLatAlign)
                deltaLatAlign(ii) = min(abs(edgLatIceIn - edgeLat(ii)));
            end

            if any(order(deltaLonAlign) > order(abs(nanmean(diff(sIceDem.lon)))) - 2)
                warning('ice_grid_init:diffLonGrid',['The longitude edges of the fine '...
                    'grid cells do not align with those of the main grid.'...
                    char(10) 'The fine grid is going to be resampled to '...
                    'align with the main grid.' char(10) 'This introduces '...
                    'error into the ice grid calculations.']);
            end
            if any(order(deltaLatAlign) > order(abs(nanmean(diff(sIceDem.lat)))) - 2)
                warning('ice_grid_init:diffLatGrid',['The latitude edges of the fine '...
                    'grid cells do not align with those of the main grid.'...
                    char(10) 'The fine grid is going to be resampled to '...
                    'align with the main grid.' char(10) 'This introduces '...
                    'error into the ice grid calculations.']);
            end

            %Find scale factor between two grids:
            nRefin = round(nanmean([abs(nanmean(diff(sHydro.lon)))/abs(nanmean(diff(sIceDem.lon))), ...
                abs(nanmean(diff(sHydro.lat)))/abs(nanmean(diff(sIceDem.lat)))]));
            %Create output grid:
            edgLonOut = nan(nRefin*numel(sHydro.lon) + 1, 1);
            edgLatOut = nan(nRefin*numel(sHydro.lat) + 1, 1);

            for ii = 1 : numel(sHydro.lon)
                lonTemp = linspace(edgeLon(ii), edgeLon(ii+1), nRefin+1);
                if ii == numel(sHydro.lon)
                    edgLonOut((ii-1)*nRefin + 1 : ii*nRefin + 1) = lonTemp;
                else
                    edgLonOut((ii-1)*nRefin + 1 : ii*nRefin) = lonTemp(1:end-1);
                end
            end
            for ii = 1 : numel(sHydro.lat)
                latTemp = linspace(edgeLat(ii), edgeLat(ii+1), nRefin+1);
                if ii == numel(sHydro.lat)
                    edgLatOut((ii-1)*nRefin + 1 : ii*nRefin + 1) = latTemp;
                else
                    edgLatOut((ii-1)*nRefin + 1 : ii*nRefin) = latTemp(1:end-1);
                end
            end

            [lonInMesh, latInMesh] = meshgrid(sIceDem.lon, sIceDem.lat);

            sIceDem.lon = edgLonOut(1 : end - 1) + 0.5*diff(edgLonOut);
            sIceDem.lat = edgLatOut(1 : end - 1) + 0.5*diff(edgLatOut);

            [lonOutMesh, latOutMesh] = meshgrid(sIceDem.lon, sIceDem.lat);

            sIceDem.data = interp2(lonInMesh, latInMesh, sIceDem.data, ...
                lonOutMesh, latOutMesh, 'linear');

        elseif regexpbl(sMeta.iceGrid, 'same')
            sIceDem.data = squeeze(sHydro.dem);
            sIceDem.lon = sHydro.lon;
            sIceDem.lat = sHydro.lat;
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
            [iceExist, iceArea] = shapefile_geo_v2(shpIce, strIceLd, sIceDem.lat, sIceDem.lon);
%             iceExist = interpshapefile_geo_v2(shpIce, sIceDem.lat, sIceDem.lon, strIceLd, 1);
            indIce = find(~isnan(iceExist));
%             disp('finished first ice shp interp')
        else
            sIceInput = read_geodata_v2(sPath.ice, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
                sIceInput.data = squeeze(sIceInput.data);
            %Check that coordinates are same as iceDEM:
            lonRes = log10(round(1/nanmean(abs(diff(sIceInput.lon)))))+1;
            latRes = log10(round(1/nanmean(abs(diff(sIceInput.lat)))))+1;
            if numel(sIceInput.lon) ~= numel(sIceDem.lon) || numel(sIceInput.lat) ~= numel(sIceDem.lat) 
                error('CCHF_backbone:diffIceres','The ice DEM and ice presence grid have different sizes but the option has been state such that the grids are expected to have the same size.');
            elseif sum(ismember(round2(sIceInput.lon, lonRes), round2(sIceDem.lon, lonRes))) < 0.9*numel(sIceDem.lon) || sum(ismember(round2(sIceInput.lat, latRes), round2(sIceDem.lat, latRes))) < 0.9*numel(sIceDem.lat)
               error('CCHF_backbone:misalignedIceDem','The ice DEM and ice presence grid do not appear to align. Maybe rounding error?');
            end

            %Grids basicaly align, so ensure they match completely.
            [indLonIcePre, indLonIceDem] = ismember(round2(sIceInput.lon, lonRes), round2(sIceDem.lon, lonRes));
            [indLatIcePre, indLatIceDem] = ismember(round2(sIceInput.lat, latRes), round2(sIceDem.lat, latRes));

            if all(indLonIcePre == 1 | indLonIcePre == 0)
                indLonIcePre = find(indLonIcePre == 1);
                indLatIcePre = find(indLatIcePre == 1);
            end
            if all(indLonIceDem == 1 | indLonIceDem == 0)
                indLonIceDem = find(indLonIceDem == 1);
                indLatIceDem = find(indLatIceDem == 1);
            end

            sIceInput.data = sIceInput.data(indLatIcePre, indLonIcePre);
            sIceInput.lon  = sIceInput.lon(indLonIcePre);
            sIceInput.lat  = sIceInput.lat(indLatIcePre);

            sIceDem.data = sIceDem.data(indLatIceDem, indLonIceDem);
            sIceDem.lon  = sIceDem.lon(indLonIceDem);
            sIceDem.lat  = sIceDem.lat(indLatIceDem);

            indIce = find(~isnan(sIceInput.data) & sIceInput.data ~= 0);
        end

        %Populate ice presence grid:
        szIce = size(sIceDem.data);
        %Initialize ice array 
        sIceInit = struct;
        %Initialize water equivalent as sparse:
        sIceInit.icwe = sparse(szIce(1), szIce(2));

        if regexpbl(sPath.ice,'.shp')
            sIceInit.icwe(indIce) = 1;
        else
            sIceInit.icwe(indIce) = sIceInput.data(indIce);
        end
        
        %Record ice area (if loaded)
        if regexpbl(sPath.ice,'.shp') && regexpbl(strIceLd, 'Area')
            sIceInit.icgrdarea = iceArea;
        end


        %Estimate fractional glacier extent within each grid cell
        subDiv = 5;
        lonSub = nan(subDiv*numel(sHydro.lon),1);
        lonSubInd = cell(numel(sHydro.lon),1);
        for ll = 1 : numel(sHydro.lon)
            cntrSub = (ll-1)*(subDiv) + 1;
            lonSubInd{ll} = (cntrSub : cntrSub + subDiv-1);

            lonSubCurr = linspace(edgeLon(ll), edgeLon(ll+1), subDiv+1);
            lonSubCurr = lonSubCurr(1:subDiv) + 0.5*(mean(abs(diff(lonSubCurr(1:subDiv)))));
            lonSub(lonSubInd{ll}) = lonSubCurr;
        end
        latSub = nan(subDiv*numel(sHydro.lat),1);
        latSubInd = cell(numel(sHydro.lat),1);
        for ll = 1 : numel(sHydro.lat)
            cntrSub = (ll-1)*(subDiv) + 1;
            latSubInd{ll} = (cntrSub : cntrSub + subDiv-1);

            latSubCurr = linspace(edgeLat(ll), edgeLat(ll+1), subDiv+1);
            latSubCurr = latSubCurr(1:subDiv) - 0.5*(mean(abs(diff(latSubCurr(1:subDiv)))));
            latSub(latSubInd{ll}) = latSubCurr;
        end

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

            tempIceBl = zeros(size(sIceInit.icwe),'single');
            tempIceBl(full(sIceInit.icwe ~= 0) & full(~isnan(sIceInit.icwe))) = 1;
            [lonMesh,latMesh] = meshgrid(sIceDem.lon, sIceDem.lat);
            [lonSubMesh,latSubMesh] = meshgrid(lonSub, latSub);
            iceSubExist = interp2(lonMesh,latMesh, tempIceBl, lonSubMesh,latSubMesh,'nearest');
        end

%         disp('Estimating fractional ice area in main grid.')
        sIceInit.icx = zeros(numel(sHydro.lat),numel(sHydro.lon));

        for ll = 1 : numel(sHydro.lat)
            for mm = 1 : numel(sHydro.lon)
                sIceInit.icx(ll,mm) = sum2d(iceSubExist(latSubInd{ll},lonSubInd{mm}))/100;
            end
        end
%         disp('finished fractional ice shp interp')

        %Note all main grid cells that contain any ice
        sIceInit.icbl = sparse(szDem(1),szDem(2));
        sIceInit.icbl(sIceInit.icx > 0) = 1;

        %Create sparse ice grid (including frame) and delete previous structures:                
        %Find indices of frame around existing ice (allows glacier morphology)

%         disp('Calculating surface flow properties for ice grid.');
        
        %Set # of frame elements based on size (want approximately 1
        %km)
            rEarth = find_att(sMeta.global,'radius_Earth');
        frIce = ceil(1000/haversine_reference(sIceDem.lon(1),sIceDem.lat(1),sIceDem.lon(2),sIceDem.lat(2), rEarth));
            frIce = max(frIce,2);
        %In for loop, set all indices within frame to 0. After for loop, find indices of all 0 elements.    
        [rIce, cIce] = ind2sub(szIce, indIce);
        indDimIce = reshape((1:(szIce(1)*szIce(2))), szIce);

        for ll = 1 : numel(indIce)
            rIceCurr = (rIce(ll)-frIce:rIce(ll)+frIce);
                rIceCurr(rIceCurr < 1) = [];
                rIceCurr(rIceCurr > szIce(1)) = [];
            cIceCurr = (cIce(ll)-frIce:cIce(ll)+frIce);
                cIceCurr(cIceCurr < 1) = [];
                cIceCurr(cIceCurr > szIce(2)) = [];

            indDimIce(rIceCurr, cIceCurr) = 0;
        end

        indIceFr = find(indDimIce == 0);

        sIceInit.icgrddem =  sparse(szIce(1), szIce(2));
        sIceInit.icgrddem(indIceFr) = sIceDem.data(indIceFr);
        sIceInit.icgrdlat = sIceDem.lat;
        sIceInit.icgrdlon = sIceDem.lon;
            clear('sIceInput','sGlaciers');
            clear('sIceDem');


        %Calculate ice surface slope and flow direction (based on surface elevation)
        [~, sIceInit.icgrdfdr, iceSurfSlope] = ...
            wflowacc(sIceInit.icgrdlon, sIceInit.icgrdlat, sIceInit.icgrddem,'type','single','edges','open','coord','geographic');

        %Calculate slope in the flow direction:
        sIceInit.icgrdslopefdr = -sparse(reshape(sum(sIceInit.icgrdfdr.*iceSurfSlope, 2), szIce));
            sIceInit.icgrdslopefdr(sIceInit.icgrddem == 0) = 0;
        %Calculate slope in all directions:
        sIceInit.icgrdslope = sparse(iceSurfSlope);
            sIceInit.icgrdslope(sIceInit.icgrddem == 0, :) = 0;
            sIceInit.icgrdslope(:, sIceInit.icgrddem == 0) = 0;


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

        sIceInit.main2ice = cell(szDem);
        sIceInit.ice2main = sparse(szIce(1),szIce(2));

        for ll = 1 : numel(sHydro.dem)
            [rC, cC] = ind2sub(szDem, ll);

            rMatch = find(sIceInit.icgrdlat < edgeLat(rC) & sIceInit.icgrdlat > edgeLat(rC+1));
            cMatch = find(sIceInit.icgrdlon > edgeLon(cC) & sIceInit.icgrdlon < edgeLon(cC+1));

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
               indCurr = indCurr(sIceInit.icgrddem(indCurr) ~= 0);

               %If glaciers (or frame around glaciers) exist in
               %current grid cell, record indices
               if ~isempty(indCurr)
                   %Sparse array in which the value indicates which main
                   %grid cell the centroid of the ice grid cell is within
                   sIceInit.ice2main(indCurr) = ll;

                   %Cell array (inherently sparse) that has same
                   %resolution as main grid and contains a list of any
                   %glacier cells within that grid:
                   sIceInit.main2ice{ll} = indCurr;
                   icBlMain(ll) = 1;
               end
            end
        end

        %Initialize array to track changes in ice depth
        sIceInit.icdwe = sparse(szDem(1),szDem(2));
        
        save(pathIceStruct, 'sIceInit', 'icBlMain');
    end
else %In this case, ice grid same as main grid
    sIceInit.icwe  = sparse(szDem(1),szDem(2));
    sIceInit.icdwe = sparse(szDem(1),szDem(2));
    sIceInit.icx = zeros(szDem,'single');
end


varLat = 'latitude';
varLon = 'longitude';
%Load debris thickness grid (if present):
if isfield(sPath, 'icdbr')
    sDebris = read_geodata_v2(sPath.icdbr, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
        sDebris.data = squeeze(sDebris.data);
    if ~isequal(size(sDebris.data), size(sIceInit.icgrddem))
        resDebris = [mean(abs(diff(sDebris.(varLat)))), mean(abs(diff(sDebris.(varLon))))];
        resIceDem = [mean(abs(diff(sIceInit.icgrdlat))), mean(abs(diff(sIceInit.icgrdlon)))];
        
        prec = -order(min(resDebris))+1;
        if isequal(round2(resDebris, prec), round2(resIceDem, prec)) || all(resDebris < resIceDem)
            sDebris.data = geodata_area_wgt(sDebris.(varLon), sDebris.(varLat), sDebris.data, sIceInit.icgrdlon, sIceInit.icgrdlat, 'nansum');
        else
            error('iceGridLoaded:mismatchIceDebris', ...
                'The ice and debris cover grids are not the same size.');
        end
    end

    sIceInit.icdbr = sDebris.data;
%     sIceInit.icdbr = sparse(szIce(1),szIce(2));
%     sIceInit.icdbr(indIce) = sDebris.data(indIce);

    disp(['Debris cover thickness grid loaded from ' sPath.icdbr]);
%     warning('ice_grid_init:debrisCover',['Debris cover thickness grid loaded. Ensure this '...
%         'is used in the ice processs representations '...
%         '(including heat calculation).']);
end

%Load glacier lake fraction grid (if present)
if isfield(sPath, 'icpndx')
    sPond = read_geodata_v2(sPath.icpndx, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
        sPond.data = squeeze(sPond.data);
    if ~isequal(size(sPond.data), size(sIceInit.icgrddem))
        resDebris = [mean(abs(diff(sPond.(varLat)))), mean(abs(diff(sPond.(varLon))))];
        resIceDem = [mean(abs(diff(sIceInit.icgrdlat))), mean(abs(diff(sIceInit.icgrdlon)))];
        
        prec = -order(min(resDebris))+1;
        if isequal(round2(resDebris, prec), round2(resIceDem, prec)) || all(resDebris < resIceDem)
            sPond.data = geodata_area_wgt(sPond.(varLon), sPond.(varLat), sPond.data, sIceInit.icgrdlon, sIceInit.icgrdlat, 'nansum');
        else
            error('CCHF_backbone:mismatchIcePond', ...
                'The ice and pond fraction grids are not the same size.');
        end
    end

    sIceInit.icpndx = sPond.data;
    sIceInit.icpndx(isnan(sIceInit.icpndx)) = 0;
%     sIceInit.icdbr = sparse(szIce(1),szIce(2));
%     sIceInit.icdbr(indIce) = sDebris.data(indIce);

    disp(['Ice pond fraction loaded from ' sPath.icpndx]);
%     warning('ice_grid_init:icPond',['Ice pond fraction loaded. Ensure this '...
%         'is actually used in the ice processs representations '...
%         '(including heat calculation).'])
end

%Load ice thickness (water equivalent) grid (if present)
if isfield(sPath, 'icwe')
    sIceWe = read_geodata_v2(sPath.icwe, 'data', nan(1,2), nan(1,2), nan(1,2), 0, 'out', 'none','no_disp', 'onefile');
        sIceWe.data = squeeze(sIceWe.data);
    if ~isequal(size(sIceWe.data), size(sIceInit.icgrddem))
        resIceWe = [mean(abs(diff(sIceWe.(varLat)))), mean(abs(diff(sIceWe.(varLon))))];
        resIceDem = [mean(abs(diff(sIceInit.icgrdlat))), mean(abs(diff(sIceInit.icgrdlon)))];
        
        prec = -order(min(resIceWe))+1;
        if isequal(round2(resIceWe, prec), round2(resIceDem, prec)) || all(resIceWe < resIceDem)
            sIceWe.data = geodata_area_wgt(sIceWe.(varLon), sIceWe.(varLat), sIceWe.data, sIceInit.icgrdlon, sIceInit.icgrdlat, 'nansum');
        else
            error('CCHF_backbone:mismatchIceWe', ...
                'The ice and ice WE grids are not the same size.');
        end
    end

    sIceInit.icwe = sIceWe.data;
%     sIceInit.icdbr = sparse(szIce(1),szIce(2));
%     sIceInit.icdbr(indIce) = sDebris.data(indIce);

    disp(['Ice water equivalent thickness loaded from ' sPath.icwe]);
%     warning('ice_grid_init:icWe',['Ice water equivalent loaded. Ensure this '...
%         'is actually used in the ice processs representations '...
%         '(including heat calculation).'])
end

disp('Finished ice grid initialization');