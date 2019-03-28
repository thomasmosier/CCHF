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

if isfield(sHydro, 'longitude')
    varLon = 'longitude';
elseif isfield(sHydro, 'lon')
    varLon = 'lon';    
else
    error('iceGridLoad:noLonField','No longitude field found in sHydro.')
end

if isfield(sHydro, 'latitude')
    varLat = 'latitude';
elseif isfield(sHydro, 'lat')
    varLat = 'lat';    
else
    error('iceGridLoad:noLatField','No latitude field found in sHydro.')
end

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

        %Get dimensions of ice grid:
        szIce = size(sIceDem.dem);
        
        %Initialize ice structure array output:
        sIceInit = struct;
        

        %Define grid that is finer than main grid for estimating fractional 
        %glacier extent within each main grid cell
        if regexpbl(sMeta.iceGrid, 'fine')
            [lonMnSub, latMnSub] = meshgrid(sIceDem.(varLon), sIceDem.(varLat));
        else
            nScl = 5;
            [lonMnSub, latMnSub] = grid_refine(edgeLonMain, edgeLatMain, nScl, 'edge');
        end

        
        %Load ice presence file (shp or grid)
%         disp('Interpolating ice shapefile to finer mesh in order to estimate fractional ice area in main grid.');
        if regexpbl(sPath.ice,'.shp')
            %Read ice shapefile:
            shpIce = m_shaperead(sPath.ice);
            shpIce = shp_MMap2Matlab(shpIce, {'Area_2','Area'},'Area');
%             %Built-in function from mapping toolbox:
%             shpIce = shaperead(sPath.ice);

            %Find ice presence at sub-grid scale (used for estimating ice fraction)
            %This part eats a lot of memory. Consider alternative
            %strategies
            iceSubExist = shapefile_geo_v2(shpIce, latMnSub, lonMnSub, 1);
        else
            %Read ice grid:
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

            %Ensure grids exctly match
            %Find indices that match
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
            %Crop based on indices:
            sIceInput.data = sIceInput.data(indLatIcePre, indLonIcePre);
            sIceInput.(varLon)  = sIceInput.(varLon)(indLonIcePre);
            sIceInput.(varLat)  = sIceInput.(varLat)(indLatIcePre);

            sIceDem.dem = sIceDem.dem(indLatIceDem, indLonIceDem);
            sIceDem.(varLon)  = sIceDem.(varLon)(indLonIceDem);
            sIceDem.(varLat)  = sIceDem.(varLat)(indLatIceDem);
            
            warning('CCHF_backbone:iceResample',['The input ice grid '...
                'is being resampled for the purpose of estimating ' ...
                'fractional ice coverage.' char(10) 'This uncertainty can be '...
                'avoided through using a shapefile instead.' char(10) 'It '...
                'should not be a problem if the glacier grid is '...
                'significantly finer resolution than the main grid.']);

            %Find ice presence at sub-grid scale (used for estimating ice fraction)
            tempIceBl = zeros(szIce,'single');
            tempIceBl(~isnan(sIceInput.data) & sIceInput.data ~= 0) = 1;
            [lonMesh,latMesh] = meshgrid(sIceDem.(varLon), sIceDem.(varLat));
            [lonSubMesh,latSubMesh] = meshgrid(lonMnSub, latMnSub);
            iceSubExist = interp2(lonMesh,latMesh, tempIceBl, lonSubMesh,latSubMesh,'nearest');
        end
    

        %Estimate fractional ice cover (main grid resolution):
%         disp('Estimating fractional ice area in main grid.')
        %Calculate scaling factor:
        nSclSq = nScl^2;
        %Initialize grid
        sIceInit.icx = sparse(szMain(1), szMain(2));
        %Loop over sub-grid to find fractional area
        for ll = 1 : numel(sHydro.(varLat))
            indLatN = find(latMnSub <= edgeLatMain(ll), 1, 'first');
            indLatS = find(edgeLatMain(ll+1) <= latMnSub, 1, 'last');
    
            for mm = 1 : numel(sHydro.(varLon))
                indLonW = find(edgeLonMain(mm) <= lonMnSub, 1, 'first');
                indLonE = find(lonMnSub <= edgeLonMain(mm+1), 1, 'last');
    
                sIceInit.icx(ll,mm) = sum2d(iceSubExist(indLatN:indLatS,indLonW:indLonE))/nSclSq;
            end
        end
        

        %Initialize main ice water equivalent grid as sparse array:
        sIceInit.igrdwe = sparse(szIce(1), szIce(2));
        
        if regexpbl(sMeta.iceGrid, 'fine')
            warning('iceGridLoad:checkFineIcePresence','Check the following section of code that finds which grid cells contain ice for the fine ice grid model state.')
            edgeLatIce = box_edg(sIceDem.(varLat));
            edgeLonIce = box_edg(sIceDem.(varLon));
               
            %Loop over sub-grid to find fractional area
            for ll = 1 : numel(sIceDem.(varLat))
                indLatN = find(sIceDem.(varLat) <= edgeLatIce(ll), 1, 'first');
                indLatS = find(edgeLatIce(ll+1) <= sIceDem.(varLat), 1, 'last');

                for mm = 1 : numel(sIceDem.(varLon))
                    indLonW = find(edgeLonIce(mm) <= sIceDem.(varLon), 1, 'first');
                    indLonE = find(sIceDem.(varLon) <= edgeLonIce(mm+1), 1, 'last');

                    sIceInit.igrdwe(ll,mm) = any2d(iceSubExist(indLatN:indLatS,indLonW:indLonE));
                end
            end
        else
            %Initalize to 1 where there is ice
            sIceInit.igrdwe(sIceInit.icx > 0) = 1;
        end
        
%         disp('finished fractional ice shp interp')

        %Assign DEM information to ice initialization structure
        sIceInit.igrddem = sIceDem.dem;
        sIceInit.igrdlat = sIceDem.(varLat);
        sIceInit.igrdlon = sIceDem.(varLon);
            clear('sIceInput','sGlaciers');
            clear('sIceDem');


        %Calculate ice surface slope and flow direction (based on surface elevation)
        [~, sIceInit.igrdfdr, iceSurfSlope] = ...
            wflowacc(sIceInit.igrdlon, sIceInit.igrdlat, sIceInit.igrddem,'type','single','edges','open','coord','geographic');

        %Calculate slope in the flow direction:
        %SLOPE = rise / run (unitless; and not angle)
        sIceInit.igrdslopefdr = -sparse(reshape(sum(sIceInit.igrdfdr.*iceSurfSlope, 2), szIce));
%             sIceInit.igrdslopefdr(sIceInit.igrddem == 0) = 0;
        %Calculate slope in all directions:
        sIceInit.igrdslope = sparse(iceSurfSlope);
%             sIceInit.igrdslope(sIceInit.igrddem == 0, :) = 0;
%             sIceInit.igrdslope(:, sIceInit.igrddem == 0) = 0;


        %%Find correspondance between glacier grid and main grid
        %Check which of main grid cell the centroid of the glacier
        %grid cell is contained within

        %Initialize two grids:
        %sparse grid = tracks which main grid cells contain 
            %glaciers 
        %cell array = Dim. of main grid; each non-empty cell-array
            %index contains the indices of the glacier-grid
            %centroids (including frame) contained with that main 
            %grid cell
        sIceInit.main2igrd = cell(szMain);
        sIceInit.igrd2main = nan(szIce(1),szIce(2), 'single');

        %Populate arrays containing correspondance between main and ice
        %grids
        for ll = 1 : numel(sHydro.dem)
            [rC, cC] = ind2sub(szMain, ll);

            rMatch = find(sIceInit.igrdlat < edgeLatMain(rC) & sIceInit.igrdlat > edgeLatMain(rC+1));
            cMatch = find(sIceInit.igrdlon > edgeLonMain(cC) & sIceInit.igrdlon < edgeLonMain(cC+1));

            %Only proceed is glacier grid cells possibly exist
            %within current main grid cell
            if ~isempty(rMatch) && ~isempty(cMatch)
               %Ensure column orientation (matters for 'combvec'/'allcomb'):
               rMatch = rMatch(:)';
               cMatch = cMatch(:)';
               
               %Get all combinations using function from File Exchange
               subCurr = allcomb(rMatch, cMatch);
               indCurr = sub2ind(szIce, subCurr(:,1), subCurr(:,2));
               
               %Build in function requiring Neural Net Toolbox
%                subCurr = combvec(rMatch,cMatch);
%                indCurr = sub2ind(szIce, subCurr(1,:), subCurr(2,:));

               %glacier dem is sparse array. Therefore, find all
               %non-zero grid cells (assumes no ice at sea-level)
               indCurr = indCurr(sIceInit.igrddem(indCurr) ~= 0);

               %If glaciers (or frame around glaciers) exist in
               %current grid cell, record indices
               if ~isempty(indCurr)
                   %Sparse array in which the value indicates which main
                   %grid cell the centroid of the ice grid cell is within
                   sIceInit.igrd2main(indCurr) = ll;

                   %Cell array (inherently sparse) that has same
                   %resolution as main grid and contains a list of any
                   %glacier cells within that grid:
                   sIceInit.main2igrd{ll} = indCurr;
%                    icBlMain(ll) = 1;
               end
            end
        end

        %Initialize arrays to track changes in ice depth
        sIceInit.icdwe = sparse(szMain(1),szMain(2));
        sIceInit.igrddwe = sparse( szIce(1), szIce(2));
        
        save(pathIceStruct, 'sIceInit');
    end
else %In this case, ice grid same as main grid
    sIceInit.igrddem = sHydro.dem; %Ice grid DEM
    sIceInit.igrdlat = sHydro.(varLat);
    sIceInit.igrdlon = sHydro.(varLon);
    sIceInit.igrdfdr      = sHydro.fdr;
    sIceInit.igrdslopefdr = sHydro.slopefdr;
    sIceInit.igrdslope    = sHydro.slope;
    sIceInit.igrdwe  = sparse(szMain(1),szMain(2)); %Ice grid water equivalent
    sIceInit.igrddwe = sparse(szMain(1),szMain(2)); %Ice grid change in water equivalent
    sIceInit.icdwe   = sparse(szMain(1),szMain(2)); %Main grid change in water equivalent
    sIceInit.icx     =  zeros(szMain,'single'); %Main grid fractional ice coverage
end


%Write ice grids to model folder
write_ESRI_v4(full(sIceInit.igrdwe), ESRI_hdr(sHydro.(varLon), sHydro.(varLat), 'corner'), fullfile(sPath.output, [sMeta.region{sMeta.siteCurr} '_glacier_initial_boolean_main_grid.asc'])         , 0);
write_ESRI_v4(full(sIceInit.icx   ), ESRI_hdr(sHydro.(varLon), sHydro.(varLat), 'corner'), fullfile(sPath.output, [sMeta.region{sMeta.siteCurr} '_glacier_initial_fractional_cover_main_grid.asc']), 2);


%Calculate grid sizes
szMain = size(sHydro.dem);
if isfield(sIceInit, 'igrddem')
    szIce = size(sIceInit.igrddem);
else
    szIce = szMain;
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
        load(pathIcweSv, 'iGrdWeTemp', 'icWeTemp');
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
                iGrdWeTemp = geodata_area_wgt(sIceWe.(varLon), sIceWe.(varLat), sIceWe.data, sIceInit.igrdlon, sIceInit.igrdlat, 'nansum');
            else
                error('CCHF_backbone:mismatchIceWe', ...
                    'The ice and ice WE grids are not the same size.');
            end
        else
            iGrdWeTemp = sIceWe.data;
        end

        if ~isequal(size(sIceWe.data), size(sHydro.dem))
            resIceWe  = [mean(abs(diff( sIceWe.(varLat)))), mean(abs(diff( sIceWe.(varLon))))];
            resMain   = [mean(abs(diff( sHydro.(varLat)))), mean(abs(diff( sHydro.(varLon))))];
            resIceGrd = [mean(abs(diff(sHydro.lat))), mean(abs(diff(sHydro.lon)))];

            prec = -order(min(resIceWe))+1;
            if isequal(round2(resIceWe, prec), round2(resIceGrd, prec)) || isequal(round2(resIceWe, prec), round2(resMain, prec)) || all(resIceWe < resMain)
                icWeTemp = geodata_area_wgt(sIceWe.(varLon), sIceWe.(varLat), sIceWe.data, sHydro.lon, sHydro.lat, 'nansum');
            else
                error('CCHF_backbone:mismatchIceWe', ...
                    'The ice and ice WE grids are not the same size.');
            end
        else
            icWeTemp = sIceWe.data;
        end
        
        iGrdWeTemp(isnan(iGrdWeTemp)) = 0;
        icWeTemp(isnan(icWeTemp)) = 0;
        
        save(pathIcweSv, 'iGrdWeTemp', 'icWeTemp');
    end

    %Populate matrices:
    sIceInit.icwe = icWeTemp;
    sIceInit.igrdwe = sparse(double(iGrdWeTemp));
    
    disp(['Ice water equivalent thickness loaded from ' ...
        fullfile(foldIcweSv, fildIcweSv) '. This grid uses the potentially finer resolution igrd.']);
%     warning('ice_grid_init:icWe',['Ice water equivalent loaded. Ensure this '...
%         'is actually used in the ice processs representations '...
%         '(including heat calculation).'])
end


%Load debris thickness grid (if present):
if isfield(sPath, 'icdbr')
    [foldDbrSv, fildDbrSv, ~] = fileparts(sPath.icdbr);
    
    pathDbrSv = fullfile(foldDbrSv, ...
        [fildDbrSv, '_', sMeta.region{sMeta.siteCurr}, '_', num2str(szMain(1)) , 'x', num2str(szMain(2)), '.mat']);
    
    if exist(pathDbrSv, 'file')
        load(pathDbrSv, 'sDebris');
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
        load(pathPndSv, 'sPond');
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


disp('Finished ice grid initialization');