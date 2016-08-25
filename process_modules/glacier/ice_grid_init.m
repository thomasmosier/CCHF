function icBlMain = ice_grid_init(sHydro, sPath, sMeta)

global sCryo
  
szDem = size(sHydro.dem); 
icBlMain = zeros(szDem, 'single');

%Load glacier coverage data and convert to depth:
%If glaciers or boolean clipping being used, create grids:
if isfield(sPath, 'ice') && ~regexpbl(sMeta.mode, 'parameter')

    %Ice grid can either be 'same' as the main model grid or can be
    %'fine'. The advtange of 'fine' is that it better captures
    %local glacier dynamics

    edgeLat = box_edg(sHydro.lat);
    edgeLon = box_edg(sHydro.lon);

    
    %Use edg to crop ice dem:
    if regexpbl(sMeta.iceGrid, 'fine')
        sIceDem = read_geodata(sPath.iceDem, sMeta, 'no_disp');
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
        shpIce = shaperead(sPath.ice);

        if isfield(shpIce,'Slope')
            strIceLd = 'Slope';
        elseif isfield(shpIce,'Aspect')
            strIceLd = 'Aspect';
        elseif isfield(shpIce,'Area')
            strIceLd = 'Area';
        else
            error('CCHF_backbone:unknownShpField',['A suitable '...
                'field for quirying the ice shapefile was not found.'])
        end
            warning('off','interpshapefile:overlap');
        iceExist = interpshapefile(shpIce, sIceDem.lat, sIceDem.lon, strIceLd);
        indIce = find(~isnan(iceExist));
    else
        sIceInput = read_geodata(sPath.ice, sMeta, 'no_disp');
            sIceInput.data = squeeze(sIceInput.data);
        %Check that coordinates are same as iceDEM:
        if sum(ismember(sIceInput.lon,sIceDem.lon)) < 0.9*numel(sIceDem.lon) || sum(ismember(sIceInput.lat,sIceDem.lat)) < 0.9*numel(sIceDem.lat)
           error('CCHF_backbone:misalignedIceDem','The ice DEM and ice presence grid do not appear to align. Maybe rounding error?');
        end

        %Grids basicaly align, so ensure they match completely.
        [indLonIcePre, indLonIceDem] = ismember(sIceInput.lon, sIceDem.lon);
        [indLatIcePre, indLatIceDem] = ismember(sIceInput.lat, sIceDem.lat);

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
    %Initialize ice array as sparse:
    sCryo.iceWE = sparse(szIce(1), szIce(2));

    if regexpbl(sPath.ice,'.shp')
        sCryo.iceWE(indIce) = 1;
    else
        sCryo.iceWE(indIce) = sIceInput.data(indIce);
    end
    
    
    %Estimate fractional glacier extent within each grid cell
    subDiv = 11;
    lonSub = nan(10*numel(sHydro.lon),1);
    lonSubInd = cell(numel(sHydro.lon),1);
    for ll = 1 : numel(sHydro.lon)
        cntrSub = (ll-1)*(subDiv-1) + 1;
        lonSubInd{ll} = (cntrSub : cntrSub + subDiv-2);

        lonSubCurr = linspace(edgeLon(ll), edgeLon(ll+1), subDiv);
        lonSubCurr = lonSubCurr(1:subDiv-1) + 0.5*(mean(abs(diff(lonSubCurr(1:subDiv-1)))));
        lonSub(lonSubInd{ll}) = lonSubCurr;
    end
    latSub = nan(10*numel(sHydro.lat),1);
    latSubInd = cell(numel(sHydro.lat),1);
    for ll = 1 : numel(sHydro.lat)
        cntrSub = (ll-1)*(subDiv-1) + 1;
        latSubInd{ll} = (cntrSub : cntrSub + subDiv-2);

        latSubCurr = linspace(edgeLat(ll), edgeLat(ll+1), subDiv);
        latSubCurr = latSubCurr(1:subDiv-1) - 0.5*(mean(abs(diff(latSubCurr(1:subDiv-1)))));
        latSub(latSubInd{ll}) = latSubCurr;
    end

    if regexpbl(sPath.ice,'.shp')
            warning('off','interpshapefile:overlap');
        iceSubExist = interpshapefile(shpIce, latSub, lonSub, strIceLd);
            iceSubExist(~isnan(iceSubExist)) = 1;
            iceSubExist(isnan(iceSubExist)) = 0;
    else
        warning('CCHF_backbone:iceResample',['The input ice grid '...
            'is being resampled for the purpose of estimating ' ...
            'fractional ice coverage.' char(10) 'This uncertainty can be '...
            'avoided through using a shapefile instead.' char(10) 'It '...
            'should not be a problem if the glacier grid is '...
            'significantly finer resolution than the main grid.']);

        tempIceBl = zeros(size(sCryo.iceWE),'single');
        tempIceBl(full(sCryo.iceWE ~= 0) & full(~isnan(sCryo.iceWE))) = 1;
        [lonMesh,latMesh] = meshgrid(sIceDem.lon, sIceDem.lat);
        [lonSubMesh,latSubMesh] = meshgrid(lonSub, latSub);
        iceSubExist = interp2(lonMesh,latMesh, tempIceBl, lonSubMesh,latSubMesh,'nearest');
    end

    sCryo.icx = zeros(numel(sHydro.lat),numel(sHydro.lon));
    for ll = 1 : numel(sHydro.lat)
        for mm = 1 : numel(sHydro.lon)
            sCryo.icx(ll,mm) = sum2d(iceSubExist(latSubInd{ll},lonSubInd{mm}))/100;
        end
    end

    %Note all main grid cells that contain any ice
    sCryo.icbl = sparse(szDem(1),szDem(2));
    sCryo.icbl(sCryo.icx > 0) = 1;

    %Create sparse ice grid (including frame) and delete previous structures:                
    %Find indices of frame around existing ice (allows glacier morphology)

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

    sCryo.iceDem =  sparse(szIce(1), szIce(2));
    sCryo.iceDem(indIceFr) = sIceDem.data(indIceFr);
    sCryo.iceLat = sIceDem.lat;
    sCryo.iceLon = sIceDem.lon;
        clear('sIceInput','sGlaciers');
        clear('sIceDem');


    %Calculate ice surface slope and flow direction (based on surface elevation)
    [~, sCryo.iceFdr, iceSurfSlope] = ...
        wflowacc(sCryo.iceLon, sCryo.iceLat, sCryo.iceDem,'type','single','edges','open','coord','geographic');

    %Calculate slope in the flow direction:
    sCryo.iceSlopeFdr = -sparse(reshape(sum(sCryo.iceFdr.*iceSurfSlope, 2), szIce));
        sCryo.iceSlopeFdr(sCryo.iceDem == 0) = 0;
    %Calculate slope in all directions:
    sCryo.iceSlope = sparse(iceSurfSlope);
        sCryo.iceSlope(sCryo.iceDem == 0, :) = 0;
        sCryo.iceSlope(:, sCryo.iceDem == 0) = 0;


%             sCryo.iceSlopeA = atand(sCryo.iceSlope);

    %Estimate thickness (force balance based on equilibrium
    %assumption from Shea et al. 2015)
    if full(all2d(sCryo.iceWE == 1 | sCryo.iceWE == 0)) %Only do if ice input is boolean field (rather than existing %thickness estimate):
        [sCryo.iceWE, iceVert] = glacier_WE_init(sCryo.iceWE, sCryo.iceSlopeFdr, sMeta);
    end

    %Estimate elevation of ice base (from dem and thickness)
    sCryo.iceBase = sCryo.iceDem - iceVert;

    %Estimate ice FDR, based on base DEM?

    %If debris layer is used, load:
    if isfield(sPath, 'debris')
        sDebris = read_geodata(sPath.debris, sMeta, 'none','no_disp');
            sDebris = squeeze(sDebris.data);
        if ~isequal(size(sDebris.data), szIce)
            error('CCHF_backbone:mismatchIceDebris', ...
                'The ice and debris cover grids are not the same size.');
        end

        %Check that coordinates line up with iceDEM:
        if sum(ismember(sDebris.lon,sCryo.iceLon)) < 0.9*numel(sDebris.lon) || sum(ismember(sDebris.lat,sCryo.iceLat)) < 0.9*numel(sDebris.lat)
           error('CCHF_backbone:misalignedIceDem','The ice DEM and ice presence grid do not appear to align. Maybe rounding error?');
        end
        sCryo.icdbr = sparse(szIce(1),szIce(2));
        sCryo.icdbr(indIce) = sDebris.data(indIce);
            clear('sDebris');
            
        warning('ice_grid_init:debrisCover',['Check that debris cover '...
            'is actually used in the ice processs representations '...
            '(including heat calculation).'])
    end


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
    
    sCryo.main2ice = cell(szDem);
    sCryo.ice2main = sparse(szIce(1),szIce(2));

    for ll = 1 : numel(sHydro.dem)
        [rC, cC] = ind2sub(szDem, ll);

        rMatch = find(sCryo.iceLat < edgeLat(rC) & sCryo.iceLat > edgeLat(rC+1));
        cMatch = find(sCryo.iceLon > edgeLon(cC) & sCryo.iceLon < edgeLon(cC+1));

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
           indCurr = indCurr(sCryo.iceDem(indCurr) ~= 0);

           %If glaciers (or frame around glaciers) exist in
           %current grid cell, record indices
           if ~isempty(indCurr)
               %Sparse array in which the value indicates which main
               %grid cell the centroid of the ice grid cell is within
               sCryo.ice2main(indCurr) = ll;

               %Cell array (inherently sparse) that has same
               %resolution as main grid and contains a list of any
               %glacier cells within that grid:
               sCryo.main2ice{ll} = indCurr;
               icBlMain(ll) = 1;
           end
        end
    end

    %Initialize array to track changes in ice depth
    sCryo.icdwe = sparse(szDem(1),szDem(2));
else %In this case, ice grid same as main grid
    sCryo.iceWE  = sparse(szDem(1),szDem(2));
    sCryo.icdwe = sparse(szDem(1),szDem(2));
    sCryo.icx = zeros(szDem,'single');
end