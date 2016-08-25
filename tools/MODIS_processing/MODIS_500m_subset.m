%Folder with MODIS data:
foldMain = '/Users/thomas/Desktop/research/geodata/MODIS/MODIS-8day-500m';
%Path to reference grid (only lat-lon matters):
pathRefGrid = '/Users/thomas/Desktop/research/hydro_model_cases/Wolverine/wolverine_dem_wc.asc';
nmReg = 'wolverine';
% pathRefGrid = '/Users/thomas/Desktop/research/hydro_model_cases/Gulkana/gulkana_dem_wc.asc';
% nmReg = 'gulkana';


outType = 1;    %1 = output uses same key as MODIS input and are resampled using nearest neighbor interpolation. 
                %~= 1, interpolates MODIS using bilinear interpolation and reports fractional snow/ice cover

if outType == 1
    interpMeth = 'nearest';
    cldThresh = 1; 
    
    nmInterp = [nmReg '_org_values'];
else
    interpMeth = 'linear';
    cldThresh = 0.6;
    
    nmInterp = [nmReg '_' interpMeth 'Interp'];
end

crd = [-155, -145, 59, 65]; %[W, E, S, N];

%If dataset is 8-day 500 m:
varMod = 'MOD_Grid_Snow_500m';
fieldMod = 'Maximum_Snow_Extent'; %'Maximum_Snow_Extent', 'Eight_Day_Snow_Cover'
%MODIS Key:
valNoDec = [0, 1, 11, 50, 254, 255]; %No decision values corresponding to 
    %'mising data', 'no decision', 'night', 'cloud', 'detector saturated', 'fill'
valOthCov = [37, 39]; %Other land cover values correspinding to 'lake', 'ocean'
valSnow = [100, 200]; %'lake ice', 'snow'
valNSnow = 25; %'no snow'

%create output folder:
foldOut = fullfile(foldMain, 'intrp2geo', nmInterp);

if ~exist(foldOut, 'dir')
    mkdir(foldOut);
end

%Load reference file:
if exist(pathRefGrid, 'file')
    sMeta = struct;
    sRefGrid = read_geodata(pathRefGrid, sMeta, 'no_disp');

    [lonTemp, latTemp] = meshgrid(sRefGrid.lon, sRefGrid.lat);
    
    lonOut = reshape(lonTemp, [], 1);
    latOut = reshape(latTemp, [], 1);
    
    lonEdg = box_edg(sRefGrid.lon);
    latEdg = box_edg(sRefGrid.lat);
    crd = [lonEdg(1), lonEdg(end), latEdg(end), latEdg(1)];
    disp('Output grid loaded from reference file');
else
    disp('Output grid loaded from vector in script preamble.');
end

%Remove folders that do not correspond to MODIS tiles:
foldTiles = dir(foldMain);
foldTiles = foldTiles(vertcat(foldTiles.isdir));
foldTiles = extractfield(foldTiles, 'name');
for ii = numel(foldTiles) : -1 : 1
    if ~regexpbl(foldTiles{ii},{'h','v'},'and')
        foldTiles(ii) = [];
    end
end

%%Load one test file from each tile:
% nLon = indLon(2) - indLon(1) + 1;
% nLat = indLat(2) - indLat(1) + 1; 




%Load one MODIS file from each tile in order to determine which tiles must
%be loaded:
for ii = numel(foldTiles) : -1 : 1
    fileCurr = dir(fullfile(foldMain, foldTiles{ii}, '*.hdf'));
    fileCurr = extractfield(fileCurr(1), 'name');
    
    pathTileTest = fullfile(foldMain, foldTiles{ii}, char(fileCurr));
   
    
    %Determine coordinates of HDF file:
    geoInfo = hdfinfo(pathTileTest, 'eos');    % HDF file info in EOS mode
        UL = geoInfo.Grid(1).UpperLeft;     % UpperLeft corner (lon, lat)
        LR = geoInfo.Grid(1).LowerRight;    % LowerRight corner (lon, lat)
    
    [lonTemp, latTemp] = MODIS_proj2geo([UL(1); LR(1)], [UL(2); LR(2)]);
        ULgeo = [lonTemp(1), latTemp(1)];
        LRgeo = [lonTemp(2), latTemp(2)];
    
    %Remove tile if entirely outside specified bounding box:
    if ~in_bounds(ULgeo, LRgeo, crd)
        foldTiles(ii) = [];
    end
    
%     Eight_Day_Snow_Cover = hdfread(filename, ...
%         '/MOD_Grid_Snow_500m/Data Fields/XDim', 'Index', {[1  1],[1  1],[2400  2400]});
%     tileTest{ii} = int(hdfread(pathTileTest{ii}, varMod, ...
%         'Fields', fieldMod)); 
%     single(hdfread(pathTileTest{ii}, varMod, ...
%         'Fields', fieldMod, 'Index',{[indLat(1) indLon(1)],[1 1], [nLat nLon]})); 
end


%%FIND GEOGRAPHIC COORDINATES OF EACH TILE WITHIN REGION OF INTEREST:
fileTiles = cell(numel(foldTiles), 1);
lonGeo = cell(numel(foldTiles), 1);
latGeo = cell(numel(foldTiles), 1);
indLon = cell(numel(foldTiles), 1);
indLat = cell(numel(foldTiles), 1);
tileDates = cell(numel(foldTiles), 1);
nTilePts = cell(numel(foldTiles), 1);
nPts = 0;

for ii = 1 : numel(foldTiles)
    fileCurr = dir(fullfile(foldMain, foldTiles{ii}, '*.hdf'));
    fileTiles{ii} = extractfield(fileCurr, 'name');
    pathTileTest = fullfile(foldMain, foldTiles{ii}, fileTiles{ii}{1});
   
    if ii == 1
        indPer = regexpi(fileTiles{ii}{1}, '\.');
        modVers = fileTiles{ii}{1}(1:indPer(1)-1); 
    end
    
    
    %Determine coordinates of HDF file:
    geoInfo = hdfinfo(pathTileTest,'eos');    % HDF file info in EOS mode
        nCol = geoInfo.Grid(1).Columns;	% columns of data
        nRow = geoInfo.Grid(1).Rows;   	% rows of data
        UL = geoInfo.Grid(1).UpperLeft;     % UpperLeft corner (lon, lat)
        LR = geoInfo.Grid(1).LowerRight;    % LowerRight corner (lon, lat)
        xRes=(LR(1)-UL(1))/nCol;
        yRes=(LR(2)-UL(2))/nRow;
    
    xVec = linspace(UL(1) + 0.5*xRes, LR(1) - 0.5*xRes, nCol);   % columns of data
    yVec = linspace(UL(2) - 0.5*yRes, LR(2) + 0.5*yRes, nRow);   % rows of data

    [xGrid, yGrid] = meshgrid(xVec, yVec);
    xGrid = reshape(xGrid, [], 1);
    yGrid = reshape(yGrid, [], 1);
    [lonGeoTemp, latGeoTemp] = MODIS_proj2geo(xGrid, yGrid);
    lonGeo{ii} = reshape(lonGeoTemp, nRow, nCol);
    latGeo{ii} = reshape(latGeoTemp, nRow, nCol);
        lonMin = min(lonGeo{ii},[],1);
        lonMax = max(lonGeo{ii},[],1);
        latMin = min(latGeo{ii},[],2);
        latMax = max(latGeo{ii},[],2);
        
    
    %Determine the indices of the lat grid within the region of interest:
    %Note that the lon in a single column will vary slightly because
    %original data are in projected coordinates.
    indW = find(lonMax < crd(1), 1, 'last');
    if isempty(indW)
        indW = 1;
    end
    indE = find(lonMin > crd(2), 1, 'first');
    if isempty(indE)
        indE = numel(lonGeo{ii}(1,:));
    end
    
    indN = find(latMin > crd(4), 1, 'last');
    if isempty(indN)
        indN = 1; 
    end
    indS = find(latMax < crd(3), 1, 'first');
    if isempty(indS)
        indS = numel(latGeo{ii}(:,1));
    end
    
    indLon{ii} = [indW, indE];
    indLat{ii} = [indN, indS];
    
    nTilePts{ii} = numel(latGeo{ii}(indLat{ii}(1):indLat{ii}(2),indLon{ii}(1):indLon{ii}(2)));
    nPts = nPts + nTilePts{ii};
    %Find available dates within folders:
    if ii == 1
        nTsMax = numel(fileTiles{ii});
        nTsMin = numel(fileTiles{ii});
    else
        nTsMax = max(nTsMax, numel(fileTiles{ii}));
        nTsMin = min(nTsMin, numel(fileTiles{ii}));
    end
    
    tileDates{ii} = nan(numel(fileTiles{ii}), 2);
    
    for jj = 1 : numel(fileTiles{ii})
        tileDates{ii}(jj,:) = MODIS_date(fileTiles{ii}{jj});
    end
end


%Check if any dates missing:
if nTsMax ~= nTsMin
    error('MODIS_500m_subset:missingTs',['The number of time-series '...
        'files in the directories is not uniform. Write code to '...
        'determine which they are.']);
end


%%RESAMPLE MODIS DATA TO STANDARD GEOGRAPHIC GRID (as in reference file)
%The approach is to do the analysis by time-series element
for ii = 1 : numel(tileDates{1}(:,1))
    %Initialize array for compiled tiles:
    dataIn = nan(nPts,1);
    %Initialize lat and lon vecs:
    if ii == 1
        latIn = nan(nPts,1);
        lonIn = nan(nPts,1);
    end
    
    %Now, loop over tiles (loading data for each)
    indCurr = 0;
    for jj = 1 : numel(foldTiles)
        %Check that current file corresponds to correct timestep
        dateCurr = MODIS_date(fileTiles{jj}{ii});
        
        if ~isequal(dateCurr, tileDates{1}(ii,:))
            error('MODIS_500m_subset:currDateWrong',['The date for the '...
                'MODIS data about to be loaded does not match the '...
                'expected date (based on loop over time series).']);
        end

        currDataIn = hdfread(fullfile(foldMain, foldTiles{jj}, fileTiles{jj}{ii}), varMod, ...
            'Fields', fieldMod);
        
%         currDataIn = single(hdfread(fullfile(foldMain, foldTiles{jj}, fileTiles{jj}{ii}), varMod, ...
%             'Fields', fieldMod, 'Index', {[indLat{jj}(1) indLon{jj}(1)],[1 1], [indLat{jj}(2) indLon{jj}(2)]})); 
        
        dataIn(indCurr + 1 : indCurr + nTilePts{jj}) = ...
            reshape(currDataIn(indLat{jj}(1):indLat{jj}(2), indLon{jj}(1):indLon{jj}(2)),[],1);

        %Record lat and lon grid in correct format:
        if ii == 1
            latIn(indCurr + 1 : indCurr + nTilePts{jj}) = ...
                reshape(latGeo{jj}(indLat{jj}(1):indLat{jj}(2), indLon{jj}(1):indLon{jj}(2)),[],1);
            lonIn(indCurr + 1 : indCurr + nTilePts{jj}) = ...
                reshape(lonGeo{jj}(indLat{jj}(1):indLat{jj}(2), indLon{jj}(1):indLon{jj}(2)),[],1);
        end
        
        indCurr = indCurr + nTilePts{jj};
    end
    
    if outType == 0
        %Transform MODIS values (1, 0, nan):
        %Other Labd cover values:
        for kk = 1 : numel(valOthCov)
            dataIn(dataIn == valOthCov(kk)) = nan;
        end
        %No decision values:
        noDec = 10000;
        for kk = 1 : numel(valNoDec)
            dataIn(dataIn == valNoDec(kk)) = 10000;
        end
        %Snow values:
        for kk = 1 : numel(valSnow)
            dataIn(dataIn ==  valSnow(kk)) = 1;
        end
        for kk = 1 : numel(valNSnow)
            dataIn(dataIn == valNSnow(kk)) = 0;
        end

        %Set all values to no data if pixels with no decision greater than 
        %threshold fraction
        dataCldTest = griddata(lonIn, latIn, dataIn, lonOut, latOut, 'nearest');
        
        if sum2d(dataCldTest == noDec) > cldThresh*numel(dataCldTest)
            dataIn(:) = nan;
        end
    end
    
    %Interpolate compiled tiles to output grid:
    if outType == 1
        %Try to preserve original values
        dataOut = griddata(lonIn, latIn, dataIn, lonOut, latOut, interpMeth);
    elseif outType == 0
        %'dataOut' has units of percentage snow cover
        dataOut = 100*griddata(lonIn, latIn, dataIn, lonOut, latOut, interpMeth);
    	dataOut(dataOut > 100) = NaN;
    end
    dataOut = reshape(dataOut, numel(sRefGrid.lat), numel(sRefGrid.lon));
    
    %Write output as NetCDF:
    if tileDates{1}(ii,2) < 10 
        strTime = [num2str(tileDates{1}(ii,1)) '00' num2str(tileDates{1}(ii,2))];
    elseif tileDates{1}(ii,2) < 100 
        strTime = [num2str(tileDates{1}(ii,1)) '0' num2str(tileDates{1}(ii,2))];
    else
        strTime = [num2str(tileDates{1}(ii,1)) num2str(tileDates{1}(ii,2))];
    end
    pathCurr = fullfile(foldOut, [modVers '.A' strTime '.' nmReg '.' interpMeth '.nc']);
    dateVec = days_2_date(tileDates{1}(ii,2)-1, [tileDates{1}(ii,1), 1, 1], 'gregorian');
    
    print_grid_NC(pathCurr, dataOut, 'sca', sRefGrid.lon, sRefGrid.lat, dateVec, -1);
end


%DATES MODIS:
% for jj = 1 : numel(filesMODIS(:)) 
%     indPer = regexpi(filesMODIS{jj},'\.');
% 
%     datesMODIS{1}(jj,:) = ...
%         days_2_date(str2double(filesMODIS{jj}(indPer(2)-3:indPer(2)-1)), ...
%         [str2double(filesMODIS{jj}(indPer(2)-7:indPer(2)-4)), 1, 1], 'gregorian');
% end
%     dateRef = [2000,1,1];
%     datesMODIS{2} = days_2_date(days_since(dateRef,datesMODIS{1}, 'gregorian') + 7, dateRef, 'gregorian');
