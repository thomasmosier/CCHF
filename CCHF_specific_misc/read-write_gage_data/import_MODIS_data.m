function sGageCurr = import_MODIS_data(path, lon, lat, mask, varargin)

%Fractional threshold for which data are flagged because too many values
%are nan (i.e. no decision):
noDecThresh = 0.6;

if ~isempty(varargin(:)) && ~isempty(varargin{1}) && regexpbl(varargin{1}, {'interp','intrp','dist', 'avg','mean','average'})
    strType = varargin{1}; 
else
	%Require user input of average or spatially distributed MODIS SCA
    %measurements:
    uiwait(msgbox(sprintf(['Enter a character string to indicate whether '...
        'to average all MODIS observations over the domain (' char(39) ...
        'avg' char(39) ') or interpolated to model domain (' char(39) ...
        'interp' char(39) ').\n']), '(Click OK to Proceed)','modal'));
    strType = input(['Load ' char(39) 'avg' char(39) ' or ' char(39) 'interp' char(39) ' MODIS images?' char(10)],'s');
end



[foldNm, fileNm, ext] = fileparts(path);

%Find all MODIS data files:
indPer = regexpi(fileNm,'\.');
filesMODIS = dir([foldNm, filesep, fileNm(1:indPer) '*', ext]);
    filesMODIS = extractfield(filesMODIS, 'name');
datesMODIS = cell(2,1); %First cell is start date and seond is end date
[datesMODIS{:}] = deal(nan(numel(filesMODIS(:)), 3));

for jj = 1 : numel(filesMODIS(:)) 
    indPer = regexpi(filesMODIS{jj},'\.');
    if numel(indPer) == 1 %Some scripts I have written to interpolate MODIS replace the periods with underscores...
        indPer = regexpi(filesMODIS{jj},'_');
    end

    datesMODIS{1}(jj,:) = ...
        days_2_date(str2double(filesMODIS{jj}(indPer(2)-3:indPer(2)-1)), ...
        [str2double(filesMODIS{jj}(indPer(2)-7:indPer(2)-4)), 1, 1], 'gregorian');
end

%Check that MODIS data type is consistent with expectation:
if regexpbl(filesMODIS{1}, 'MOD10CM') %Monthly, geogeraphic coordinates
    datesMODIS{2}(:,1:2) = datesMODIS{1}(:,1:2);
    datesMODIS{2}(:,3) = eomday(datesMODIS{1}(:,1), datesMODIS{1}(:,2));
    varMod = 'MOD_CMG_Snow_5km';
    fieldMod = 'Snow_Cover_Monthly_CMG';
    fieldCloud = '';
    evalType = 'max';
elseif regexpbl(filesMODIS{1}, 'MOD10C1') %Daily, geographic coordinates
    datesMODIS{2} = datesMODIS{1};
    varMod = 'MOD_CMG_Snow_5km';
    fieldMod = 'Day_CMG_Snow_Cover';
    fieldCloud = 'Day_CMG_Cloud_Obscured';
    evalType = 'max';
    %'Day_CMG_Confidence_Index'
elseif regexpbl(filesMODIS{1}, 'MOD10C2') %8-Day, geographic coordinates
    dateRef = [2000,1,1];
    datesMODIS{2} = days_2_date(days_since(dateRef,datesMODIS{1}, 'gregorian') + 7, dateRef, 'gregorian');
    varMod = 'MOD_CMG_Snow_5km';
    fieldMod = 'Eight_Day_CMG_Snow_Cover';
    fieldCloud = 'Eight_Day_CMG_Cloud_Obscured'; %Must scale by cloud
    evalType = 'max';
    %'Eight_Day_CMG_Confidence_Index'
elseif regexpbl(filesMODIS{1}, 'MOD10A2')
    dateRef = [2000,1,1];
    datesMODIS{2} = days_2_date(days_since(dateRef,datesMODIS{1}, 'gregorian') + 7, dateRef, 'gregorian');
    varMod = 'MOD_Grid_Snow_500m';
    fieldMod = 'Maximum_Snow_Extent'; 
    evalType = 'max';
else
    %Use import wizard to figure out field names for different HDF
    %datasets
    indPer = regexpi(filesMODIS{1},'\.');
    warning('read_gagedata:MODIStype',['The selected MODIS data '...
        'are type ' char(39) filesMODIS{1}(1:indPer(1)-1) ...
        char(39) ', which is currently not supported. '...
        'Therefore the data will not be loaded. Consider adding this dataset.']);
    return
end

%Determine file type:
if regexpbl(filesMODIS{1},'hdf')
    strFile = 'hdf';
elseif regexpbl(filesMODIS{1}, {'nc','netCDF'})
    strFile = 'nc';
else
    error('import_MODIS_data:unknownFileType', [filesMODIS{1} ' is of an unknown filetype.'])
end




%         %Keep only temporal subset of files:
%         if ~isempty(time)
%             daysMODIS = days_since([1800, 1, 1], datesMODIS, 'gregorian');
%             daysKeep  = days_since([1800, 1, 1],       time, 'gregorian');
%                 daysKeep = sort(daysKeep);
%             indKeep = find(daysMODIS >= daysKeep(1) & daysMODIS <= daysKeep(2));
%             
%             filesMODIS = filesMODIS(indKeep);
%             datesMODIS = datesMODIS(indKeep,:);
%         end

if ~isempty(mask)
    szMODIS = [numel(datesMODIS{1}(:,1)), size(mask)];
    dataMODIS = nan(szMODIS);
    
    indNanMask = find(isnan(mask));
    if ~isempty(indNanMask)
       [rNan, cNan] = ind2sub(szMODIS(2:3), indNanMask); 
    end
else
    warning('read_gagedata:maskNeeded',['Loading MODIS data '...
        'requires a mask as an optional input argument. The ' ...
        'mask must be the sptial grid of the model and nan at '...
        'locations where the model will not be evaluated.']);
    
    indNanMask = [];
end



fill = zeros(numel(datesMODIS{1}(:,1)),1,'single');
% if regexpbl(strType, {'avg','mean','average'})
%     fill = zeros(numel(datesMODIS{1}(:,1)),1,'single');
% elseif regexpbl(strType, {'interp','intrp','dist'})
%     fill = zeros(szMODIS,'single');
% end

% [lonInterp, latInterp] = meshgrid(lon,lat);
lonInterp = lon;
latInterp = lat;

%Loop over all MODIS files to load:
for jj = 1 : numel(datesMODIS{1}(:,1))
    if regexpbl(strFile, 'hdf')
        if jj == 1 %Calculate MODIS GRID
            lonEdgMODIS = linspace(-180,180,7201);
            latEdgMODIS = linspace(90,-90,3601);

            nBord = 5; %add border around MODIS (necessary for proper interpolation:
            indLon = [find(lonEdgMODIS >= lon(1),1,'first')-1 - nBord, find(lonEdgMODIS <= lon(end),1,'last') + nBord];
            indLat = [find(latEdgMODIS <= lat(1),1,'first')-1 - nBord, find(latEdgMODIS >= lat(end), 1,'last') + nBord];
            indLon(indLon < 1) = 1;
            indLon(indLon > 7200) = 7200;
            indLat(indLat < 1) = 1;
            indLat(indLat > 3600) = 3600;

            dxy = 0.05;
            lonMODIS = linspace(-180+dxy/2, 180-dxy/2, 7200);
            latMODIS = linspace(  90-dxy/2, -90+dxy/2, 3600);

            [lonMODIS, latMODIS] = meshgrid(lonMODIS(indLon(1):indLon(2)), latMODIS(indLat(1):indLat(2)));
        end

        %Read each MODIS file 
    %             MODISTemp3 = single(hdfread(fullfile(foldNm,filesMODIS{jj}), varMod, 'Fields', 'Eight_Day_CMG_Confidence_Index', 'Index',{[indLat(1) indLon(1)],[1 1],[nLat nLon]})); 
    %             MODISTemp2 = single(hdfread(fullfile(foldNm,filesMODIS{jj}), varMod, 'Fields', 'Eight_Day_CMG_Cloud_Obscured', 'Index',{[indLat(1) indLon(1)],[1 1],[nLat nLon]})); 
        nLon = indLon(2) - indLon(1) + 1;
        nLat = indLat(2) - indLat(1) + 1; 
        MODISTemp = single(hdfread(fullfile(foldNm,filesMODIS{jj}), varMod, ...
            'Fields', fieldMod, 'Index',{[indLat(1) indLon(1)],[1 1],[nLat nLon]})); 

        MODISTemp(MODISTemp > 100 | MODISTemp < 0) = nan;

        %If cloud field (is case for all images with temporal
        %resolution less than 1 month):
        if ~isempty(fieldCloud)
            MODISCld = single(hdfread(fullfile(foldNm,filesMODIS{jj}), varMod, ...
                'Fields', fieldCloud, 'Index',{[indLat(1) indLon(1)],[1 1],[nLat nLon]}));  

           MODISTemp = MODISTemp ./ (100 - MODISCld);
           MODISTemp(MODISCld >= 0.3) = NaN;
        end
    elseif regexpbl(strFile,'nc')
        fileC = fullfile(foldNm,filesMODIS{jj});
        
%         sDataCurr.time = ncread(fileC,'time');
%         sDataCurr.attTime = ncinfo(fileC,'time');
%             sDataCurr.attTime = squeeze(struct2cell(sDataCurr.attTime.Attributes))';
%         indUnit = strcmpi(sDataCurr.attData,'units');
%             [rUnit, cUnit] = find(indUnit == 1);
%             dataUnit = sDataCurr.attData{rUnit, cUnit+1};

        fInfo = ncinfo(fileC);
            varAvail = extractfield(fInfo, 'Variables');
            if iscell(varAvail) && isstruct(varAvail{1})
               varAvail = extractfield(varAvail{1},'Name'); 
            end
            
            varAvail(strcmpi(varAvail,'lon')) = [];
            varAvail(strcmpi(varAvail,'lat')) = [];
            varAvail(strcmpi(varAvail,'time')) = [];
            varAvail(strcmpi(varAvail,'time_bnds')) = [];
        if numel(varAvail(:)) == 1
            varLd = char(varAvail{1});
        else
            error('import_MODIS_data:multipleDataFields', ['There appear '...
                'to be multiple data fields in the current NetCDF file.']);
        end
        
        lonMODIS = ncread(fileC, 'lon');
            lonMODIS = lonMODIS(:)';
        latMODIS = ncread(fileC, 'lat');
            latMODIS = latMODIS(:);
        MODISTemp = ncread(fileC, varLd);
        
        attData = ncinfo(fileC, varLd);
            attData = squeeze(struct2cell(attData.Attributes))';
        [rNoVal, ~] = ind2sub(size(attData),find(strcmpi(attData, 'missing_value') == 1));
        valNoData = attData{rNoVal,2};
        MODISTemp(MODISTemp == valNoData) = nan;
    else
        error('import_MODIS_data:unknownFileType',['The MODIS file type is ' ...
            strType ', which is an unknown type.']);
    end
    
    %MODIS KEY: 
    %0-100=percent of snow in cell, 
    %211=night, 250=cloud, 253=no decision, 254=water mask, 255=fill

    if numel(size(MODISTemp)) == 3 
        if numel(MODISTemp(:,1,1)) == 1
            MODISTemp = squeeze(MODISTemp);
        else
            error('import_MODIS_data:temporalArrayLoaded', ...
                ['The loaded MODIS grid has a non-singleton temporal '...
                'dimension, which is not expected (i.e. has not been '...
                'programmed for.']);
        end
    end
     
    %Interpolate MODIS data to model grid:
    if ~isequal(round2(lon,3), round2(lonMODIS,3)) || ~isequal(round2(lat,3), round2(latMODIS,3)) 
        warning('off','all');
        dataMODIS(jj,:,:) = round(PCHIP_2D(lonMODIS, latMODIS, single(MODISTemp), lonInterp, latInterp));
    %         dataMODIS(jj,:,:) = round(interp2(lonMODIS,latMODIS,single(MODISTemp),lonInterp, latInterp,'nearest'));
        warning('on','all');
    else
        dataMODIS(jj,:,:) = single(MODISTemp);
    end
    
    %Set interpolated data to nan if mask present:
    if ~isempty(indNanMask)
       indNan3 = jj + (rNan-1)*szMODIS(1) + (cNan-1)*szMODIS(2)*szMODIS(1);
       dataMODIS(indNan3) = nan;
    end
    
    %Flag timesteps where more than given percentage of data are nan
    if sum2d(isnan(MODISTemp)) >= noDecThresh*numel(MODISTemp)
        fill(jj) = 1;
    end
%     if regexpbl(strType,{'avg','mean','average'})
%         if sum2d(isnan(MODISTemp)) >= noDecThresh*numel(MODISTemp)
%             fill(jj) = 1;
%         end
%     elseif regexpbl(strType,{'interp','intrp','dist'})
%         [rFl, cFl] = ind2sub(szMODIS(2:3), find(isnan(squeeze(dataMODIS(jj,:,:))))); 
%         indFl3 = jj + (rFl-1)*szMODIS(1) + (cFl-1)*szMODIS(2)*szMODIS(1);
%         fill(indFl3) = 1;
%     end
end


%Create output from MODIS time-series (either averaged or distributed):
if regexpbl(strType,{'avg','mean','average'})
    dataGage = nan(numel(datesMODIS{1}(:,1)),1);
    for jj = 1 : numel(datesMODIS{1}(:,1))
        dataGage(jj) = round(mean2d(dataMODIS(jj,:,:)));
    end

     sGageCurr = struct(...
        'lon',  'avg',...
        'lat',  'avg',...
        'casi', single(dataGage),... %Stands for Covered Area Snow and Ice
        'flag', fill, ...
        'dateStart', datesMODIS{1}, ...
        'dateEnd', datesMODIS{2}, ...
        '', evalType);
elseif regexpbl(strType,{'interp','intrp','dist'})
     sGageCurr = struct(...
        'lon',  'all',...
        'lat',  'all',...
        'casi', single(dataMODIS),... %Stands for Covered Area Snow and Ice
        'path_casi', path,...
        'flag', fill, ...
        'dateStart', datesMODIS{1}, ...
        'dateEnd', datesMODIS{2}, ...
        'type_casi', evalType); 
else
   error('import_MODIS_data:unknownInterpretationType', [strType ...
       ' is an unknown MODIS data interpretation type.']);
end