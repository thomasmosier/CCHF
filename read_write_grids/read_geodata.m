function sData = read_geodata(pathData, sMeta, varargin)
% Copyright 2013, 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the GlobalClimateData Downscaling Package.
% 
% The GlobalClimateData Downscaling Package is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The Downscaling Package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%Read NetCDF or ESRI ASCII geo-referenced data and write output
%INPUT:
%'path' refers to a geo-referenced NetCDF or ESRI formatted ASCII file
%'sMeta' is a structure array with the following fields (among others):
    %sMeta.currVar = variable to load
    %sMeta.yrsOut = [first year of output data, last year of output data]
    %sMeta.yrsClim = [first year of climatology, last year of climatology]
    %sMeta.mnths = list of months requested
    %sMeta.crd = [lonW, lonE, latS, latN]
    %sMeta.currTime = [year, month] (optional: [year, month, day])

%'varargin' used to specify if a climatology is going to be produced.  
%If varargin contains 'c', the years loaded are 'sMeta.yrsClim' and the
%resulting data are averaged.  Otherwise, data loaded correspond to 
%'sMeta.yrsOut'. 

%Output: 
%'sData' isstructure array with format:
    %sData.info = Cell array of general dataset information 
    %sData.lat = vector of latitudes
    %sData.attLat = cell array with attribute list for latitude coordinate
    %sData.lon = vector of longitudes
    %sData.attLon = cell array with attribute list for longitude coordinate
    %sData.time = vector of times corresponding to each data grid loaded
    %sData.attTime = cell array with attribute list for time coordinate
    %sData.data = 2 or 3 dimensional data array with climate data
    %sData.attData = cell array with attribute list for data
    %sData.var = Name of the climate variable in the structure
    
    

%Choose which time-elemnts to load based upon 'varargin' (default = 'ts')
if ~isempty(varargin(:)) 
    if regexpbl(varargin{1},'disp') && numel(varargin(:)) > 1
        [outTyp, yrsLd, mnthsLd] = output_type(varargin{2}, sMeta);
    elseif ~regexpbl(varargin{1},'disp')
        [outTyp, yrsLd, mnthsLd] = output_type(varargin{1}, sMeta);
    else
        [outTyp, yrsLd, mnthsLd] = output_type('', sMeta);
    end
else
    [outTyp, yrsLd, mnthsLd] = output_type('', sMeta);
end


if ~regexpbl(varargin,'no_disp')
    if isfield(sMeta,'currTime')
        disp([mnth_str(sMeta.currTime(2)) ' data are now being read from ' char(39) ...
            pathData char(39) '.']);
    else
        disp(['Data are now being read from ' char(39) ...
            pathData char(39) '.']);
    end
end

if ~isnumeric(outTyp) && regexpbl(outTyp,'all')
    warning('read_geodata:outTypAll',['The case for ' outTyp ...
        ' has not been coded yet.  Once coded, it will load all months '...
        'and all years specified by metadata structure']);
end

%Initialize time vector:
% if ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'none'))
%     aryTime = [];
% elseif isnumeric(outTyp)
%     aryTime = nan(numel(outTyp(:,1)),3);
% else
%     aryTime = nan((yrsLd(2)-yrsLd(1)+1)*length(mnthsLd),3);
% end
climTime = [];
    
%Initialize Structure:
sData = struct('lat', NaN, 'attLat', cell(1,1), 'lon', NaN, ...
    'attLon', cell(1,1), 'time', NaN, 'attTime', cell(1,1), ...
    'data', NaN, 'attData', cell(1,1));


%Find type of file being loaded:
if isdir(pathData)
    %Determine the extension of the files in the path directory:
    filesTemp = dir(pathData);
    for ii = 1 : length(filesTemp)
        [~,~,ext] = fileparts(filesTemp(ii).name);
        if ~isempty(regexpi(ext,'asc'))
            break
        end
    end
    
else
    [~,~, ext] = fileparts(pathData);
    filesTemp = cell(1,1);
    filesTemp{1} = pathData;
end

%Convert filesTemp to cell array
if isstruct(filesTemp)
    filesTemp = struct2cell(filesTemp);
    filesTemp = filesTemp(1,3:end);
end

if isfield(sMeta,'crd')
    %Ensure bounds are in correct order:
    if sMeta.crd(3) > sMeta.crd(4)
        warning('rego_geodata:lonBounds','Longitudinal bounds being flipped because order appears incorrect');
        sMeta.crd(3:4) = fliplr(sMeta.crd(3:4));
    end
    if sMeta.crd(1) > sMeta.crd(2)
        warning('rego_geodata:latBounds','Latitudinal bounds being flipped because order appears incorrect');
        sMeta.crd(1:2) = fliplr(sMeta.crd(1:2));
    end
end


%Initialize aryTime, which is necessary for non-temporal data:
aryTime = nan(1,2);

% flagDay = 0; %Assume monthly data until proven otherwise

%%IDENTIFY TYPE OF DATA:

%Load NetCDF file format:
if regexpbl(ext,{'nc'}) && ~regexpbl([filesTemp{1},pathData],'gpcc') %% 2nd clause necessary because GPCC uses odd NetCDF format
    %Find file to load:
    if isdir(pathData)
        fileNcTemp = dir(fullfile(pathData,'*.nc'));
            fileNcTemp = extractfield(fileNcTemp, 'name');
            
        if isempty(fileNcTemp)
            error('read_geodata:noNcFiles1',['No NetCDF data were found in ' pathData '.']);
        end
    elseif regexpbl(pathData,'.nc')
        fileNcTemp{1} = pathData;
        [yrsGcm, mnthsGcm, ~] = CMIP5_time(pathData);
    else
        error('read_geodata:noNcFiles3',['No NetCDF data were found in ' pathData '.']);
    end
    
    
    %Test type of NetCDf file (try to find function that can read time from name:
    testNc = CMIP5_time(fileNcTemp{1});
    if all2d(isnan(testNc))
        warning('read_geodata:unknownNcTime',['The time of the current '...
            'NetCDF file is unknown. This is due to the format of the file string.']);
    end

    %Find files to load
    if ~regexpbl(pathData,'.nc') && isnumeric(outTyp)
        warning('read_geodata:ncFileFind', ['Find current NetCDF file. '...
            'This is slow, consider providing specific file to read_geodata function.']);
        
        indTest = numel(outTyp);
        for ii = 1 : numel(fileNcTemp)
          	testDate = CMIP5_time(fileNcTemp{ii});
            if ii == 1
                testDate = testDate(:,~isnan(testDate(1,:)));
                indTest = min(testDate(1,:),indTest); 
            end
            
            if all([outTyp(1:indTest) >= testDate(1,1:indTest), outTyp(1:indTest) <= testDate(1,1:indTest)])
                pathData = fullfile(pathData,fileNcTemp{ii});
                yrsGcm(1) = -9999;
                break
            end
        end
    end
    

    %Compare time-series requested to those available in specified NetCDF
    %file (using CMIP5 file naming convention)
    gcmUse = [];
    if sum(~isnan(yrsGcm)) > 0 && yrsGcm(1) ~= -9999 %Check that the time format of the GCM is recognized (based upon file name).
        %If not all time-series elements present, search folder for analogous
        %files
        if min2d(yrsLd) < yrsGcm(1) || max2d(yrsLd) > yrsGcm(2) || any(~ismember(mnthsLd, (mnthsGcm(1):mnthsGcm(2))))
            [foldNc, fileNc, extNc] = fileparts(pathData);
            indRt = regexpi(fileNc, '_');
            if ~isempty(indRt)
                fileNc = fileNc(1:indRt(end));
            end
            filesNc = dir([foldNc, filesep, fileNc, '*', extNc]);
            filesNc = extractfield(filesNc,'name');

            %Choose which files to load and in which order (ensure data are
            %continuous)
            gcmUse = zeros(length(filesNc(:)),4);
            for ii = 1 : length(filesNc(:))
                gcmTemp = CMIP5_time(fullfile(foldNc, filesNc{ii}));
                gcmUse(ii,:) = [gcmTemp(1,1), gcmTemp(2,1), gcmTemp(1,2), gcmTemp(2,2)];
                
                %Do not load if ...
                if gcmUse(ii,2) < min2d(yrsLd) || gcmUse(ii,1) > max2d(yrsLd) %GCM ends before requested years or begins after requested years.
                    gcmUse(ii,:) = NaN;
                elseif gcmUse(ii,2) == min2d(yrsLd) && any(~ismember(mnthsLd,(gcmUse(ii,3):gcmUse(ii,4)))) %GCM ends in first year to load and does not contain requested months
                    gcmUse(ii,:) = NaN;
                elseif gcmUse(ii,1) == max2d(yrsLd) && any(~ismember(mnthsLd,(gcmUse(ii,3):gcmUse(ii,4)))) %GCM begins in last year to load and does not contain requested months
                    gcmUse(ii,:) = NaN;
                end
            end

            %Delete files that do not contain necessary time-series elements:
            filesNc(isnan(gcmUse(:,1))) = [];
            gcmUse(isnan(gcmUse(:,1)),:) = [];

            %Find temporal ordering:
            [~, ordStrtGcm] = sort(gcmUse(:,1));
            [~, ordEndGcm] = sort(gcmUse(:,2));
            if ~isequal(ordStrtGcm,ordEndGcm)
                warning('read_geodata:GcmOrder',['The temporal ordering of '...
                    'GCMs to load depends on whether the start date or end '...
                    'dates are compared.  The ordering using starts dates '...
                    'will be used.'])
            end

            %Order files based upon dates:
            filesNc = filesNc(ordStrtGcm);
            gcmUse = gcmUse(ordStrtGcm,:);

            %Check that the data appear continuous:
            if length(filesNc(:)) > 1 && ~isequal(gcmUse(1:end-1,2)+1, gcmUse(2:end,1))
                error('read_geodata:GcmOrder',['The temporal coverage of '...
                    'the GCM data appear to not be continuous.'])
            end
            
            %Loop over each of file to load
            for ii = 1 : length(filesNc(:))
                if ii == 1 && length(filesNc(:)) > 1
                    warning('read_geodata:splice',['The NetCDF data requested '...
                        'are contained within multiple files, which will be '...
                        'automatically spliced together by the program.']);
                end

                warning('off','gcm_load:startYr');
                warning('off','gcm_load:EndYr');
                warning('off','NC_monthly_time_use:notAvail');
                [sDataTemp] = NC_load_v2(fullfile(foldNc, filesNc{ii}), sMeta, outTyp);
                warning('on');

                if ii == 1
                    sData = sDataTemp;
                else
                    %Convert time to common reference:
                    [gcmRefDateTemp, ~] = NC_time_units(sDataTemp.attTime);
                    calTemp =  NC_cal(sDataTemp.attTime);
                    [gcmRefDate, ~] = NC_time_units(sData.attTime);
                    cal =  NC_cal(sDataTemp.attTime);
                    
                    dateVecTemp = days_2_date(sDataTemp.time, gcmRefDateTemp, calTemp);
                    sDataTemp.time = days_since(gcmRefDate, dateVecTemp, cal);
                    sData.time = [sData.time; sDataTemp.time];
                    sData.data = cat(1, sData.data, sDataTemp.data);
                end
                
                %check that all data loaded:
                if ii == length(filesNc(:))
                    [gcmRefDate, gcmUnits] = NC_time_units(sData.attTime);
                    cal =  NC_cal(sData.attTime);
                    if regexpbl(gcmUnits,'days')
                        dateGCMVec = days_2_date(sData.time, gcmRefDate, cal);
                    elseif regexpbl(gcmUnits,'hours')
                        dateGCMVec = days_2_date(sData.time/24, gcmRefDate, cal);
                    else
                        error('read_geodata:unknownUnits',['The units '...
                            gcmUnits ' have not been programmed for.']);
                    end
                    
                    if ~all2d(isnan(dateGCMVec))
                        if isempty(outTyp)
                            if numel(filesTemp) == 1
                                aryTime = dateGCMVec;
                            else
                                error('read_geodata:ambigFiles',...
                                    ['Multiple files are in the '...
                                    'currently selected directory but '...
                                    'no information has been provided '...
                                    'on which to load.']);
                            end
                        elseif isnumeric(outTyp)
                            nYrsLd = numel(outTyp(:,1));

                            for kk = 1 : numel(mnthsLd)
                               aryTime((kk-1)*nYrsLd+1:kk*nYrsLd,:) = [yrsLd(:,1),mnthsLd(kk)*ones(nYrsLd,1)];
                            end
                        elseif all(~isnan(yrsLd)) && all(~isnan(mnthsLd))
                            nYrsLd = yrsLd(end)-yrsLd(1)+1;

                            for kk = 1 : numel(mnthsLd)
                               aryTime((kk-1)*nYrsLd+1:kk*nYrsLd,:) = [(yrsLd(1):yrsLd(end))',mnthsLd(kk)*ones(nYrsLd,1)];
                            end
                        end

                        [dateBoth,~,~] = intersect(dateGCMVec(:,1:2),aryTime(:,1:2),'rows');
                        nMissing = numel(unique(aryTime(:,1:2),'rows')) - numel(dateBoth(:,1:2));
                        if nMissing ~= 0
                           warning('read_geodata:NcSpliceTsX', ...
                               [num2str(nMissing) ' requested NetCDF '...
                               'time-series elements were not loaded from ' ...
                               char(39) pathData char(39) '.' char(10) ...
                               'This may cause an error later in the script.']) 
                        end 
                    end
                end
            end
        else
            [sData] = NC_load_v2(pathData,sMeta,outTyp);
        end
    elseif isnan(yrsGcm(1)) || yrsGcm(1) == -9999 %Flag that data are not CMIP5
        [sData] = NC_load_v2(pathData,sMeta,outTyp);
    else
        error('read_geodata:yrsGcmUnknown',['The years present in the '...
            'current GCM file are unknown.']);
    end
    
    %None of files
    if ((min2d(yrsLd) < yrsGcm(1) && max2d(yrsLd) < yrsGcm(2)) || (min2d(yrsLd) > yrsGcm(1) && max2d(yrsLd) > yrsGcm(2))) && isempty(gcmUse)
        return
    end
    

    %Change fill values and missing values to NaN:
    indMissVal = strcmpi(sData.attData,'missing_value');
    if sum(sum(indMissVal)) > 0
        [misValRow, misValCol] = find(indMissVal == 1);
        sData.data(sData.data == sData.attData{misValRow, misValCol+1}) = NaN;
    end
    indMissVal = strcmpi(sData.attData,'_FillValue');
    if sum(sum(indMissVal)) > 0
        [misValRow, misValCol] = find(indMissVal == 1);
        sData.data(sData.data == sData.attData{misValRow, misValCol+1}) = NaN;
    end
    
    if ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'c'))
        [dateRef, tUnits] = NC_time_units(sData.attTime);

        if ~isempty(regexpi(tUnits,'days since'))
            cal =  NC_cal(sData.attTime);

            dateStrt = days_2_date(sData.time(1), dateRef, cal);
            dateEnd = days_2_date(sData.time(end), dateRef, cal);
        elseif ~isempty(regexpi(tUnits,'ka BP'))
            cal =  NC_cal(sData.attTime);
            
            dateStrt = time_kaBP_2_standard(sData.time(1), dateRef,cal);
            dateEnd = time_kaBP_2_standard(sData.time(end), dateRef,cal);
        else
            disp(['The GCM' char(39) 's time units are' tUnits '.']);
            error('NC_time:refUnknown',['A case has not been '...
                'written to deal with the current time '...
                'reference frame.']);
        end
        climTime = [dateStrt(1), dateEnd(1)];
        
        sData.attTime = [sData.attTime; {'type', 'climatology '}];
        for ii = 1 : length(mnthsLd)
            sData.attTime{end, 2}  = [sData.attTime{end, 2}, ...
               num2str(climTime(1)) '-' num2str(mnthsLd(ii)) 'thru' ...
            num2str(climTime(2)) '-' num2str(mnthsLd(ii))];
            if ii ~= length(mnthsLd)
                sData.attTime{end, 2}  = [sData.attTime{end, 2} '; '];
            end
        end
        
    else
        sData.attTime = [sData.attTime; {'type', 'time-series'}];
    end
    
    
    [gcmRef, gcmRefUnit] = NC_time_units(sData.attTime);
    if regexpbl(gcmRefUnit,{'days','since'},'and')
        aryTime = days_2_date(sData.time, gcmRef, NC_cal(sData.attTime));
    elseif regexpbl(gcmRefUnit,{'hours','since'},'and')
        aryTime = days_2_date(sData.time/24, gcmRef, NC_cal(sData.attTime));
%         %Convert from '0' hour to '24':
%         aryTime(aryTime(:,4) == 0,4) = 24;
    elseif regexpbl(gcmRefUnit,'unknown')
        aryTime = nan(3,1);
    else
        error('read_geodata:populateAryTime',['GCM units of ' gcmRefUnits ' are not recognized.']);
    end
    
    %Surmize timestep in data:
%     [~, gcmTimeUnits] = NC_time_units(sData.attTime);
%         timeTest = ncread(path,'time');
        dTime = find_att(sData.attTime,'time_step');
    
%     %Populate time array:
%     if numel(yrsLd) == 2
%         nYrsLd  = yrsLd(2) - yrsLd(1) + 1;
%     else
%         nYrsLd = numel(unique(yrsLd));
%     end
    
%     mnthsUniq = unique(mnthsLd);
%     nMnths = numel(mnthsUniq);
%     aryTime = nan(nYrsLd*nMnths,2);
%         
%     for ii = 1 : nMnths
%         aryTime((ii-1)*nYrsLd+ii : ii*nYrsLd,1) = (yrsLd(1):yrsLd(end))';
%         aryTime((ii-1)*nYrsLd+ii : ii*nYrsLd,2) = mnthsUniq(ii);
%     end

%Load WC binary files:
elseif ~isempty(regexpi(ext,'bil')) || ~isempty(regexpi(ext,'hdr')) %WorldClim
    if ~isempty(regexpi(sMeta.currVar,'tmp'))   %Necessary because WorldClim uses 'tmean' instead of 'tmn'.
        wcVar = 'tmean';
    elseif ~isempty(regexpi(sMeta.currVar,'tmn')) 
        wcVar = 'tmin';
    elseif ~isempty(regexpi(sMeta.currVar,'tmx')) 
        wcVar = 'tmax';
    elseif regexpi(sMeta.currVar,'pre')
        wcVar = 'prec';
    elseif regexpbl(sMeta.currVar,{'alt','elev','DEM'})
        wcVar = 'alt';
    else
        error('read_geodata:unknownWCvar',[sMeta.currVar ' is an '...
            'unknown variable for the WorldClim dataset. Program for this if correct.']);
    end
    
    %Check if WordClim path includes file or is general directory
    [root,file,extTemp] = fileparts(pathData);
    %If includes file, remove so that it searches directory 
    %(necessary because there are both hdr and data files):
    if ~isempty(extTemp)
        pathWC = root;
    else
        pathWC = fullfile(root,file);
    end
    
    %Check if this WorldClim dataset includes hyphen or not:
    if ~isempty(dir(fullfile(pathWC, strcat(wcVar, '_', '*.hdr' ) ) ))
        sepWC = '_';
    else
        sepWC = '';
    end
    
    %Find file for current month:
    if regexpbl(wcVar,'alt')
        fileWcHdr = dir(fullfile(pathWC, strcat(wcVar, sepWC, '.hdr' ) ) );
        fileWcData = dir(fullfile(pathWC, strcat(wcVar, sepWC, '.bil' ) ) );
        
        dTime = 'NA';
    else
        fileWcHdr = dir(fullfile(pathWC, strcat(wcVar, sepWC, num2str(mnthsLd), '.hdr' ) ) );
        fileWcData = dir(fullfile(pathWC, strcat(wcVar, sepWC, num2str(mnthsLd), '.bil' ) ) );
        
        dTime = 'climatology';
    end
    if ~isempty(fileWcHdr) && ~isempty(fileWcData)
        [hdrESRI, metaESRI] = read_wc_bin_hdr(fullfile(pathWC,fileWcHdr(1).name));
        if isfield(sMeta,'crd')
            [sData.data, hdrESRI, ~, ~] = read_wc_bin_data_v2(fullfile(pathWC,fileWcData(1).name), hdrESRI, sMeta.crd, sMeta.currVar, 'in');
        else
            [sData.data, ~, ~, ~] = read_wc_bin_data_v2(fullfile(pathWC,fileWcData(1).name), hdrESRI, [], sMeta.currVar, []);
        end
    else
        error('read_geodata:wcFileMissing',['Either the WordClim .hdr or .bil file associated with ' ...
            char(39) fileWcHdr char(39) ' was not found.']);
    end
    aryTime = [NaN, mnthsLd, NaN];
    climTime = [1950,2000];
    sData.info = {'source','WorldClim'};


elseif ~isempty(regexpi(ext,'asc'))
    %Read data file:
    [~,~,ext]= fileparts(pathData);
    if ~isempty(ext) %Path is file; load only this one
        [sData.data, hdrESRI, metaESRI] = read_ESRI(pathData);
        
        %Check if data loaded has time embedded in filename
        %(*_year_month.asc)
        indNum = regexpi(pathData,'\d');
        if ~isempty(indNum)
            indTime = indNum(([0,diff(indNum)] ~= 1));
            indUnd = regexpi(pathData,'[_.]');
            if numel(indTime) > 1 && numel(indUnd) > 1
                %ASSUMES MONTHLY DATA BECAUSE ONLY LOOKS FOR TWO NUMBERS:
                aryTime = [str2double(pathData(indTime(end-1):indUnd(end - 1)-1)), str2double(pathData(indTime(end):indUnd(end)-1))];
            end
        end
    else %Find all times in directory provided
        %Find files:
        if isnumeric(outTyp)
            [filesUse, aryTime, hdrESRI, metaESRI] = ASC_file_find(pathData, outTyp);
        else
            [filesUse, aryTime, hdrESRI, metaESRI] = ASC_file_find(pathData, yrsLd, mnthsLd);
        end

       %Load data:
       if ~isempty(filesUse)
           for ii = 1 : length(filesUse)
                [dataCurr, ~, ~] = read_ESRI(filesUse{ii});

                if ii == 1 %Initialize data array
                    sData.data = NaN(length(filesUse(:)),length(dataCurr(:,1)),length(dataCurr(1,:)));  
                end
                sData.data(ii,:,:) = dataCurr;
           end
       else
           sData.data = single.empty([0,0,0]);
           sData.lat = single.empty(0);
           sData.lon = single.empty(0);
           sData.time = single.empty(0);
           return
       end
    end

    if regexpbl(pathData,'cru')
        sData.info = {'source','CRU'};
        warning('read_geodata:CRUScale',['The current data are in a folder with CRU.  ' ...
            char(10) 'CRU Ascii files are scaled by a factor of ten.  The current data are NOT being scaled.'])
%         sData.data = sData.data/10; %All CRU data scaled by factor of ten
    elseif regexpbl(pathData,'worldclim')
        sData.info = {'source','WorldClim'};
        if regexpbl(sMeta.currVar,{'tmp','tas','tmin','tmax','tmn','tmx'})
            sData.data = sData.data/10; %WorldClim temperature data scaled by factor of ten
        end
    elseif regexpbl(pathData,'delta')
        sData.info = {'source','Downscaled data'};
    else
        sData.info = {'source','Unknown'};
    end
    
    dTime = 'month';
elseif isfield(sMeta,'hisTS') && ~isempty(regexpi(sMeta.hisTS,'will'))
    if ~isdir(pathData)
        warning('read_geodata:WMfiles','The path for W&M data is a file, which is being converted to a directory.')
        pathData = fileparts(pathData);
    end
    
    [sData.data, aryTime, hdrESRI, metaESRI] = read_Wil_Mat_data(pathData, yrsLd, mnthsLd);
    sData.info = {'source', 'Willmott and Matsuura'};

    dTime = 'month';
elseif isfield(sMeta,'hisTS') && ~isempty(regexpi(sMeta.hisTS,'gpcc'))
    if ~isdir(pathData)
        warning('read_geodata:GPCCfiles','The path for the GPCC data is a file, which is being converted to a directory.')
        pathData = fileparts(pathData);
    end
    
    [filesUse, aryTime] = GPCC_file_find(pathData, yrsLd, mnthsLd);
    
    for ii = 1 : length(filesUse(:))
        fileCurr = fullfile(pathData,filesUse{ii});
        
        [dataCurr, hdrESRI, metaESRI, ~, ~] = read_GPCC_ASCII(fileCurr);

        if ii == 1
            sData.data = NaN(length(filesUse(:)),length(dataCurr(:,1)),length(dataCurr(1,:)),'single'); 
        end
        sData.data(ii,:,:) = dataCurr;

    end
    sData.info = {'source','GPCC'};
    
    dTime = 'month';
elseif ~isempty(regexpi(ext,'tif')) %Load geotif data
    %Often GeoTIFFs are provided as regional tiles. Therefore, search for
    %other tiles in input:
    if isdir(pathData)
        foldData = pathData;
    else
        [foldData, ~, ~] = fileparts(pathData);
    end
    
    
    filesTemp = dir([foldData, filesep, '*', ext]);
    filesTemp = extractfield(filesTemp,'name');
    
    
    I = cell(numel(filesTemp),1);
    tempLon = [];
    tempLat = [];
    for ii = 1 : numel(filesTemp)
        I{ii} = GEOTIFF_READ(fullfile(foldData,filesTemp{ii}));
            I{ii}.y = I{ii}.y';
        tempLon = [tempLon, I{ii}.x];
        tempLat = [tempLat; I{ii}.y];
    end
    
    sData.lon = unique(tempLon);
    sData.lat = flipud(unique(tempLat));
    
    if isfield(sMeta,'crd') 
        sData.lon(sData.lon < sMeta.crd(1)) = [];
        sData.lon(sData.lon > sMeta.crd(2)) = [];
        sData.lat(sData.lat < sMeta.crd(3)) = [];
        sData.lat(sData.lat > sMeta.crd(4)) = [];
    end
    
    valNoData = -9999;
    sData.data = valNoData*ones([numel(sData.lat), numel(sData.lon)], 'like', I{1}.z);
    
    for ii = 1 : numel(I)
        %Check the ordering of the loaded geoTif:
        blLon = 0; %0 = ordered, 1 = flip, 2 = other
        blLat = 0; %0 = ordered, 1 = flip, 2 = other
        if ~issorted(I{ii}.x)
            if issorted(fliplr(I{ii}.x))
                blLon = 1;
            else
                blLon = 2;
            end
        end
        if ~issorted(flipud(I{ii}.y))
            if issorted(I{ii}.y)
                blLat = 1;
            else
                blLat = 2;
            end
        end
        
        
        %Find indices:
        indLonCurr = find(ismember(sData.lon, I{ii}.x));
        indLatCurr = find(ismember(sData.lat, I{ii}.y));
        
        if numel(indLonCurr) == numel(I{ii}.x) && numel(indLatCurr) == numel(I{ii}.y)
            if blLat == 0 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = I{ii}.z;
            elseif blLat == 1 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = flipud(I{ii}.z);
            elseif blLat == 0 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = fliplr(I{ii}.z);
            elseif blLat == 1 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = rot90(I{ii}.z,2);
            elseif blLon == 2 || blLat == 2
                error('read_geodata:geoTifOrder','This geoTif indice ordering has not been programmed for.');
            end
        elseif ~isempty(indLonCurr) && ~isempty(indLatCurr) 
            [indLonCurr, indLonDataCurr] = ismember(sData.lon, I{ii}.x);
            [indLatCurr, indLatDataCurr] = ismember(sData.lat, I{ii}.y);
            
            indLonCurr = find(indLonCurr == 1);
            indLatCurr = find(indLatCurr == 1);
            indLonDataCurr = find(indLonDataCurr == 1);
            indLatDataCurr = find(indLatDataCurr == 1);
            
            if blLat == 0 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = I{ii}.z(indLatDataCurr, indLonDataCurr);
            elseif blLat == 1 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = flipud(I{ii}.z(indLatDataCurr, indLonDataCurr));
            elseif blLat == 0 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = fliplr(I{ii}.z(indLatDataCurr, indLonDataCurr));
            elseif blLat == 1 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = rot90(I{ii}.z(indLatDataCurr, indLonDataCurr),2);
            elseif blLon == 2 || blLat == 2
                error('read_geodata:geoTifOrder','This geoTif indice ordering has not been programmed for.');
            end
        end
    end
    
    sData.data(sData.data == valNoData) = nan;

    outTyp = 'none';
    
elseif regexpbl(ext,{'img'}) %Load USGS IMG files 
    if ~regexpbl(filesTemp{1},'USGS')
       error('read_geodata:unknownIMG',['Currently, the only IMG ' ...
           'files that can be read are those from the USGS. ' char(10) ...
           'The reason is that it is known that the USGS file names ' ...
           'provide the position and resolution information.']);
    end
    
    %Find all files in folder:
    if isdir(pathData)
        foldData = pathData;
    else
        [foldData, ~, ~] = fileparts(pathData);
    end
    
    
    filesTemp = dir([foldData, filesep, '*', ext]);
    filesTemp = extractfield(filesTemp,'name');
    
    %LOAD ALL FILES IN FOLDER
    I = cell(numel(filesTemp),1);
    iLon = cell(numel(filesTemp),1);
    iLat = cell(numel(filesTemp),1);
    
    tempLon = [];
    tempLat = [];
    for ii = 1 : numel(filesTemp)
        indUnd = regexpi(filesTemp{ii}, '_');
        strCrd = filesTemp{ii}(indUnd(end-1)+1:indUnd(end)-1);
        
        if numel(strCrd) < 6
            strCrd = filesTemp{ii}(indUnd(end)+1:end-4);
        end
        
        res = filesTemp{ii}(indUnd(2)+1:indUnd(3)-1);
        if strcmpi(res, '2') %Image resolution is 2 arc seconds
            delta = 2/3600;
        else
            error('read_geodata:unknownImgRes',['The resolution ' res ...
                ' has not been programmed for.']);
        end
        
        indQuad = regexp(strCrd,'\D*');
        indOff = regexp(strCrd,'\d*');
        
        %Define latitude grid
        if strcmpi(strCrd(indQuad(1)), 'n')
            edgN = str2double(strCrd(indOff(1):indQuad(2)-1));
            edgS = edgN - 1; %1 degree tiles;
            iLat{ii} = (edgN - 0.5*delta : -delta : edgS + 0.5*delta)';
        elseif strcmpi(strCrd(indQuad(1)), 's')
            edgN = -str2double(strCrd(indOff(1):indQuad(2)-1));
            edgS = edgN - 1; %1 degree tiles;
            iLat{ii} = (edgN - 0.5*delta : -delta : edgS + 0.5*delta)';
        else
            error('read_geodata:unknownQuad',['The quadrant nomenclature ' ...
                strCrd(indQuad(1)) ' has not been programmed for.']);
        end
        
        %Define longitude grid:
        if strcmpi(strCrd(indQuad(2)), 'w')
            edgW = -str2double(strCrd(indOff(2):end));
            edgE = edgW + 1; %1 degree tiles;
            iLon{ii} = (edgW + 0.5*delta : delta : edgE - 0.5*delta);
        elseif strcmpi(strCrd(indQuad(2)), 'e')
            edgW = str2double(strCrd(indOff(2):end));
            edgE = edgW + 1; %1 degree tiles;
            iLon{ii} = (edgW + 0.5*delta : delta : edgE - 0.5*delta);
        else
            error('read_geodata:unknownQuad',['The quadrant nomenclature ' ...
                strCrd(indQuad(2)) ' has not been programmed for.']);
        end
        
        fid = fopen(fullfile(foldData, filesTemp{ii})); 
        I{ii} = fread(fid, [numel(iLat{ii}), numel(iLon{ii})]); 
        fclose(fid); 
        
        tempLon = [tempLon, iLon{ii}];
        tempLat = [tempLat; iLat{ii}];
    end
%     keyboard
    
    %Find all unique lat and lon points
    sData.lon = unique(tempLon);
    sData.lat = flipud(unique(tempLat));
    
    %Crop lat and lon pts to 
    if isfield(sMeta,'crd') 
        sData.lon(sData.lon < sMeta.crd(1)) = [];
        sData.lon(sData.lon > sMeta.crd(2)) = [];
        sData.lat(sData.lat < sMeta.crd(3)) = [];
        sData.lat(sData.lat > sMeta.crd(4)) = [];
    end
    
    valNoData = -9999;
    sData.data = valNoData*ones([numel(sData.lat), numel(sData.lon)], 'like', I{1});
        
    for ii = 1 : numel(I)
        %Check the ordering of the loaded geoTif:
        blLon = 0; %0 = ordered, 1 = flip, 2 = other
        blLat = 0; %0 = ordered, 1 = flip, 2 = other
        if ~issorted(iLon{ii})
            if issorted(fliplr(iLon{ii}))
                blLon = 1;
            else
                blLon = 2;
            end
        end
        if ~issorted(flipud(iLat{ii}))
            if issorted(iLat{ii})
                blLat = 1;
            else
                blLat = 2;
            end
        end
        
        
        %Find indices:
        indLonCurr = find(ismember(sData.lon, iLon{ii}));
        indLatCurr = find(ismember(sData.lat, iLat{ii}));
        
        if numel(indLonCurr) == numel(iLon{ii}) && numel(indLatCurr) == numel(iLat{ii})
            if blLat == 0 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = I{ii};
            elseif blLat == 1 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = flipud(I{ii});
            elseif blLat == 0 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = fliplr(I{ii});
            elseif blLat == 1 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = rot90(I{ii},2);
            elseif blLon == 2 || blLat == 2
                error('read_geodata:geoTifOrder','This geoTif indice ordering has not been programmed for.');
            end
        elseif ~isempty(indLonCurr) && ~isempty(indLatCurr) 
            [indLonCurr, indLonDataCurr] = ismember(sData.lon, iLon{ii});
            [indLatCurr, indLatDataCurr] = ismember(sData.lat, iLat{ii});
            
            if all(indLonCurr == 1 | indLonCurr == 0)
                indLonCurr = find(indLonCurr == 1);
                indLatCurr = find(indLatCurr == 1);
            end
            if all(indLonDataCurr == 1 | indLonDataCurr == 0)
                indLonDataCurr = find(indLonDataCurr == 1);
                indLatDataCurr = find(indLatDataCurr == 1);
            end
            
            if blLat == 0 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = I{ii}(indLatDataCurr, indLonDataCurr);
            elseif blLat == 1 && blLon == 0
                sData.data(indLatCurr, indLonCurr) = flipud(I{ii}(indLatDataCurr, indLonDataCurr));
            elseif blLat == 0 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = fliplr(I{ii}(indLatDataCurr, indLonDataCurr));
            elseif blLat == 1 && blLon == 1
                sData.data(indLatCurr, indLonCurr) = rot90(I{ii}(indLatDataCurr, indLonDataCurr),2);
            elseif blLon == 2 || blLat == 2
                error('read_geodata:geoTifOrder','This geoTif indice ordering has not been programmed for.');
            end
        end
    end
    
    sData.data(sData.data == valNoData) = nan;
    
    outTyp = 'none';
    
else
	error('read_data:noType',['The data in ' pathData ' is of an unknown type.']); 
end

%If aryTime is only month and year, add day (mid-month)
if ~regexpbl(outTyp, 'none') %Only do if time data being processed
    if numel(aryTime(1,:)) == 2 && sum2d(isnan(aryTime)) == 0
        aryTime = [aryTime(:,1:2), round(eomday(aryTime(:,1),aryTime(:,2))/2)];
%     elseif numel(aryTime(1,:)) ~= 3 && sum2d(isnan(aryTime)) == 0
%         error('read_geodata:nAryTime',['The array ' char(39) 'aryTime' char(39) ...
%             ' has ' num2str(numel(aryTime(1,:))) ' entries, which has not been coded for.' ]);
    end
end

if isempty(regexpi(ext,'nc')) && isempty(regexpi(ext,'grb'))
    if isfield(sMeta,'crd')
        %Clip Data:
        if isempty(regexpi(ext,'bil')) && isempty(regexpi(ext,'hdr')) && isempty(regexpi(ext,'tif'))  && isempty(regexpi(ext,'img')) %WorldClim data already clipped
            [~, bndsInd] = adj_bounds(hdrESRI, sMeta.crd, 'out');
            hdrESRI = adj_grid_hdr(hdrESRI, bndsInd);
            if length(size(sData.data)) == 3
                sData.data = sData.data( : ,bndsInd(4) : bndsInd(3) , ...
                    bndsInd(1) : bndsInd(2));
            else
                sData.data = sData.data(bndsInd(4) : bndsInd(3) , ...
                    bndsInd(1) : bndsInd(2));
            end
        end
    end

    %Write sData.lat:
    if (~isfield(sData,'lon') || (isfield(sData,'lon') && all(isnan(sData.lon)))) || (~isfield(sData,'lat') || (isfield(sData,'lat') && all(isnan(sData.lat))))
        [sData.lat, sData.lon] = ESRI_hdr2geo(hdrESRI,metaESRI);
    end
    sData.attLat = {'long_name', 'latitude'; 'units', 'degrees_north'};
    sData.attLon = {'long_name', 'longitude'; 'units', 'degrees_east'};

    %Write sData.time:
    if isnumeric(outTyp) || (~isnumeric(outTyp) && isempty(regexpi(outTyp,'none')))
        strCal = 'gregorian';
        daysSince = [1900,1,1];
        sData.attTime = {'long_name', 'time'; 'units', ...
            ['days since ' num2str(daysSince(1)) '-' num2str(daysSince(2)) '-' num2str(daysSince(3))]; ...
            'calendar', strCal; 'type','time-series'; 'time_step', dTime};
        if ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'c')) %Demarcates climatology
            if isempty(climTime)
               climTime = [aryTime(1,1), aryTime(end,1)];
            end

            sData.time = nan(1,1);
            sData.attTime{end,2} = 'climatology ';
            for ii = 1 : length(mnthsLd)
                sData.attTime{end, 2}  = [sData.attTime{end, 2}, ...
                   num2str(climTime(1)) '-' num2str(mnthsLd(ii)) 'thru' ...
                num2str(climTime(2)) '-' num2str(mnthsLd(ii))];
                if ii ~= length(mnthsLd)
                    sData.attTime{end, 2}  = [sData.attTime{end, 2} '; '];
                end
            end
        else %Demarcates time-series
            if ~issorted(aryTime,'rows')
                [aryTime,indTime] = sortrows(aryTime);
                sData.data = sData.data(indTime,:,:);
            end
            if sum2d(~isnan(aryTime)) == 0
                sData.time = nan(numel(aryTime(:,1)),1);
            else
                sData.time = days_since(daysSince, aryTime, strCal);
            end
        end
    end

    %Write Attributes:
    if isfield(sMeta,'currVar')
        if regexpbl(sMeta.currVar,'pre')
            units = 'mm';
        elseif regexpbl(sMeta.currVar,{'tmp','tmn','tmx'})
            units = 'Celsius';
        elseif regexpbl(sMeta.currVar,{'alt','DEM','elev'})
            units = 'm';
        else
            units = 'Unknown';
        end
        sData.var = sMeta.currVar;
        sData.attData = {'long_name', sMeta.currVar; 'units', units; ...
            '_FillValue', -9999; 'missing_value', -9999};  
    elseif regexpbl(pathData,{'pr','pre','tmp','tmean','tas','tmn','tmin','tasmin','tmx','tmax','tasmax'})
        if regexpbl(pathData,{'pr','pre'})
            units = 'mm';
            sData.var = 'pre';
        elseif regexpbl(pathData, {'tmp','tmean','tas','tmn','tmin','tasmin','tmx','tmax','tasmax'})
            units = 'Celsius';
            if regexpbl(pathData, {'tmp','tmean','tas'})
                sData.var = 'tmp';
            elseif regexpbl(pathData, {'tmn','tmin','tasmin'})
                sData.var = 'tmn';
            elseif regexpbl(pathData, {'tmx','tmax','tasmax'})
                sData.var = 'tmx';
            end
        end
        sData.attData = {'long_name', sData.var; 'units', units; ...
            '_FillValue', -9999; 'missing_value', -9999}; 
        
    else
        sData.attData = {'long_name', 'unknown'; 'units', 'unknown'; ...
            '_FillValue', -9999; 'missing_value', -9999}; 
    end
end

sizeData = size(sData.data);
if length(sizeData) == 2
    sizeData = [1, sizeData];
end

%Try to check that all output is as expected:
if ~isnumeric(outTyp) && isempty(regexpi(outTyp,'none')) && ~all(isnan(yrsLd)) && ~all(isnan(mnthsLd))
    if regexpbl(dTime,'month')
        nTsLd = (yrsLd(2) -  yrsLd(1) + 1)*length(mnthsLd);

        %If Climatology requested, check size and avg if necessary
        if ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'c'))
            %If data loaded is not already a climatology ensure correct number of elements loaded and average
            if sizeData(1) > 1 
                days_w_dates_cmpr(sData,yrsLd,mnthsLd); 
                sData.data = squeeze(mean(sData.data,1));
                sData.time = nan(1,1);
            end
        elseif ~isnumeric(outTyp) && isempty(regexpi(outTyp,'c'))  %Check if correct number of time-series elements loaded:
            sData.data = squeeze(sData.data);
            if nTsLd ~= 1 && sizeData(1) == 1 && ~isnan(nTsLd)
                error('read_geodata:wrongDim',['The data matrix loaded has '...
                    num2str(length(size(sData.data))) ' dimensions but 3 '...
                    'are expected.']);
            end
            days_w_dates_cmpr(sData,yrsLd,mnthsLd); 
        end 
    else
        dateFilled = date_vec_fill([yrsLd(1), mnthsLd(1), 1], [yrsLd(end), mnthsLd(end), eomday(yrsLd(end), mnthsLd(end))], NC_cal(sData.attTime));
        if numel(dateFilled(:,1)) > numel(sData.time)
            warning('read_geodata:daysMissing',['The number of dates '...
                'requested appears to total ' ...
                num2str(numel(dateFilled(:,1))) ' but the number '...
                'loaded equals ' num2str(numel(sData.time))]);
        end
    end
end

if regexpbl(class(sData.data), 'double')
    sData.data = single(sData.data);
end


%Check .lon column vector and .lat row vector:
if length(sData.lon(:,1)) > 1 && length(sData.lon(1,:)) == 1 
    sData.lon = sData.lon';
end
if length(sData.lat(1,:)) > 1 && length(sData.lat(:,1)) == 1 
    sData.lat = sData.lat';
end

%%Check data oriented as if looking at map (North on top, East on right)
if ~issorted(sData.time)
    [sData.time, indSortTime] = sort(sData.time);
    sData.data = sData.data(indSortTime,:,:);
end
if ~issorted(flipud(sData.lat))
    [sData.lat, indSortLat] = sort(sData.lat,'descend');
    sData.data = sData.data(:,indSortLat,:);
end
if ~issorted(sData.lon)
    [sData.lon, indSortLon] = sort(sData.lon);
    sData.data = sData.data(:,:,indSortLon);
end


%%Check dates loaded with dates requested:
if ~isempty(sData.time) && ~isnan(sData.time(1))
    %Find dates loaded:
    [dateRefChck, tUnits] = NC_time_units(sData.attTime);

    if regexpbl(tUnits,'day')
        dateLoad = days_2_date(sData.time, dateRefChck, NC_cal(sData.attTime));
    elseif regexpbl(tUnits,'hour')
        dateLoad = days_2_date(sData.time/24, dateRefChck, NC_cal(sData.attTime));
    else
        error('read_geodata:unknownTUnit', [char(39) tUnits char(39) ' are unknown time units. Code for this.']);
    end
    
    if regexpbl(dTime,{'month'})
        datePresent = intersect(aryTime(:,1:2),dateLoad(:,1:2),'rows');
    elseif regexpbl(dTime,{'day'})
        datePresent = intersect(aryTime(:,1:3),dateLoad(:,1:3),'rows');
    elseif regexpbl(dTime,'hour')
        datePresent = intersect(aryTime(:,1:4),dateLoad(:,1:4),'rows');
    else
        error('read_geodata:unknownTStep', [char(39) dTime char(39) ' is an unknown time step. Code for this.']);
    end
    
    if numel(datePresent(:,1)) ~= numel(dateLoad(:,1)) || numel(datePresent(:,1)) ~= numel(aryTime(:,1)) 
        error('read_geodata:wrngTime',['The data loaded have ' num2str(numel(datePresent(:,1))) ...
            ' entries, which do not exactly match the ' ...
            num2str(numel(datePresent(:,1))) ' entries requested.' char(10) ...
            'Check that the data present correspond with the downscaling '...
            'dates requested.']);
    end
end


if ~regexpbl(varargin,'no_disp')
    %Display success message:
    if isfield(sMeta,'currTime')
        disp([mnth_str(sMeta.currTime(2)) ' data have finished being loaded from ' char(39) ...
            pathData char(39) '.' char(10)]);
    else
        disp(['Data have finished being loaded from ' char(39) ...
            pathData char(39) '.' char(10)]);
    end
end

