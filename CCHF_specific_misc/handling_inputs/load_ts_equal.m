function load_ts_equal(sPath, sHydro, dateCurr, sMeta)
%loads any month of data whose filename embedded dates matches
%tsDate and clips the grid according to the input hdr data.

global sAtm

%INPUTS:
%dataPath = path of data to load
%tsDate = vector of date [year, month, day]
%sHydro = structure of with entries
    %sHydro.lat = latitudes to load
    %sHydro.lon = longitudes to load

%OUTPUTS:
%sDataO = geodata structure of form:
    %sDataO.info = Cell array of general dataset information 
    %sDataO.lat = vector of latitudes
    %sDataO.attLat = cell array with attribute list for latitude coordinate
    %sDataO.lon = vector of longitudes
    %sDataO.attLon = cell array with attribute list for longitude coordinate
    %sDataO.time = vector of times corresponding to each data grid loaded
    %sDataO.attTime = cell array with attribute list for time coordinate
    %sDataO.data = 2 or 3 dimensional data array with climate data
    %sDataO.attData = cell array with attribute list for data
    %sDataO.var = Name of the climate variable in the structure
   


%Exit function if all requested data already loaded:
if ~isempty(dateCurr) && ~all(isnan(dateCurr))
    for ll = 1 : numel(sMeta.varLd) %Loop over all variables to load
        %Don't read data if no path to current variable 
        if ~isfield(sPath,sMeta.varLd{ll}) || isempty(sPath.(sMeta.varLd{ll}))
            continue
        end    
          
        %READ DATA:
        %Use date indexing if defined previously (faster than finding 
        %file within function by providing function the date).
        if isfield(sPath, [sMeta.varLd{ll} 'File']) %Field written with function 'path_data_ind'

            %If climate data aligns with grid, special field will have been written for sPath (inside of "backbone_v*") 
            if isfield(sPath,'varNm')
                if regexpbl(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr},'.nc')
                    sDataCurr.lon = sHydro.lon;
                    sDataCurr.lat = sHydro.lat;
                    sDataCurr.time = ncread(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr},'time');
                    sDataCurr.attTime = ncinfo(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr},'time');
                        sDataCurr.attTime = squeeze(struct2cell(sDataCurr.attTime.Attributes))';
                    sDataCurr.attData = ncinfo(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr},sPath.varNm{ll});
                        sDataCurr.attData = squeeze(struct2cell(sDataCurr.attData.Attributes))';
                    indUnit = strcmpi(sDataCurr.attData,'units');
                        [rUnit, cUnit] = find(indUnit == 1);
                        dataUnit = sDataCurr.attData{rUnit, cUnit+1};
                    sDataCurr.data = ncread(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr},sPath.varNm{ll});

                    if strcmpi(dataUnit,'mm')
                        sDataCurr.data = sDataCurr.data / 1000;
                    elseif strcmpi(dataUnit,'K') || strcmpi(dataUnit,'kelvin')
                        sDataCurr.data = sDataCurr.data - 273.15;
                    end
                elseif regexpbl(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr},{'.txt','.asc'})
                    sDataCurr.lon = sHydro.lon;
                    sDataCurr.lat = sHydro.lat;
                    
                    %Look for marker of time-res
                    if isfield(sPath, [sMeta.varLd{ll} 'Res'])
                        if regexpbl(sPath.([sMeta.varLd{ll} 'Res']), {'day', 'daily'})
                            testDigit = 3;
                        elseif regexpbl(sPath.([sMeta.varLd{ll} 'Res']), {'month'})
                            testDigit = 2;
                        end
                    else
                        %Figure out time resolution
                        indUnd = regexpi(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}, '_');
                        testDigit = zeros(numel(indUnd) + 1, 1);
                        for ii = 1 : numel(testDigit)
                            if ii == 1 
                                if ~isnan(str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(1:indUnd(1)-1)))
                                    testDigit(ii) = 1; 
                                end
                            elseif ii == numel(testDigit) 
                                indExt = regexpi(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}, '\.');
                                if ~isnan(str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(indUnd(end)+1:indExt(end)-1)))
                                    testDigit(ii) = 1; 
                                end
                            else
                                if ~isnan(str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(indUnd(ii-1)+1:indUnd(ii)-1)))
                                    testDigit(ii) = 1;
                                end
                            end
                        end
                    end
                    
                    
                    indUnd = regexpi(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}, '_');
                    indExt = regexpi(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}, '\.');
                    if sum(testDigit) == 2
                        dateCurr = ...
                            [str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(indUnd(end-1)+1:indUnd(end)-1)), ...
                            str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(indUnd(end)+1:indExt(end)-1))];
                    elseif sum(testDigit) == 3
                        dateCurr = ...
                            [str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(indUnd(end-2)+1:indUnd(end-1)-1)), ...
                            str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(indUnd(end-1)+1:indUnd(end)-1)), ...
                            str2double(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}(indUnd(end)+1:indExt(end)-1))];
                    end
                    
                    %Create NC style time att (necessary for calendar, etc.)
                    sDataCurr.attTime = {'calendar', 'Gregorian'; 'units', ['days_since 1000-01-01']};
                    sDataCurr.time = days_since([1000,1,1],dateCurr, 'Gregorian');
                    
                    sDataCurr.data = nan([1, numel(sHydro.lat), numel(sHydro.lon)],'single');
                    [sDataCurr.data(1,:,:), ~, ~] = read_ESRI(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr});
                    
                    %If precipitation, assumes units are mm:
                    if regexpbl(sMeta.varLd{ll}, 'pr')
                        if isequal(dateCurr, sMeta.dateRun(1,:))
                            warning('load_ts:preConvert',['Units of the '...
                                'precipitation input are assumed to be '...
                                'mm, which are being converted to m.']);
                        end
                        sDataCurr.data = sDataCurr.data / 1000;
                    end
                        
                else
                    error('load_ts:unknownFileType',['The current file type is not recognized. '...
                        'Best option is to not use ' char(39) ...
                        'path_find_files' char(39) ', which will initiate use of read_geodata']);
                end
            else
                sDataCurr = read_geodata(sPath.([sMeta.varLd{ll} 'File']){sMeta.indCurr}, sMeta, 'no_disp');
            end
        else
            sDataCurr = read_geodata(sPath.(sMeta.varLd{ll}), sMeta, dateCurr, 'no_disp');
        end
            
            
        
        %Check to see that data loaded and not nan
        if isempty(sDataCurr.data) || (all2d(isnan(sDataCurr.data)) && all(isnan(sDataCurr.time) == 0))
            if all2d(isnan(sDataCurr.data))
               warning('load_ts:allNaN','Data were loaded, but they are all NaN');
            end

%             %If previous data provided as input, assign it to output:
%             if ~isempty(varargin) && ~isempty(varargin{1})
%                 if ll == 1 %Create array without data, then re-populate it
%                     sInput = rmfield_x(sLoaded,{'data','pr','pre','tmp','tmn','tmx','tas','tmin','tmax'});
%                 end
%                 sInput.(sMeta.varLd{ll}) = sLoaded.(sMeta.varLd{ll});
%             else
%                 sInput.(sMeta.varLd{ll}) = nan;
%             end
            
            continue
        end

        %Transform data to lat-lon grid used with other datasets:
        if ~isequal(round2(sDataCurr.lat,4),round2(sHydro.lat,4)) || ~isequal(round2(sDataCurr.lon,4),round2(sHydro.lon,4 ))
            warning('one_month_data_load:regrid', ...
                ['Data loaded from ' char(39) ...
                strrep(sPath.(sMeta.varLd{ll}),'\',':') char(39) ...
                ' are being regridded.  This adds processing time.']);

            %Initialize new grid:
            sDataCurr.dataTemp = nan(squeeze(numel(sDataCurr.data(:,1,1))), numel(sHydro.lat), numel(sHydro.lon));

            %Regrid:
            for ii = 1 : numel(sDataCurr.data(:,1,1))
                sDataCurr.dataTemp(ii,:,:) = area_int_2D_v2(sDataCurr.lon, sDataCurr.lat, squeeze(sDataCurr.data(ii,:,:)), sHydro.lon, sHydro.lat,'nansum');
            end

            %Transfer regular data field:
            sDataCurr.data = sDataCurr.dataTemp;
            sDataCurr = rmfield(sDataCurr,'dataTemp');

            %Replace lat and lon:
            sDataCurr.lon = sHydro.lon;
            sDataCurr.lat = sHydro.lat;
        end

        %Check units:
        if isfield(sDataCurr,'attData')
            [units, rowUnit] = NC_units(sDataCurr.attData);
            %Change units to metric:
            if isfield(sDataCurr,'var') && regexpbl(sDataCurr.var,{'pr'}) %Convert pre units to meters 
                if strcmpi(units,'mm') 
                    sDataCurr.attData{rowUnit,2} = 'm';
                    sDataCurr.data = sDataCurr.data / 1000;
                elseif ~strcmpi(units,'m')
                   error('month_load:preUnits',['The precipitation units are ' units ', which have not been coded for.']) 
                end    
            elseif isfield(sDataCurr,'var') && regexpbl(sDataCurr.var,{'tmp','tmn','tmx'}) && ~regexpbl(units,{'celsius','kelvin'})
                warning('month_load:tempUnits',['Temperature units are ' units ' which are not recognized.']);
            end
        end


        [refDateNew, ~] = NC_time_units(sDataCurr.attTime);
        calNew =  NC_cal(sDataCurr.attTime);
        
%         if isfield()
%         dateOld = 
        
        %%Splice with previous data:
        if ~isfield(sAtm,'attTime')
            sAtm.attTime = sDataCurr.attTime;
            sAtm.time = sDataCurr.time;
        else
            %Find reference times for both already loaded and presently
            %loaded datasets:
            [refDateOld, ~] = NC_time_units(sAtm.attTime);

            %Find time as date vectors for both already loaded and presently
            %loaded datasets
            datesNew = days_2_date(sDataCurr.time, refDateNew, calNew);

            %Convert dates of previously loaded time-series from original calendar to calendar of new data:
            sAtm.time = days_since(refDateOld, datesNew, calNew);
            sAtm.datecurr = datesNew;
        end
        
        if ll == 1 && ~isfield(sAtm,'lat')
            sAtm.lat = sDataCurr.lat;
            sAtm.lon = sDataCurr.lon;
        end

        %Transfer data from current variable to combined data structure
        %array:
        sAtm.(sMeta.varLd{ll}) = sDataCurr.data;
%         if isequal(sAtm.time, sDataCurr.time)
%             
%         else
%             keyboard
%             error('load_ts:timeMismatch', ['The time entries for the '...
%                 'present loaded time-series to not match. Possibly code '...
%                 'case for this (i.e. use subset of time indices?']);
%         end
    end 
else
    return 
end
    
         