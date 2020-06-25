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

function load_ts_equal(sPath, sHydro, dateCurr, sMeta)
%loads any month of data whose filename embedded dates matches
%tsDate and clips the grid according to the input hdr data.

global sAtm

%INPUTS:
%dataPath = path of data to load
%tsDate = vector of date [year, month, day]
%sHydro = structure of with entries
    %sHydro.(varLat) = latitudes to load
    %sHydro.(varLon) = longitudes to load

%OUTPUTS:
%sDataO = geodata structure of form:
    %sDataO.info = Cell array of general dataset information 
    %sDataO.(varLat) = vector of latitudes
    %sDataO.attLat = cell array with attribute list for latitude coordinate
    %sDataO.(varLon) = vector of longitudes
    %sDataO.attLon = cell array with attribute list for longitude coordinate
    %sDataO.time = vector of times corresponding to each data grid loaded
    %sDataO.attTime = cell array with attribute list for time coordinate
    %sDataO.data = 2 or 3 dimensional data array with climate data
    %sDataO.attData = cell array with attribute list for data
    %sDataO.var = Name of the climate variable in the structure
   
varLon = 'longitude';
varLat = 'latitude';

datesUse = sMeta.dateRun;
if iscell(datesUse)
    if isfield(sMeta, 'siteCurr')
        datesUse = datesUse{sMeta.siteCurr};
    else
        error('loadTsEqual:noSiteCurr', ['The date field is a cellarray. '...
            'This requires the presence of a siteCurr field to determine which index to use.']);
    end
end


%Exit function if all requested data already loaded:
if ~isempty(dateCurr) && ~all(isnan(dateCurr))
    for ll = 1 : numel(sMeta.varLd) %Loop over all variables to load
        
        varCurr = sMeta.varLd{ll};
        
        %Don't read data if no path to current variable 
        if ~isfield(sPath,varCurr) || isempty(sPath.(varCurr))
            continue
        end  
        
        %Check if sAtm already contains needed data
        %If so, skip rest of this function (do not re-load data)
        if isfield(sAtm, varCurr) && isfield(sAtm, ['date' varCurr])
            indDateLded = ismember(sAtm.(['date' varCurr]), datesUse(sMeta.indCurr,:),'rows');
            if any(indDateLded == 1) %Expected time-series for full duration of model run
                %Update time indice for current variable and skip to next variable
                sAtm.(['ind' varCurr]) = find(indDateLded == 1, 1, 'first');
                continue
            elseif sAtm.(['date' varCurr])(sAtm.(['ind' varCurr]), 1) ~= datesUse(sMeta.indCurr, 1) %Check abnormal case (such as one year of representative time-series that are repeated)
                [~, fileCurr, ~] = fileparts(sPath.([varCurr 'File']){sMeta.indCurr});
                timeCheck = CMIP5_time(fileCurr);
                
                %If the file that is supposed to be loaded is for a
                %different year, but has same other time elements as
                %current model run time
                if ~any(timeCheck(:,1) == datesUse(sMeta.indCurr, 1)) && any(timeCheck(:,1) == sAtm.(['date' varCurr])(1, 1))
                    indDateLded = ismember(sAtm.(['date' varCurr])(:,2:end), datesUse(sMeta.indCurr,2:end),'rows');
                    if any(indDateLded == 1) 
                        %Update time indice for current variable and skip to next variable
                        sAtm.(['ind' varCurr]) = find(indDateLded == 1, 1, 'first');
                        continue
                    end
                end
            end
        end
        
          
        %READ DATA:
        %Use date indexing if defined previously (faster than finding 
        %file within function by providing function the date).
        if isfield(sPath, [varCurr 'File']) %Field written with function 'path_data_ind'
            pathCurr = sPath.([varCurr 'File']){sMeta.indCurr};
            %If climate data aligns with grid, special field will have been written for sPath (inside of "engine_v*") 
            if regexpbl(pathCurr, '.nc')
                nmsTemp = extract_field(ncinfo(pathCurr), 'Variables');
                varInFile = extract_field(nmsTemp{1}, 'Name');
                
                varRead = varCurr;
                if ~any(strcmpi(varInFile, varCurr))
                    if regexpbl(varCurr, 'pre')
                        if any(strcmpi(varInFile, 'pr'))
                            varRead = 'pr';
                        else
                            error('load_ts_equal:noPrePresent',['Precip ' ...
                                varCurr ' alternative names have not been programmed for. Add in statement above.']);
                        end
                    elseif strcmpi(varCurr, 'pr')
                        if any(strcmpi(varInFile, 'pre'))
                            varRead = 'pre';
                        else
                            error('load_ts_equal:noPrePresent',['Precip ' ...
                                varCurr ' alternative names have not been programmed for. Add in statement above.']);
                        end
                    elseif regexpbl(varCurr, 'tmean')
                        if any(strcmpi(varInFile, 'tas'))
                            varRead = 'tas';
                        else
                            error('load_ts_equal:noPrePresent',['Precip ' ...
                                varCurr ' alternative names have not been programmed for. Add in statement above.']);
                        end
                    elseif strcmpi(varCurr, 'tas')
                        if any(strcmpi(varInFile, 'tmean'))
                            varRead = 'tmean';
                        elseif any(strcmpi(varInFile, 'tmp'))
                            varRead = 'tmp';
                        else
                            error('load_ts_equal:noPrePresent',['Temperature ' ...
                                varCurr ' alternative names have not been programmed for. Add in statement above.']);
                        end
                    elseif strcmpi(varCurr, 'tasmax')
                        if any(strcmpi(varInFile, 'tasmax'))
                            varRead = 'tasmax';
                        elseif any(strcmpi(varInFile, 'tmax'))
                            varRead = 'tmax';
                        elseif any(strcmpi(varInFile, 'tmx'))
                            varRead = 'tmx';
                        else
                            error('load_ts_equal:noPrePresent',['Temperature ' ...
                                varCurr ' alternative names have not been programmed for. Add in statement above.']);
                        end
                    elseif strcmpi(varCurr, 'tasmin')
                        if any(strcmpi(varInFile, 'tasmin'))
                            varRead = 'tasmin';
                        elseif any(strcmpi(varInFile, 'tmin'))
                            varRead = 'tmin';
                        elseif any(strcmpi(varInFile, 'tmn'))
                            varRead = 'tmn';
                        else
                            error('load_ts_equal:noPrePresent',['Temperature ' ...
                                varCurr ' alternative names have not been programmed for. Add in statement above.']);
                        end
                    elseif strcmpi(varCurr, 'bcdep')
                        if any(strcmpi(varInFile, 'bcdep'))
                            varRead = 'bcdep';
                        elseif any(strcmpi(varInFile, 'totaldeposition'))
                            varRead = 'totaldeposition';
                        elseif any(strcmpi(varInFile, 'bc'))
                            varRead = 'bc';
                        elseif any(strcmpi(varInFile, 'oc'))
                            varRead = 'oc';
                        elseif any(strcmpi(varInFile, 'dust'))
                            varRead = 'dust';
                        else
                            error('load_ts_equal:noPrePresent',['Temperature ' ...
                                varCurr ' alternative names have not been programmed for. Add in statement above.']);
                        end
                    else
                        error('load_ts_equal:noVarPresent',['Variable ' ...
                            varCurr ' alternative names have not been programmed for. Add in statement above.']);
                    end
                end
                
                
                %Assign supporting fields from hydrologic grid
                sDataCurr.(varLon) = sHydro.(varLon);
                sDataCurr.(varLat) = sHydro.(varLat);
                %Assign supporting fields from NetCDf file
                sDataCurr.time    = ncread(pathCurr,'time');
                sDataCurr.attTime = ncinfo(pathCurr,'time');
                    sDataCurr.attTime = squeeze(struct2cell(sDataCurr.attTime.Attributes))';
                sDataCurr.attData = ncinfo(pathCurr, varRead);
                    sDataCurr.attData = squeeze(struct2cell(sDataCurr.attData.Attributes))';
                %Assign data field from NetyCDF file
                sDataCurr.data = single(ncread(pathCurr,varRead));
                
                %Bug fix for specific case where leap day written
                %incorrectly:
                if numel(sDataCurr.time) == 29 && all(diff(sDataCurr.time(1:28)) == 1) && diff(sDataCurr.time(28:29)) > 10000
                    %nYrs = nYrs (same as last value)
                    sDataCurr.time(29) = sDataCurr.time(28)+1;
                end
            elseif regexpbl(pathCurr, {'.txt','.asc'})
                sDataCurr.(varLon) = sHydro.(varLon);
                sDataCurr.(varLat) = sHydro.(varLat);

                %Look for marker of time-res
                if isfield(sPath, [varCurr 'Res'])
                    if regexpbl(sPath.([varCurr 'Res']), {'day', 'daily'})
                        testDigit = 3;
                    elseif regexpbl(sPath.([varCurr 'Res']), {'month'})
                        testDigit = 2;
                    end
                else
                    %Figure out time resolution
                    indUnd = regexpi(sPath.([varCurr 'File']){sMeta.indCurr}, '_');
                    testDigit = zeros(numel(indUnd) + 1, 1);
                    for ii = 1 : numel(testDigit)
                        if ii == 1 
                            if ~isnan(str2double(sPath.([varCurr 'File']){sMeta.indCurr}(1:indUnd(1)-1)))
                                testDigit(ii) = 1; 
                            end
                        elseif ii == numel(testDigit) 
                            indExt = regexpi(sPath.([varCurr 'File']){sMeta.indCurr}, '\.');
                            if ~isnan(str2double(sPath.([varCurr 'File']){sMeta.indCurr}(indUnd(end)+1:indExt(end)-1)))
                                testDigit(ii) = 1; 
                            end
                        else
                            if ~isnan(str2double(sPath.([varCurr 'File']){sMeta.indCurr}(indUnd(ii-1)+1:indUnd(ii)-1)))
                                testDigit(ii) = 1;
                            end
                        end
                    end
                end


                indUnd = regexpi(sPath.([varCurr 'File']){sMeta.indCurr}, '_');
                indExt = regexpi(sPath.([varCurr 'File']){sMeta.indCurr}, '\.');
                if sum(testDigit) == 2
                    dateCurr = ...
                        [str2double(sPath.([varCurr 'File']){sMeta.indCurr}(indUnd(end-1)+1:indUnd(end)-1)), ...
                        str2double(sPath.([varCurr 'File']){sMeta.indCurr}(indUnd(end)+1:indExt(end)-1))];
                elseif sum(testDigit) == 3
                    dateCurr = ...
                        [str2double(sPath.([varCurr 'File']){sMeta.indCurr}(indUnd(end-2)+1:indUnd(end-1)-1)), ...
                        str2double(sPath.([varCurr 'File']){sMeta.indCurr}(indUnd(end-1)+1:indUnd(end)-1)), ...
                        str2double(sPath.([varCurr 'File']){sMeta.indCurr}(indUnd(end)+1:indExt(end)-1))];
                end

                %Create NC style time att (necessary for calendar, etc.)
                sDataCurr.attTime = {'calendar', sMeta.cal; 'units', 'days_since 1000-01-01'};
                sDataCurr.time = days_since([1000,1,1], dateCurr, sMeta.cal);

                sDataCurr.data = nan([1, numel(sHydro.(varLat)), numel(sHydro.(varLon))],'single');
                [sDataCurr.data(1,:,:), ~, ~] = read_ESRI(sPath.([varCurr 'File']){sMeta.indCurr});

                %If precipitation and units unknown, assume they are mm but print warning:
                if regexpbl(varCurr, 'pr')
                    if isequal(dateCurr, datesUse(1,:))
                        warning('load_ts:preConvert',['Units of the '...
                            'precipitation input are assumed to be '...
                            'mm, which are being converted to m.']);
                    end
                    sDataCurr.data = sDataCurr.data / 1000;
                end

            else
                sMeta.dateCurr = datesUse(sMeta.indCurr,:);
                sDataCurr = read_geodata(sPath.([varCurr 'File']){sMeta.indCurr}, sMeta, 'no_disp');
                sDataCurr.data = single(sDataCurr.data);
            end
        else
            sDataCurr = read_geodata(sPath.(varCurr), sMeta, dateCurr, 'no_disp');
            sDataCurr.data = single(sDataCurr.data);
        end
          
        %Set variable field of sDataCurr struct:
        sDataCurr.var = varCurr;
            
        %Check to see that data loaded and not nan
        if isempty(sDataCurr.data) || (all(all2d(isnan(sDataCurr.data))) && all(isnan(sDataCurr.time) == 0))
            if all2d(isnan(sDataCurr.data))
               warning('load_ts:allNaN','Data were loaded, but they are all NaN');
            end

%             %If previous data provided as input, assign it to output:
%             if ~isempty(varargin) && ~isempty(varargin{1})
%                 if ll == 1 %Create array without data, then re-populate it
%                     sInput = rmfield_x(sLoaded,{'data','pr','pre','tmp','tmn','tmx','tas','tmin','tmax'});
%                 end
%                 sInput.(varCurr) = sLoaded.(varCurr);
%             else
%                 sInput.(varCurr) = nan;
%             end
            
            continue
        end

        %Transform data to lat-lon grid used with other datasets:
        if nanmean(abs(sDataCurr.(varLat)-sHydro.(varLat))) > sMeta.spEps || nanmean(abs(sDataCurr.(varLon)-sHydro.(varLon))) > sMeta.spEps
        %if ~isequal(round2(sDataCurr.(varLat),4),round2(sHydro.(varLat),4)) || ~isequal(round2(sDataCurr.(varLon),4),round2(sHydro.(varLon),4 ))
            warning('one_month_data_load:regrid', ...
                ['Data loaded from ' char(39) ...
                strrep(sPath.(varCurr),'\',':') char(39) ...
                ' are being regridded.  This adds processing time.']);
 
            %Initialize new grid:
            sDataCurr.dataTemp = nan(squeeze(numel(sDataCurr.data(:,1,1))), numel(sHydro.(varLat)), numel(sHydro.(varLon)));

            %Regrid:
            for ii = 1 : numel(sDataCurr.data(:,1,1))
                sDataCurr.dataTemp(ii,:,:) = area_int_2D_v2(sDataCurr.(varLon), sDataCurr.(varLat), squeeze(sDataCurr.data(ii,:,:)), sHydro.(varLon), sHydro.(varLat),'nansum');
            end

            %Transfer regular data field:
            sDataCurr.data = sDataCurr.dataTemp;
            sDataCurr = rmfield(sDataCurr,'dataTemp');

            %Replace lat and lon:
            sDataCurr.(varLon) = sHydro.(varLon);
            sDataCurr.(varLat) = sHydro.(varLat);
        end

        
        %Change units to metric:
        if isfield(sDataCurr,'attData')
            [units, rowUnit] = NC_units(sDataCurr.attData);
            
            if isfield(sDataCurr,'var') && regexpbl(sDataCurr.var,{'pr'}) %Convert pre units to meters
                if strcmpi(units, 'mm') 
                    sDataCurr.attData{rowUnit,2} = 'm';
                    sDataCurr.data = sDataCurr.data / 1000;
                elseif strcmpi(units, 'mm/month') 
                    if regexpbl(sMeta.dt, 'month')
                        sDataCurr.attData{rowUnit,2} = 'm';
                        sDataCurr.data = sDataCurr.data / 1000;
                    else
                        error('month_load:preUnitsMnth',['The precipitation units are ' units ' and the model time step is ' sMeta.dt '. This time mismatch has not been programmed for.']) 
                    end
                elseif strcmpi(units, 'mm/day')
                    if regexpbl(sMeta.dt, {'day','daily'})
                        sDataCurr.attData{rowUnit,2} = 'm';
                        sDataCurr.data = sDataCurr.data / 1000;
                    else
                        error('month_load:preUnitsDay',['The precipitation units are ' units ' and the model time step is ' sMeta.dt '. This time mismatch has not been programmed for.']) 
                    end
                elseif ~strcmpi(units,'m')
                   error('month_load:preUnits',['The precipitation units are ' units ', which have not been coded for.']) 
                end    
            elseif isfield(sDataCurr,'var') && regexpbl(sDataCurr.var,{'tmp','temp','tmean','tas','tmn','tmin','tasmin','tmx','tmax','tasmax'}) 
                if regexpbl(units,'kelvin') || strcmpi(units, 'K')
                    sDataCurr.data = sDataCurr.data - 273.15;
                    sDataCurr.attData{rowUnit, 2} = 'Celsius';
                elseif ~regexpbl(units,'celsius')
                    error('month_load:tempUnits',['Temperature units are ' units ' which are not recognized.']);
                end
                
            end
        end


        [refDateNew, ~] = NC_time_units(sDataCurr.attTime);
        calNew =  NC_cal(sDataCurr.attTime);
        
        
%         if isequal(dateCurr(2:3), [2,1]) || isequal(dateCurr(2:3), [3,1])
%             keyboard
%         end
        %%Splice with previous data:
        if ~isfield(sAtm,'attTime')
            sAtm.attTime = sDataCurr.attTime;
            sAtm.time = sDataCurr.time;
        end
        
        %Find reference times for both already loaded and presently
        %loaded datasets:
        [refDateOld, ~] = NC_time_units(sAtm.attTime);

        %Find time as date vectors for both already loaded and presently
        %loaded datasets
        datesNew = days_2_date_v2(sDataCurr.time, refDateNew, calNew);

        if regexpbl(sMeta.dt, {'day', 'daily'})
            datesNew = floor(datesNew(:,1:3));
        elseif regexpbl(sMeta.dt, 'month')
            datesNew = floor(datesNew(:,1:2));
        end
        
        %Convert dates of previously loaded time-series from original calendar to calendar of new data:
        sAtm.time = days_since(refDateOld, datesNew, calNew);

        sAtm.(['date' varCurr]) = datesNew;
        
        if ll == 1 && ~isfield(sAtm,'lat')
            sAtm.(varLat) = sDataCurr.(varLat);
            sAtm.(varLon) = sDataCurr.(varLon);
        end

        %Transfer data from current variable to combined data structure
        %array:
        sAtm.(varCurr) = sDataCurr.data;
        indDateLded = ismember(sAtm.(['date' varCurr]), datesUse(sMeta.indCurr,:), 'rows');
        if any(indDateLded)
            sAtm.(['ind' varCurr]) = find(indDateLded == 1, 1, 'first');
        else %This logic in the else statement is meant to catch time-averaged climate input variables (e.g. bcdeposition)
            filesUniqCurr = unique(sPath.([varCurr 'File']));
            filesUniqPr = unique(sPath.('prFile'));
            
            blTimeAvg = 0;
            if (numel(filesUniqCurr) < numel(filesUniqPr) && mode(sAtm.(['date' varCurr])(:,1)) ~= datesUse(sMeta.indCurr, 1))
                blTimeAvg = 1;
            elseif ~isequal(unique(sAtm.(['date' varCurr])(1,:)), unique(sAtm.('datepr')(1,:))) && strcmpi(varCurr, 'bcdep')
                blTimeAvg = 1;
            end
 
            if blTimeAvg == 1
                indDateLded = ismember(sAtm.(['date' varCurr])(:,2:end), datesUse(sMeta.indCurr,2:end), 'rows');
                if any(indDateLded)
                    sAtm.(['ind' varCurr]) = find(indDateLded == 1, 1, 'first');
                else
                    error('load_ts_equal:noIndAvg',['At time indice ' num2str(sMeta.indCurr) ...
                   ' no ' varCurr ' date was found. This is for the case of an average climate variable.']);
                end
            else 
                error('load_ts_equal:noIndTs',['At time indice ' num2str(sMeta.indCurr) ...
                   ' no ' varCurr ' date was found. This is for the case of a time-series climate variable.']);
            end 
        end
        
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
    
         