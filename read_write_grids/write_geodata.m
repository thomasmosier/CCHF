function write_geodata(pathData,sData,sMeta,prec,wrtTyp,varargin)
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


%Input: %Designed to work with function 'read_geodata.m'
%'path' = full path for output file
%'sData' isstructure array with format:
    %sData.var = Name of the climate variable in the structure
    %sData.info = Cell array of general dataset information 
    %sData.lat = vector of latitudes
    %sData.attLat = cell array with attribute list for latitude coordinate
    %sData.lon = vector of longitudes
    %sData.attLon = cell array with attribute list for longitude coordinate
    %sData.time = vector of times corresponding to each data grid loaded
    %sData.attTime = cell array with attribute list for time coordinate
    %sData.data = 2 or 3 dimensional data array with climate data
    %sData.attData = cell array with attribute list for data
%'sMeta' is a structure array describing data properties:
    %sMeta.currVar = variable to load
    %sMeta.hisTS = 'GPCC', 'CRU', or 'Willmott'
    %sMeta.yrsOut = [first year of output data, last year of output data]
    %sMeta.yrsClim = [first year of climatology, last year of climatology]
    %sMeta.mnths = list of months requested
    %sMeta.crd = [lonW, lonE, latS, latN]
    %sMeta.currTime = [year, month]
%'prec' = number of decimals to use in writing data
%'wrtTyp' = ASCII or netCDF

if ~regexpbl(varargin,'no_disp')
    disp([mnth_str(sMeta.currTime(2)) ' data are now being written to ' ...
        char(39) pathData char(39) '.']);
end

%Create output directory if it doesn't already exist:
[dirData, file, ext] = fileparts(pathData);
if ~exist(dirData,'dir')
    mkdir(dirData);
end

if ~isempty(regexpi(wrtTyp,'asc'))
    [root,file,~] = fileparts(pathData);
    pathData = fullfile(root,[file,'.asc']);
    
    if isempty(prec)
        if ~isempty(regexpi(char(sData.var),'pre'))
            prec = 0;
        else
            prec = 1;
        end
    end
    
    if length(size(sData.data)) > 2
        nIter = length(sData.data(:,1,1));
    else
        nIter = 1;
    end
    
    for ii = 1 : nIter
        if nIter > 1
            %Find current time:
            if sum(size(sData.time(ii,:)) == 2) == 1 %time data likely stored as [year, month]
                vecTCurr = sData.time(ii,:);
                if vecTCurr(1) < vecTCurr(2)
                    vecTCurr(1:2) = fliplr(vecTCurr(1:2));
                end
            elseif sum(size(sData.time(ii,:)) == 1) == 2 %Likely the time data is recorded as 'days since ...'
                [dateRef, tUnits] = NC_time_units(sData.attTime);
                
                if ~isempty(regexpi(tUnits,'days since'))
                    cal =  NC_cal(sData.attTime);
                    
                    vecTCurr = days_2_date(sData.time(ii), dateRef, cal);
                else
                    disp(['The GCM' char(39) 's time units are' tUnits '.']);
                    error('NC_time:refUnknown',['A case has not been '...
                        'written to deal with the current time '...
                        'reference frame.']);
                end
            else
                error('write_geodata:unknownTime',...
                    ['The dimensions of the time vector are ' ...
                    num2str(length(sData.time(:,1))) ' rows by ' ...
                    num2str(length(sData.time(1,:))) ' columns.'  char(10) ...
                    'A case to interpret this has not been written.'])
            end
            
            indDate = regexpi(file,'_');
            if ~isempty(indDate)
                file = [file(1:indDate(end-1)), num2str(vecTCurr(1)), ...
                    '_' num2str(vecTCurr(2)), ext];
            end
            pathCurr = fullfile(dirData, file);
        else
            pathCurr = pathData;
        end
        
        %Check that step sizes are all equal.  
        step = nanmean(abs(horzcat(reshape(diff(sData.lon),1,[]), reshape(diff(sData.lat),1,[]))));
        %If stepsizes not equal, resample data:

        if round(100*(step - abs(diff(sData.lat(1:2))))/step) ~= 0 || round(100*(step - abs(diff(sData.lon(1:2))))/step) ~= 0
            %Display warning alerting about resampling 
            if ii == 1
                warning('write_geodata:avgStep',['The data currently being '...
                    'written have a non-uniform step size and will '...
                    'therefore be resampled to a uniform grid.  This will cause inaccurate '...
                    'representation.  To write data using non-uniform '...
                    'step-size, write in NetCDF format.']);
                
                %ADJUST LONGITUDE
                %Create longitude vector and parameters for:
                nLonPts = numel(sData.lon);
                stepAlign = 0.001*step;
                lonEdgO = box_edg(sData.lon);

                nLonAlignPts = round(2*step/stepAlign);
                lonLap = nan(nLonAlignPts,1);
                if sData.lon(1) < sData.lon(end)
                    for jj = 1 : nLonAlignPts
                        lonR1 = sData.lon(1) - step + jj*stepAlign;
                        lonLap(jj) = sum(abs(lonEdgO - box_edg(lonR1 : step : lonR1 + (nLonPts-1)*step)));  
                    end
                elseif sData.lon(1) > sData.lon(end)
                    for jj = 1 : nLonAlignPts
                        lonR1 = sData.lon(1) + step - jj*stepAlign;
                        lonLap(jj) = sum(abs(lonEdgO - box_edg(lonR1 : - step : lonR1 - (nLonPts-1)*step)));  
                    end
                else
                    error('write_geodata:lonOne','Case with unequal step size and only one indice not written for.')
                end

                %Find minimum displacement between original and new longitude
                %grids:
                [valLon, indLonLap] = min(lonLap);
                if valLon > nLonPts*0.5*step
                    warning('write_geodata:lonResample',['The total difference between longitude centroids values is ' num2str(valLon) ', which is larger than expected.']);
                end

                %Define new longitude grid:
                if sData.lon(1) < sData.lon(end)
                    lonR1 = sData.lon(1) + step - indLonLap*stepAlign;
                    lonR = lonR1 : step : lonR1 + (nLonPts-1)*step;
                elseif sData.lon(1) > sData.lon(end)
                    lonR1 = sData.lon(1) - step + indLonLap*stepAlign;
                    lonR = lonR1 : -step : lonR1 - (nLonPts-1)*step;
                else
                    error('write_geodata:lonOne','Case with unequal step size and only one indice not written for.');
                end

                %ADJUST LATITUDE
                %Create latitude vector:
                nLatPts = numel(sData.lat);
                latEdgO = box_edg(sData.lat);
                nLatAlignPts = round(2*step/stepAlign);
                latLap = nan(nLatAlignPts,1);

                if sData.lat(1) < sData.lat(end)
                    for jj = 1 : nLatAlignPts
                        latR1 = sData.lat(1) - step + jj*stepAlign;
                        latLap(jj) = sum(abs(latEdgO - box_edg(latR1 : step : latR1 + (nLatPts-1)*step)'));  
                    end
                elseif sData.lat(1) > sData.lat(end)
                    for jj = 1 : nLatAlignPts
                        latR1 = sData.lat(1) + step - jj*stepAlign;
                        latLap(jj) = sum(abs(latEdgO - box_edg(latR1 : -step : latR1 - (nLatPts-1)*step)'));  
                    end
                else
                    error('write_geodata:lonOne','Case with unequal step size and only one indice not written for.')
                end

                %Find minimum displacement between original and new longitude
                %grids:
                [valLat, indLatLap] = min(latLap);
                if valLat > nLatPts*0.5*step
                    warning('write_geodata:lonResample',['The total difference between longitude centroids values is ' num2str(valLat) ', which is larger than expected.']);
                end

                if sData.lat(1) < sData.lat(end)
                    latR1 = sData.lat(1) - step + indLatLap*stepAlign;
                    latR = (latR1 : step : latR1 + (nLatPts-1)*step)';
                elseif sData.lat(1) > sData.lat(end)
                    latR1 = sData.lat(1) + step - indLatLap*stepAlign;
                    latR = (latR1 : -step : latR1 - (nLatPts-1)*step)';

                else
                    error('write_geodata:latOne','Case with unequal step size and only one indice not written for.')
                end
            end
            
            
            
            %RESAMPLE DATA TO NEW GRID:
            if length(size(sData.data)) == 3 && numel(latR) > 1 && numel(lonR) > 1
                dataCurr = area_int_2D_v2(sData.lon,sData.lat,sData.data(ii,:,:),lonR,latR);
            elseif length(size(sData.data)) == 2 && numel(latR) > 1 && numel(lonR) > 1
                dataCurr = area_int_2D_v2(sData.lon,sData.lat,sData.data(:,:),lonR,latR);
            else
                dataCurr = NaN;
            end
            
            hdr = ESRI_hdr(lonR, latR, 'corner');
        else %Step sizes equal
            if length(size(sData.data)) == 3
                dataCurr = squeeze(sData.data(ii,:,:));
            elseif length(size(sData.data)) == 2
                dataCurr = sData.data;
            else
                error('write_geodata:noData',['The number of grids found ' ...
                    'to write is ' num2str(nIter) '.  A case has not been '...
                    'written for this.']);
            end
            
            hdr = ESRI_hdr(sData.lon, sData.lat, 'corner');
        end
        
        write_ESRI_v4(squeeze(dataCurr), hdr, pathCurr, prec);
    end
elseif regexpbl(wrtTyp,{'netCDF','nc'})
    if isempty(ext) %If extension empty, create filename
        %Remove existing numbers from data path, if any, and replace with dates 
        %corresponding to all output:

        extNC = '.nc';
        indDate = regexpi(file,'_');
        if ~isempty(indDate)
            file = file(1:indDate(end-1)-1);
        end

        %Determine if data are climatology:
        [blClim, yrsClim,~] = NC_isclim(sData);
        if blClim == 1
            yrsWrt = yrsClim;
            file = [file, '_clim'];
        else
            yrsWrt = sMeta.yrsOut;
        end

        %Ensure month string is two digits:
        if sMeta.mnths(1) < 10
            mnthStart = ['0' num2str(sMeta.mnths(1))];
        else
            mnthStart = num2str(sMeta.mnths(1));
        end
        if sMeta.mnths(end) < 10
            mnthEnd = ['0' num2str(sMeta.mnths(end))];
        else
            mnthEnd = num2str(sMeta.mnths(end));
        end

        %Create full path
        file = [file, '_' num2str(yrsWrt(1)), mnthStart '-', ...
            num2str(yrsWrt(2)), mnthEnd, extNC];
        pathCurr = fullfile(dirData, file);
    else
        pathCurr = pathData;
    end
    
    %Write data to NetCDF file:
    write_NC(pathCurr, sData, prec);
    
    %%THIS SECTION CASUES MATLAB TO CRASH EVERY TIME
%     %If current time includes ending time-series element, ensure NetCDF
%     %data are in chronological order:
%     if sMeta.currTime(2) == sMeta.mnths(end)
%         %Read data from NetCDf file:
%         vecTime = ncread(pathCurr, 'time');
%         if ~issorted(vecTime)
%             keyboard
%             %Open NetCDF file and retrieve variable attributes:
%             fidNC = netcdf.open(pathCurr,'WRITE');
%             vidData = netcdf.inqVarID(fidNC,char(sData.var));
%             vifTime = netcdf.inqVarID(fidNC,'time');
%             aryData = netcdf.getVar(fidNC, vidData);
% 
%             %Sort time and data:
%             [vecTime, ordSrt] = sort(vecTime);
%             aryData = aryData(ordSrt,:,:);
% 
%             %Write re-ordered time and data vectors:
%             netcdf.putVar(fidNC,vifTime, 0, numel(vecTime), vecTime);
%             netcdf.putVar(fidNC,vidData, [0, 0, 0], size(aryData), aryData);
%             netcdf.close(fidNC);
%         end
%     end
else
   error('write_data:noType',['The data output format chosen is ' ...
       wrtTyp ', which is an unkown type.']); 
end


if ~regexpbl(varargin,'no_disp')
    disp([mnth_str(sMeta.currTime(2)) ' data have finished being written '...
    'to ' char(39) pathData char(39) '.' char(10)]);
end


end