function write_NC(pathData, sData, varargin)
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


%Creates NetCDF geo-referenced climate data file using data and metadata
%provided by inputs.

%INPUTS:
%path = full filepath (including extension) for desired output .nc file
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

    
%DATA ACTUALLY NEED TO HAVE THREE DIMENIONS
% %Ensure data does not have extra dimensions:
% sData.data = squeeze(sData.data);
    
%Select data type to use (default is single precision):
if iscell(varargin) && ~isempty(varargin)
    prec = varargin{1};
    if prec == 0
        dataTyp = 'int16';
        filVal = intmax(dataTyp); %Missing data value (all NaN values go to this):
    elseif prec > 0 && prec < 5
        dataTyp = 'single';
        filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
    else
        dataTyp = 'double';
        filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
    end
elseif ~isempty(regexpi(char(sData.var),{'pre','pr'}))
    dataTyp = 'uint16';
    filVal = intmax(dataTyp); %Missing data value (all NaN values go to this):
else
    dataTyp = 'single';
    filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
end    



% %Initialize:
% ovrWrt = '';

% %If file exists and dimensions are different, ask if it should be
% %overwritten:
% if exist(pathData,'file')
%     
%     ds = ncdataset(pathData);
%     [gcmLat, ~] = NC_lat(ds);
%     [gcmLon, ~] = NC_lon(ds);
%     if ~isequal(sData.lon, gcmLon) || ~isequal(sData.lat,gcmLat)
%         pathDisp = strrep(pathData, filesep, [filesep filesep]);
%     	ovrWrt = input(['The file ' pathDisp ' exists with different ' ...
%             'dimensions than the data being written to it.' char(10) ...
%             'Overwrite this file (Y or N)'],'s'); %#ok<NASGU>
%         if isempty(regexpi(ovrWrt,'y'))
%             disp(['The file ' pathData ' is not being overwritten and the data ' ...
%                 'is therefore not being stored to a file.']);
%             return
%         elseif ~isempty(regexpi(ovrWrt,'y'))
%             clear ds
%             delete(pathData)
%         end
%     end
% end
    



%Define NetCDF File Variables or indices for writing extra data:
if ~exist(pathData,'file')
    strtTime = 1;
    strtData = [1,1,1];
%     strtData = ones(numel(size(sData.data)),1);
    
    strFormat = 'netcdf4';     
    %Not sure why to use 'format' = 'classic' vs. 'format' = 'netcdf4' 
    %(if using the latter, 'ChunkSize' = length of each dimension)
    
    %Define attributes:
    nccreate(pathData, char(sData.var), 'Dimensions', ...
        {'time' Inf 'lat' length(sData.lat(:)) 'lon' length(sData.lon(:))}, ...
        'Format',strFormat, 'Datatype', dataTyp, 'FillValue', filVal); 
    
%     %Redefine chunksize:
%     ncid = netcdf.open(pathData,'NC_WRITE');
%     varid = netcdf.inqVarID(ncid,char(sData.var));
%     netcdf.defVarChunking(ncid,varid,'CONTIGUOUS'); 
%     netcdf.close(ncid);
    
    nccreate(pathData,'time','Dimensions',{'time' Inf}, ...
        'Format',strFormat, 'Datatype', 'single');
    
    if isfield(sData,'time_bnds')
        nccreate(pathData,'time_bnds','Dimensions',{'ntb' 2 'time' Inf}, ...
            'Format',strFormat, 'Datatype', 'single');
    end
    
    nccreate(pathData,'lat', 'Dimensions',{'lat' length(sData.lat(:))}, ...
        'Format',strFormat, 'Datatype', 'double');
    nccreate(pathData,'lon','Dimensions',{'lon' length(sData.lon(:))}, ...
        'Format',strFormat, 'Datatype', 'double'); 
else
    %Find current number of datapoints in NetCDF file:
    infoData = ncinfo(pathData, char(sData.var));
    strtData = infoData.Size;
    
    infoTime = ncinfo(pathData, 'time');
    strtTime = infoTime.Size;
    
    strtTime = strtTime + 1;
    strtData = [strtData(1) + 1, 1, 1];
%     strtData = [strtData(1) + ones(numel(size(sData.data)),1)];
end

%Determine if data are climatology
[blClim, yrsClim, mnthClim] = NC_isclim(sData);
%If climatology, sData.time = NaN

if blClim == 1
    aryTime = NaN(1,3);
        aryTime(1,1) = int8(mean(yrsClim));
        aryTime(1,2) = int8(mnthClim);
        aryTime(1,3) = eomday(aryTime(1), aryTime(2));
        vecTime = days_since([1900,1,1], [aryTime(:,1), aryTime(:,2), round(aryTime(:,3)/2)], 'gregorian');
else
    %Create time vector and define time reference (if unit conversion occurs)
    if sum(size(sData.time) == 1) == 0
        if size(sData.time(1,:)) == 2
            aryTime = NaN(length(sData.time(:,1)),3);
            aryTime(:,1:2) = int8(sData.time);
            aryTime(:,3) = eomday(aryTime(:,1),aryTime(:,2));
            vecTime = days_between([1900,1,1], [aryTime(:,1), aryTime(:,2), round(aryTime(:,3)/2)], 'gregorian');
            %Write reference date in time attributes:
            timeUntTemp = strcmpi(sData.attTime,'units');
            if sum(timeUntTemp) == 0
                sData.attTime{end+1,1} = {'units', 'days since 1900-01-01'};
            else
                [untRow, untCol] = find(timeUntTemp == 1);
                sData.attTime{untRow, untCol+1} = 'days since 1900-01-01';
            end
            %Write calendar used to time attributes:
            timeCalTemp = strcmpi(sData.attTime,'calendar');
            if sum(timeCalTemp) == 0
                sData.attTime{end+1,1} = {'calendar', 'gregorian'};
            else
                [calRow, calCol] = find(timeCalTemp == 1);
                sData.attTime{calRow, calCol+1} = 'gregorian';
            end

        else
            error('write_NC:unknownTime',['The dimensions of the time vector are ' ...
                num2str(length(sData.time(:,1))) ' rows by ' ...
                num2str(length(sData.time(1,:))) ' columns.'  char(10) ...
                'A case to interpret this has not been written.'])
        end
    elseif sum(size(sData.time) == 1) == 1 || sum(size(sData.time) == 1) == 2
        vecTime = sData.time;
    else
        error('write_NC:unknownTime',['The dimensions of the time vector are ' ...
                num2str(length(sData.time(:,1))) ' rows by ' ...
                num2str(length(sData.time(1,:))) ' columns.'  char(10) ...
                'A case to interpret this has not been written.'])
    end 
end

%Check if time elements are nan:
if sum(isnan(vecTime)) > 0
   warning('write_NC:nanTime',[num2str(sum(isnan(vecTime))) ' NaN time '...
       'values have been detected.  This may cause the program to crash.']); 
end

%Only write attributes if file is new:
if strtTime == 1
    %Write global attributes:
    if isfield(sData,'info')
        for ii = 1 : length(sData.info(:,1))
            ncwriteatt(pathData,'/',sData.info{ii,1},sData.info{ii,2}) %'/' indicates global attribute
        end
    end

    %Write time variable attributes:
    for ii = 1 : length(sData.attTime(:,1))
        ncwriteatt(pathData,'time',sData.attTime{ii,1},sData.attTime{ii,2})
    end

    %Write latitude variable attributes:
    for ii = 1 : length(sData.attLat(:,1))
        ncwriteatt(pathData,'lat',sData.attLat{ii,1},sData.attLat{ii,2})
    end

    
    %Write longitude variable attributes:
    for ii = 1 : length(sData.attLon(:,1))
        ncwriteatt(pathData,'lon',sData.attLon{ii,1},sData.attLon{ii,2})
    end

    %Write data variable attributes:
    %Remove '_FillValue' (it's written above) from data attributes
    indFil = nan;
    for ii = 1 : length(sData.attData)
        if ~isempty(regexpi(sData.attData{ii,1},'fill'))
           indFil = ii; 
        end
    end

    if ~isnan(indFil)
        sData.attData(indFil,:) = [];
    end

    %change 'missing_value' in data attributes to be the same as the fill value
    indMis = nan;
    for ii = 1 : length(sData.attData)
        if ~isempty(regexpi(sData.attData{ii,1},'missing'))
           indMis = ii; 
        end
    end

    if ~isnan(indMis)
        sData.attData{indMis,2} = filVal;
    end

    %Write data attributes:
    for ii = 1 : length(sData.attData(:,1))
        ncwriteatt(pathData,char(sData.var),sData.attData{ii,1},sData.attData{ii,2}) %'/' indicates global attribute
    end
end

%Change NaN data to fill value:
sData.data( isnan(sData.data) ) = filVal;

%%Create and write variables to NetCDF file:
%Write time vector:
ncwrite(pathData,'time', vecTime, strtTime);

if isfield(sData,'time_bnds')
    ind2 = find(size(sData.time_bnds) == 2);
    if numel(ind2) > 0
        if numel(ind2) == 1 && ind2 == 2
            sData.time_bnds = sData.time_bnds';
        end
        ncwrite(pathData,'time_bnds', sData.time_bnds, [1, strtTime]); 
    else
       warning('write_NC:dimT_bnds',['Variable ' char(39) 'time_bnds' ...
           char(39) ' is not being written because it does not have proper dimensions.']); 
    end

end
    
%Only write attributes if file is new:
if strtTime == 1
    %Latitude:
    ncwrite(pathData,'lat',sData.lat);
    
    %Longitude:
    ncwrite(pathData,'lon',sData.lon);
end

%Data needs to have 3 indices (if has 2, add 1)
if length(size(sData.data)) == 2
    dataTemp = nan([1,size(sData.data)]);
    dataTemp(1,:,:) = sData.data;
    sData.data = dataTemp;
end

%Write Data:
if regexpbl(dataTyp,'int16')
    sData.data = int16(sData.data);
elseif regexpbl(dataTyp,'single')
    sData.data = single(sData.data);
end

ncwrite(pathData, char(sData.var), sData.data, strtData, [1,1,1]);