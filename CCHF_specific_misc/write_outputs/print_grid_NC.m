function print_grid_NC(pathData, data, var, lon, lat, timeVec, varargin)
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
    %info = Cell array of general dataset information 
    %lat = vector of latitudes
    %attLat = cell array with attribute list for latitude coordinate
    %lon = vector of longitudes
    %attLon = cell array with attribute list for longitude coordinate
    %timeVec = vector of times corresponding to each data grid loaded
    %attTime = cell array with attribute list for time coordinate
    %data = 2 or 3 dimensional data array with climate data
    %attData = cell array with attribute list for data
    %var = Name of the climate variable in the structure

    
%DATA ACTUALLY NEED TO HAVE THREE DIMENIONS
% %Ensure data does not have extra dimensions:
% data = squeeze(data);

%Select data type to use (default is single precision):
if iscell(varargin) && ~isempty(varargin)
    prec = varargin{1};
    if prec == -1
        dataTyp = 'uint16';
        filVal = intmax(dataTyp); %Missing data value (all NaN values go to this):
    elseif prec == 0
        dataTyp = 'int16';
        filVal = intmax(dataTyp); %Missing data value (all NaN values go to this):
    elseif prec > 0 && prec < 5
        dataTyp = 'single';
        filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
    else
        dataTyp = 'double';
        filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
    end
elseif ~isempty(regexpi(char(var),{'pre','pr'}))
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
%     if ~isequal(lon, gcmLon) || ~isequal(lat,gcmLat)
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
    


% keyboard

%Define NetCDF File Variables or indices for writing extra data:
if ~exist(pathData,'file')
    [pathCrt, ~, ~]= fileparts(pathData);
    if ~exist(pathCrt,'dir')
        mkdir(pathCrt);
    end
    
    strtTime = 1;
    strtData = [1,1,1];
%     strtData = ones(numel(size(data)),1);
    
    strFormat = 'netcdf4';     
    %Not sure why to use 'format' = 'classic' vs. 'format' = 'netcdf4' 
    %(if using the latter, 'ChunkSize' = length of each dimension)
    
    %Define attributes:
    nccreate(pathData, var, 'Dimensions', ...
        {'time' Inf 'lat' length(lat(:)) 'lon' length(lon(:))}, ...
        'Format',strFormat, 'Datatype', dataTyp, 'FillValue', filVal); 
    
%     %Redefine chunksize:
%     ncid = netcdf.open(pathData,'NC_WRITE');
%     varid = netcdf.inqVarID(ncid,char(var));
%     netcdf.defVarChunking(ncid,varid,'CONTIGUOUS'); 
%     netcdf.close(ncid);
    
    nccreate(pathData,'time','Dimensions',{'time' Inf}, ...
        'Format',strFormat, 'Datatype', 'single');
    
    nccreate(pathData,'time_bnds','Dimensions',{'ntb' 2 'time' Inf}, ...
        'Format',strFormat, 'Datatype', 'single');     
    
    nccreate(pathData,'lat', 'Dimensions',{'lat' length(lat(:))}, ...
        'Format',strFormat, 'Datatype', 'double');
    nccreate(pathData,'lon','Dimensions',{'lon' length(lon(:))}, ...
        'Format',strFormat, 'Datatype', 'double'); 
else
    %Find current number of datapoints in NetCDF file:
    infoData = ncinfo(pathData, char(var));
    strtData = infoData.Size;
    
    infoTime = ncinfo(pathData, 'time');
    strtTime = infoTime.Size;
    
    strtTime = strtTime + 1;
    strtData = [strtData(1) + 1, 1, 1];
%     strtData = [strtData(1) + ones(numel(size(data)),1)];
end


switch numel(timeVec(1,:))
    case 1
        refDate = [1900,1,1];
    case 2
        refDate = [1900,1,1];
    case 3
        refDate = [1900,1,1];
    case 4
        refDate = [1900,1,1,0];
    otherwise
        error('print_grid_NC:timeVec',[num2str(numel(timeVec(1,:))) ...
            ' is an unexpected number of elements in the time vector.']);
end
    
timeDays = days_since(refDate, timeVec, 'gregorian');


attData = {'missing_value', filVal ...
    };

attTime = {'units', 'days since 1900-01-01'; ...
    'calendar', 'gregorian'; ...
    'bounds', 'time_bnds'; ...
    'axis', 'T'; ...
    'long_name', 'time'; ...
    'standard_name', 'time' ...
    };

if numel(refDate) == 4
    attTime{1,2} = 'days since 1900-01-01-0';
end

attLon = {'units', 'degrees_east'; ...
    'axis', 'X'; ...
    'long_name', 'longitude'; ...
    'standard_name', 'longitude'};

attLat = {'units', 'degrees_north'; ...
    'axis', 'Y'; ...
    'long_name', 'latitude'; ...
    'standard_name', 'latitude' ...
    };


%Check if time elements are nan:
if sum(isnan(timeDays)) > 0
   warning('write_NC:nanTime',[' NaN time '...
       'values have been detected.  This may cause the program to crash.']); 
end

%Only write attributes if file is new:
if strtTime == 1
    %Write global attributes:
    
    info = {'Created in','Conceptual Cryosphere Hydrology Framework'; ...
            'Created on', datestr(datetime('now','Format','d-MMM-y'))};
    for ii = 1 : length(info(:,1))
        ncwriteatt(pathData,'/',info{ii,1},info{ii,2}) %'/' indicates global attribute
    end

    %Write time variable attributes:
    for ii = 1 : length(attTime(:,1))
        ncwriteatt(pathData,'time',attTime{ii,1},attTime{ii,2})
    end

    
    %Write latitude variable attributes:
    for ii = 1 : length(attLat(:,1))
        ncwriteatt(pathData,'lat',attLat{ii,1},attLat{ii,2})
    end

    
    %Write longitude variable attributes:
    for ii = 1 : length(attLon(:,1))
        ncwriteatt(pathData,'lon',attLon{ii,1},attLon{ii,2})
    end
% 
%     %Write data variable attributes:
%     %Remove '_FillValue' (it's written above) from data attributes
%     indFil = nan;
%     for ii = 1 : length(attData)
%         if ~isempty(regexpi(attData{ii,1},'fill'))
%            indFil = ii; 
%         end
%     end

%     if ~isnan(indFil)
%         attData(indFil,:) = [];
%     end


    %Write data attributes:
    for ii = 1 : length(attData(:,1))
        ncwriteatt(pathData,char(var),attData{ii,1},attData{ii,2}) %'/' indicates global attribute
    end
end

%Change NaN data to fill value:
% data( isnan(data) ) = filVal;

%%Create and write variables to NetCDF file:
%Write time vector:
ncwrite(pathData,'time', timeDays, strtTime);

ncwrite(pathData,'time_bnds', [timeDays; timeDays+1], [1, strtTime]); 

    
%Only write lat and lon variables if file is new:
if strtTime == 1
    %Latitude:
    ncwrite(pathData,'lat',lat);
    
    %Longitude:
    ncwrite(pathData,'lon',lon);
end

%Data needs to have 3 indices (if has 2, add 1)
if length(size(data)) == 2
    dataTemp = nan([1,size(data)]);
    dataTemp(1,:,:) = data;
    data = dataTemp;
end

%Convert any NaN's to filVal
data(isnan(data)) = filVal;

%Write Data:
if regexpbl(dataTyp,'uint16')
    data = uint16(data);
elseif regexpbl(dataTyp,'int16')
    data = int16(data);
elseif regexpbl(dataTyp,'single')
    data = single(data);
end

ncwrite(pathData, char(var), data, strtData, [1,1,1]);