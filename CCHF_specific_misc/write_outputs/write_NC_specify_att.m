function write_NC_specify_att(pathData, sData, varargin)
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
    %.var. = all variables to write (must have lat, lon, time, and corresponding georeferenced variable
    %.att. = all attributes to write
    
  
fields = fieldnames(sData.var);
varAll = fields;
fields(strcmpi(fields,'lon')) = [];
    fields(strcmpi(fields,'longitude')) = [];
    fields(strcmpi(fields,'lat')) = [];
    fields(strcmpi(fields,'latitude')) = [];
    fields(strcmpi(fields,'time')) = [];
    fields(strcmpi(fields,'time_bnds')) = [];
if numel(fields) == 1
    var = fields{1};
else
    error('write_NC_specify_att:autoVarFail','Cannot auto-detect variable to load.');
end

varAll(strcmpi(varAll,'time')) = [];
varAll(strcmpi(varAll,'time_bnds')) = [];
varAll(strcmpi(varAll,var)) = [];  

nmAtt = fieldnames(sData.att);
    

%Find if filVal exists:
filVal = find_att(sData.att.(var),'FillValue', 'no_warning');
if ~isempty(filVal)
    indFil = strcmpi(sData.att.(var),'FillValue');
    sData.att.(var)(indFil,:) = [];
end
if isempty(filVal)
    filVal = find_att(sData.att.(var),'_FillValue', 'no_warning');
    if ~isempty(filVal)
        indFil = strcmpi(sData.att.(var),'_FillValue');
        sData.att.(var)(indFil,:) = [];
    end
end



%DATA ACTUALLY NEED TO HAVE THREE DIMENIONS
% %Ensure data does not have extra dimensions:
% data = squeeze(data);

%Select data type to use (default is single precision):
if iscell(varargin) && ~isempty(varargin)
    prec = varargin{1};
    if prec == 0
        dataTyp = 'int16';
        if ~isempty(filVal)
            filVal = intmax(dataTyp); %Missing data value (all NaN values go to this):
        end
    elseif prec > 0 && prec < 5
        dataTyp = 'single';
        if ~isempty(filVal)
            filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
        end
    else
        dataTyp = 'double';
        if ~isempty(filVal)
            filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
        end
    end
elseif ~isempty(regexpi(char(var),{'pre','pr'}))
    dataTyp = 'uint16';
    if ~isempty(filVal)
        filVal = intmax(dataTyp); %Missing data value (all NaN values go to this):
    end
else
    dataTyp = 'single';
    if ~isempty(filVal)
        filVal = flintmax(dataTyp); %Missing data value (all NaN values go to this):
    end
end    


%keyboard

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
        {'time' Inf 'lat' length(sData.var.lat(:)) 'lon' length(sData.var.lon(:))}, ...
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
    
    nccreate(pathData,'lat', 'Dimensions',{'lat' length(sData.var.lat(:))}, ...
        'Format',strFormat, 'Datatype', 'double');
    nccreate(pathData,'lon','Dimensions',{'lon' length(sData.var.lon(:))}, ...
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

% attData = sData.att.(var);
% attTime = sData.att.time;
% 
% attLon = sData.att.lon;
% 
% attLat = sData.att.lat;
  
% 
% timeDays = days_since([1900,1,1], timeVec, 'gregorian');
% 
% %Check if time elements are nan:
% if sum(isnan(timeDays)) > 0
%    warning('write_NC:nanTime',[num2str(sum(isnan(vecTime))) ' NaN time '...
%        'values have been detected.  This may cause the program to crash.']); 
% end

%Only write attributes if file is new:
if strtTime == 1
    %Write global attributes:
    
    info = {'Created in','Conceptual Cryosphere Hydrology Framework'; ...
            'Created on', datestr(datetime('now','Format','d-MMM-y'))};
    for ii = 1 : length(info(:,1))
        ncwriteatt(pathData,'/',info{ii,1},info{ii,2}) %'/' indicates global attribute
    end

    
    for jj = 1 : numel(nmAtt)
        if ~isempty(sData.att.(nmAtt{jj}))
            for kk = 1 : numel(sData.att.(nmAtt{jj})(:,1))
                ncwriteatt(pathData,nmAtt{jj}, sData.att.(nmAtt{jj}){kk,1}, sData.att.(nmAtt{jj}){kk,2})
            end
        end
    end
    
    
%     %Write time variable attributes:
%     for ii = 1 : length(attTime(:,1))
%         ncwriteatt(pathData,'time',attTime{ii,1},attTime{ii,2})
%     end
% 
%     
%     %Write latitude variable attributes:
%     for ii = 1 : length(attLat(:,1))
%         ncwriteatt(pathData,'lat',attLat{ii,1},attLat{ii,2})
%     end
% 
%     
%     %Write longitude variable attributes:
%     for ii = 1 : length(attLon(:,1))
%         ncwriteatt(pathData,'lon',attLon{ii,1},attLon{ii,2})
%     end
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


%     %Write data attributes:
%     for ii = 1 : length(attData(:,1))
%         ncwriteatt(pathData,char(var),attData{ii,1},attData{ii,2}) %'/' indicates global attribute
%     end
end

%Change NaN data to fill value:
sData.var.(var)( isnan(sData.var.(var)) ) = filVal;

%Data needs to have 3 indices (if has 2, add 1)
if length(size(sData.var.(var))) == 2
    dataTemp = nan([1,size(sData.var.(var))]);
    dataTemp(1,:,:) = sData.var.(var);
    sData.var.(var) = dataTemp;
end

%Write Data:
if regexpbl(dataTyp,'int16')
    sData.var.(var) = int16(sData.var.(var));
elseif regexpbl(dataTyp,'single')
    sData.var.(var) = single(sData.var.(var));
end

ncwrite(pathData, char(var), sData.var.(var), strtData, [1,1,1]);


%%Create and write variables to NetCDF file:
%Write time vector:
ncwrite(pathData,'time', sData.var.time, strtTime);

if isfield(sData.var,'time_bnds')
    %Write all fields:
    ncwrite(pathData,'time_bnds', [sData.var.time; sData.var.time+1], [1, strtTime]);
end

for ii = 1 : numel(varAll)
    ncwrite(pathData,varAll{ii}, sData.var.(varAll{ii}));
end


