function [sData] = NC_load_v2(path,sMeta,varargin)
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


%Designed to load data from NetCDF file using CMIP5 conventions.

%Input:
%'path' refers to a geo-referenced NetCDF or ESRI formatted ASCII file
%'sMeta' is a structure array with the following fields (among others):
    %sMeta.currVar = variable to load
    %sMeta.lRHist = 'GPCC', 'CRU', or 'Willmott'
    %sMeta.yrsOut = [first year of output data, last year of output data]
    %sMeta.yrsClim = [first year of climatology, last year of climatology]
    %sMeta.mnths = list of months requested
    %sMeta.crd = [lonW, lonE, latS, latN]
    %sMeta.currTime = [year, month]
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


path = char(path);

sNC = ncinfo(path);
nVar = numel(sNC.Variables);
nmVar = cell(nVar,1);
for ii = 1 : nVar
    nmVar{ii} = sNC.Variables(ii).Name;
end

%Get longitude:
if any(strcmpi(nmVar,'longitude'))
    ncLon = ncread(path,'longitude');
    sData.attLon = ncinfo(path,'longitude');
    sData.attLon = squeeze(struct2cell(sData.attLon.Attributes))';
elseif any(strcmpi(nmVar,'lon'))
    ncLon = ncread(path,'lon');
    sData.attLon = ncinfo(path,'lon');
    sData.attLon = squeeze(struct2cell(sData.attLon.Attributes))';
else
    error('NC_load:noLonVar','The name of the longitude variable could not be located.');
end
if numel(ncLon(:,1)) > numel(ncLon(1,:))
    ncLon = ncLon';    
end

%Get Latitude:
if any(strcmpi(nmVar,'latitude'))
    ncLat = ncread(path,'latitude');
    sData.attLat = ncinfo(path,'latitude');
    sData.attLat = squeeze(struct2cell(sData.attLat.Attributes))';
elseif any(strcmpi(nmVar,'lat'))
    ncLat = ncread(path,'lat');
    sData.attLat = ncinfo(path,'lat');
    sData.attLat = squeeze(struct2cell(sData.attLat.Attributes))';
else
    error('NC_load:noLatVar','The name of the latitude variable could not be located.');
end
if numel(ncLat(:,1)) < numel(ncLat(1,:))
    ncLat = ncLat';    
end

%Get time:
if any(strcmpi(nmVar,'time'))
    ncTime = ncread(path,'time');
    sData.attTime = ncinfo(path,'time');
    sData.attTime = squeeze(struct2cell(sData.attTime.Attributes))';
else
    error('NC_load:noTimeVar','The name of the time variable could not be located.');
end

%Get variable name for data (don't read):
if isfield(sMeta,'currVar')
    metVar = char(sMeta.currVar);
    %%Load data corresponding to time, lat, and lon of interest from file:
    if ~isempty(regexpi(metVar,'pre'))
        dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'pre'));

        if sum(dataInd) == 0
            dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'pr'));
        end
    elseif ~isempty(regexpi(metVar,'tmp'))   
        dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'tmp'));

        if sum(dataInd) == 0
            dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'tas'));
        end
    elseif ~isempty(regexpi(metVar,'tmn'))
        dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'tmn'));
        if sum(dataInd) == 0
            dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'tasmin'));
        end
    elseif ~isempty(regexpi(metVar,'tmx'))
        dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'tmx'));
        if sum(dataInd) == 0
            dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,'tasmax'));
        end
    else
        dataInd = cellfun(@(x) ~isempty(x), regexpi(nmVar,sMeta.currVar));
    end

    if sum(dataInd) > 1 %This prevents the wrong selection from being 
                        %automatically made and prompts the user to select the 
                        %correct variable
        dataInd = 0;
    end
else
    dataInd = 0;
    metVar = 'Unknown';
end

%Try to find variable if (1) not present in sMeta or if not found above.
if sum(dataInd) == 0
    nmVarTemp = nmVar;
    nmVarTemp(strcmpi(nmVarTemp,'lon')) = [];
    nmVarTemp(strcmpi(nmVarTemp,'longitude')) = [];
    nmVarTemp(strcmpi(nmVarTemp,'lat')) = [];
    nmVarTemp(strcmpi(nmVarTemp,'latitude')) = [];
    nmVarTemp(strcmpi(nmVarTemp,'time')) = [];
    nmVarTemp(strcmpi(nmVarTemp,'time_bnds')) = [];
    
    if numel(nmVarTemp) == 1
        dataInd = find(strcmpi(nmVar, nmVarTemp{1}) == 1);
    else
        %Display GUI for user to find variable if the automated system was unable to find it. 
        dataInd = NC_man_sel(dataInd, nmVar, metVar);
    end
end
if sum(dataInd) == 0
    error('gcm_load:metVarUnknown',[char(39) sMeta.currVar char(39)...
    'is either not present in the NetCDF data or is an unknown '...
    'name.  Can code for this.']);
end
nmData = char(nmVar(dataInd));


%Selects full years! 
%- Must later subtract out undesired months.
%Choose years and months to load based upon 'varargin':
[outTyp, ~, ~] = output_type(varargin{1}, sMeta);

%Find time indices to use:
[sData.time, timeInd] = NC_time_use_v2(ncTime, sData.attTime, sMeta, outTyp);

if all(isnan(timeInd))
    timeGcm = ncread(path,'time');
    timeInd = 1 : numel(timeGcm);
end


%Ensure input coordinates in proper order:
if isfield(sMeta,'crd')
    if sMeta.crd(3) > sMeta.crd(4)
       sMeta.crd(3:4) = fliplr(sMeta.crd(3:4));
    end
    if sMeta.crd(1) > sMeta.crd(2)
       sMeta.crd(1:2) = fliplr(sMeta.crd(1:2));
    end
end


%%Find latitude indices to use:
[sData.lat, latInd] = NC_lat_use_v3(ncLat, sData.attLat, sMeta);

%%Find longitude indices to use:
[sData.lon, lonInd] = NC_lon_use_v3(ncLon, sData.attLon, sMeta);


%%Load GCM data:
[sData.data, sData.attData, sData.info, sData.var] ...
    = NC_data_use_v2(path, nmData, timeInd, latInd, lonInd);


%Make .lon column vector and .lat row vector:
if length(sData.lon(:,1)) > 1 && length(sData.lon(1,:)) == 1 
    sData.lon = sData.lon';
end
if length(sData.lat(1,:)) > 1 && length(sData.lat(:,1)) == 1 
    sData.lat = sData.lat';
end


%Sort data and ensure elements correspond to map directions (i.e. top =
%North and right = East):
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


%Check for 'units' paramter in metadata
dataUnitInd = strcmpi(sData.attData,'units');
%Initialize units variable:
gcmDataUnit = 'unknown';
%Locate unit metadata entry:
if sum(sum(dataUnitInd)) > 0
    [dataUnitRow, dataUnitCol] = find(dataUnitInd == 1);
    gcmDataUnit = sData.attData{dataUnitRow, dataUnitCol+1};
else
    error('NC_load:unkownUnits',['The units of ' char(39) path char(39) ...
        ' were not located due to the NetCDF format.  Please examine ' ...
        'this case and code for format.']);
end


%Find units and calendar:
[gcmRef, gcmTimeUnits] = NC_time_units(sData.attTime);
calGcm =  NC_cal(sData.attTime);
%Surmize timestep in data:
timeTest = ncread(path,'time');
tStep = t_res_geodata(timeTest,gcmTimeUnits);
if regexpbl(tStep,'unknown')
    tStep = find_att(sData.attTime,'time_step','no_warning');
    if isempty(tStep)
        error('NC_load_v2:tStep','Data have unknown time step.');
    end
end

if ~any(strcmpi(sData.attTime,'time_step'))
    sData.attTime(end+1,:) = {'time_step', tStep};
end


% if isempty(sData.time)
%    keyboard 
% end

%If 'time_bounds' present, transform 'sData.time' to be mid-point of bounds. 
indTBnds = find(strcmpi(nmVar, 'time_bnds'));
if ~isempty(indTBnds) && ~isempty(sData.time)
    tBnds = double(ncread(path,nmVar{indTBnds}));
    sData.('time_bnds') = tBnds;
    
    %'time_bnds' does not have attributes.
%    	sData.attTime_bnds = ncinfo(path,'time_bnds');
%     sData.attTime_bnds = squeeze(struct2cell(sData.attTime_bnds.Attributes))';
%     
%     sData.('attTime_bnds') = 
    
    if numel(tBnds(1,:)) == 2 %Transpose to put in correct orientation for 'diff()'
        tBnds = tBnds';
    end
    dTime = diff(tBnds)';
    if sData.time(1) == tBnds(1,1) 
        sData.time = sData.time + 0.5*dTime;
    elseif sData.time(1) == tBnds(2,1) 
        sData.time = sData.time - 0.5*dTime;
    end
else
	if regexpbl(gcmTimeUnits,'hour')
        %Change hr = 1 : hr = 0.5 and add tBounds
        dHr = diff(sData.time);
        if numel(dHr) ~= 0
           if all(dHr == dHr(1))
               sData.('time_bnds') = [sData.time' - dHr(1)*ones(1,numel(sData.time)); sData.time'];
               sData.time = sData.time - 0.5*dHr(1);
           else
                warning('NC_load_v2:unknownTOffset','Entering keyboard to diagnose situation.')
                keyboard
           end
           
            if ~regexpbl(gcmTimeUnits,{'day','month'})
                warning('NC_load_v2:unknownTOffset','Entering keyboard to diagnose situation.')
                keyboard
            end
        end
	end
end





%This finds if the units of precipitation are likely flux and converts to
%total monthly mm of precipitation:
if ~isempty(regexpi(gcmDataUnit,'kg m-2 s-1')) || ~isempty(regexpi(gcmDataUnit,'m/s')) 
    scalePreDay = 3600*24;
    
    %'kg m-2 s-1' to 'mm' OR 'm/s' to 'm'
    if regexpbl(tStep,'day')
        sData.data = scalePreDay*sData.data;
    elseif regexpbl(tStep,'hour')
%             if numel(gcmRef) == 3
%                gcmRef = [gcmRef, 1]; 
%             end

            if ~isempty(indTBnds)
                dHr = diff(tBnds);
            else
                dHr = diff(sData.time);
            end
            if regexpbl(gcmTimeUnits,{'days','since'},'and')
                dHr = diff(sData.time)*24;
            end
            
            
            if all(dHr == dHr(1))
                dHr = dHr(1);
                sData.data = sData.data*3600*dHr;
            else
                if numel(dHr) == numel(sData.time)
                    for ii = 1 : length(sData.time)
                        sData.data(ii,:,:) = sData.data(ii,:,:)*3600*dHr(ii);
                    end
                else
                    error('NC_load:unevenDiffHr',['The number of hours in '...
                        'each time step varies, which has not been coded for.']);
                end
            end
    elseif regexpbl(tStep,'month')
        if regexpbl(gcmTimeUnits,{'days','since'},'and')
            dateOut = days_2_date(sData.time, gcmRef, calGcm);
        elseif regexpbl(gcmTimeUnits,'ka BP')
            dateOut = time_kaBP_2_standard(sData.time,gcmRef,calGcm);
        else
            error('NC_load:timeUnit',['Not coded for time units: ' gcmTimeUnits]);
        end

        dateOut(:,3) = days_in_month(dateOut,calGcm);

        for ii = 1 : length(dateOut(:,1))
            sData.data(ii,:,:) = sData.data(ii,:,:)*scalePreDay*dateOut(ii,3);
        end
    else
        error('NC_load:timeUnit',['Not coded for time units: ' gcmTimeUnits]);
    end
    
    
    %If 'm/s' to 'm', convert to mm:
    if regexpbl(gcmDataUnit,'m/s')
        sData.data = 1000*sData.data;
    end

    %Change units label:
    sData.attData{dataUnitRow, dataUnitCol+1} = 'mm';
elseif strcmpi(sData.var,'pr') && ~regexpbl(sData.attData{dataUnitRow, dataUnitCol+1},'mm')
    error('NC_data_use:PreUnits',['The stated units of this ' ...
        'netCDF file are ' gcmDataUnit ' but the required units '...
        'are mm.  Write a clause to convert this to a standard unit.']);
end

%Convert from Kelvin to Celsius:
if ~isempty(regexpi(gcmDataUnit,'kelvin')) || strcmpi(gcmDataUnit,'k') %GCM data seems to use 'K' instead of writing out 'kelvin'
    %Convert from Kelvin to Celsius:
    sData.data = sData.data - 273.15;
    %Change units label:
    sData.attData{dataUnitRow, dataUnitCol+1} = 'Celsius';
end

end