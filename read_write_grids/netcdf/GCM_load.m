function [gcmTs, gcmDataUnit, gcmLat, gcmLon, yrsOutTemp] = GCM_load(gcmPath, sMeta, clim)
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


%GCM time-series is output in mm of precipitation and Celsius (provided that loaded file is in pre flux or Kelvin)

%Format for 'clim' option:

% %Test selection parameters:
% metVar = 'pre';
% mnthsLoad = (1:12);
% yrsLoad = [2001; 2030]
% 
% strucHdr.crd = [-120,-100,45,60];
% %IPCC AR4 (GRB)

metVar = sMeta.currVar;

%Load GCM NetCDF of GRB file:
ds = ncdataset(gcmPath);

if ~exist('clim','var') || ~iscell(clim)
    clim = cell(3,1);
end

%If years given as third argument, use them, else use years from 'strucHdr'
if isempty(clim{2})
    yrsLoad = sMeta.yrsOut;
elseif ~isempty(clim{2})
    yrsLoad = clim{2};
else
    error('gcm_load:yrsLoad','No year vector detected designating which GCM data to load.');
end

%Use to check if requested downscaled years are altered later in this
%script:
yrsOutTemp = yrsLoad;


%%Find indices for time of interest in NetCDF or GRB file:
[gcmStrtDate, gcmEndDate] = NC_time(ds);

%Compares user requested time-series elements to those available; amends
%'yrsLoad' as necessary; selects indices for all 12 months
[yrsLoad, nYrsUse, gcmTimeUseInd] = NC_yrs_use(sMeta, gcmStrtDate, gcmEndDate, clim);

%%Find which latitude indices to use:
[gcmLat, gcmLatUseInd] = NC_lat_use(ds, sMeta);

%%Find which longitude indices to use:
[gcmLon, gcmLonUseInd] = NC_lon_use(ds, sMeta);

%%Load and crop GCM data:
[gcmTs, gcmDataUnit] = NC_data_load(ds, sMeta, metVar, gcmTimeUseInd, gcmLatUseInd, gcmLonUseInd, nYrsUse, clim);

%Check that years to load have not been altered through the script.  If
%this has happened, warn user
if sum(yrsOutTemp == yrsLoad) ~= 2
    warning('gcm_load:yrsChanged',['The original range of years ' ...
        'selected to be downscaled was ' num2str(yrsOutTemp(1)) '->' ...
        num2str(yrsOutTemp(2)) ' but has been changed to ' num2str(yrsLoad(1)) ...
        '->' num2str(yrsLoad(2)) ' because of time-series elements present ' ...
        'in the GCM data.']);
    yrsOutTemp = yrsLoad;    
end

%Check that number of time-series elements loaded matches those requested
if length(gcmTs(:,1,1)) ~= yrsLoad*length(sMeta.mnths)
    error('gcmload:wrongts','An incorrect number of time-series elements were loaded.');
end

%If option chosen, make climatology
if ~isempty( regexpi( char(clim{1}),'clim') )
    if nYrsUse == length(gcmTs(:,1,1))
        gcmClimTemp = zeros(size(gcmTs(1,:,:)));
        for ii = 1 : nYrsUse
            gcmClimTemp = gcmClimTemp + gcmTs(ii,:,:);
        end
        
        gcmTs = gcmClimTemp / nYrsUse;
    else
        error('gcm_load:climMnth',['The function is in climatology '...
            'mode, although there is an error with the number of '...
            'time-series elements present.']);
    end
gcmTs = squeeze(gcmTs);
end


end




%FOR TESTING:
%     yrsLoad
%     gcmStrtDate
%     gcmEndDate
%     gcmTimeUseInd
%     gcmLat
%     gcmLon
%     (yrsLoad(2) - yrsLoad(1) + 1)*12
