function days_w_dates_cmpr(sData,yrs,mnths)
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


%This function checks:
%-Number of times present in metadata equal number of data grids present
%-First and last time entries for data to first and last entries of
%   'yrs' and 'mnths' to check that they match.  Returns an error if they
%   don't, but returns nothing otherwise.  
%-That the number of data elements present equals the number of dates
%   present.
%Warning: Only designed to work with time-series data, not climatologies
%(skips if type unknown or is climatology).


%Find if data is time-series or climatology:
indTyp = strcmpi(sData.attTime,'type');
if sum(sum(indTyp)) > 0
    [rowTyp, colTyp] = find(indTyp == 1);
    typeData = sData.attTime(rowTyp,colTyp + 1);
end

if sum(sum(indTyp)) > 0 && ~isempty(regexpi(typeData,'clim'))
    return
elseif sum(sum(indTyp)) == 0
    warning('days_w_dates_cmpr:noType',['The type of data present is '...
    'not known, causing the time fidelity check to be skipped.'])
    return
end
    
%Calculate number of time-series elements present in data:
if length(size(sData.data)) == 3
    nTS = length(sData.data(:,1,1));
elseif length(size(sData.data)) == 2
    nTS = 1;
end

%Check that number of time-series elements matches number of data entries:
if length(sData.time(:,1)) ~= nTS
    error('daysWDatesCmpr:TsVTime',['The number of time entries present '...
        'is ' num2str(length(sData.time(:,1))) ' but the number of data '...
        'grids present is ' num2str(nTS) '.']);
end

%Calculate start and end times from metadata:
calData =  NC_cal(sData.attTime);
[dateData, ~] = NC_time_units(sData.attTime);
dateStrtData = days_2_date(  sData.time(1), dateData, calData);
dateEndData  = days_2_date(sData.time(end), dateData, calData);

%check that first dates match:
if ~isequal(dateStrtData(1:2),[yrs(1),mnths(1)]) 
    error('daysWDatesCmpr:fTsLd',['The time metadata for the '...
        ' first entry is ' num2str(dateStrtData(2)) '/' ...
        num2str(dateStrtData(1)) ' but the first time-series element '...
        'requested is ' num2str(mnths(1)) '/' num2str(yrs(1)) '.']);
end

if ~isequal(dateEndData(1:2),[yrs(2),mnths(end)])
    error('daysWDatesCmpr:lTsLd',['The time metadata for the '...
        ' last entry is ' num2str(dateEndData(2)) '/' ...
        num2str(dateEndData(1)) ' but the first time-series element '...
        'requested is ' num2str(mnths(1)) '/' num2str(yrs(1)) '.']);
end

%Calculate number of requested data entries:
nReq = (yrs(2) - yrs(1) + 1)*length(mnths);

if nTS ~= nReq
    error('daysWDatesCmpr:TsVReq',['The number of grids that should '...
        'have been loaded is ' num2str(nReq) ' but ' num2str(nTS) ...
        ' were actually loaded.']);
end

end