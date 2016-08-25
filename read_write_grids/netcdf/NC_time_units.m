function [gcmRef, gcmUnits] = NC_time_units(attTime)
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


%Find reference date used in GCM, following standard CMIP5 metadata 
%convention
%Input:
%'attTime' = time metadata following conventions used in CMIP5 GCMs

%Output:
%'gcmRef' = [year, month, day] using in GCM as reference
%'gcmUnits' = method used in GCM to define time entry relative to reference

tUnitInd = strcmpi(attTime,'units');
[tUntRow, tCol] = find(tUnitInd == 1);
if isempty(tUntRow) || isempty(tCol)
    gcmRef = 'Unknown';
    gcmUnits = 'Unknown';
    return
end

tUntInfo = attTime{tUntRow, tCol+1};
strtDateInd = regexpi(tUntInfo,'-');

if ~isempty(strtDateInd)
    indSpace = regexpi(tUntInfo,' ');
    indEnd = find(indSpace > strtDateInd(end), 1, 'first');
    if ~isempty(indEnd)
        indEnd = indSpace(indEnd);
    else
       indEnd = numel(tUntInfo);
    end
elseif isempty(strtDateInd) && regexpbl(tUntInfo,{'ka','bp'},'and')
    gcmUnits = tUntInfo;
    gcmRef  = [2000, 1, 1];
    disp(['Reference units of ka years before present are used in ' ...
        'GCM dataset.  Reference time of ' num2str(gcmRef(3)) '/' ...
        num2str(gcmRef(2)) '/' num2str(gcmRef(1)) ' are being chosen to represent present.']);
    
    return
else
    error('NC_time_units:uknownUnit',[char(39) tUntInfo char(39) ...
        ' are unknown units and have not been programmed for.']);
end

gcmUnitsInd = regexpi(tUntInfo,'\d');
if ~isempty(gcmUnitsInd)
    gcmUnits = tUntInfo(1:gcmUnitsInd(1)-2);
else
    gcmUnits = tUntInfo;
end

%Find reference date used:
if ~isempty(strtDateInd)
    if numel(strtDateInd) == 1
        gcmRef = [  str2double(tUntInfo(strtDateInd(1)-4:strtDateInd(1)-1)), ...
                    str2double(tUntInfo(strtDateInd(1)+1:end)), 1]; 
    else
        gcmRef = [  str2double(tUntInfo(strtDateInd(1)-4:strtDateInd(1)-1)), ...
                    str2double(tUntInfo(strtDateInd(1)+1:strtDateInd(2)-1)), ...
                    str2double(tUntInfo(strtDateInd(2)+1:indEnd))];
    end
end

if regexpbl(gcmUnits,'hour') && numel(gcmRef) == 3
   gcmRef = [gcmRef, 1]; 
end