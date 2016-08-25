function [gcmTime, gcmDateVec, attTime] = NC_avail_time(gcmTime, attTime)
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


%Finds start and end dates for NetCDF climate data file.
%'strtDate' = [year; month]


%Read time units from GCM file:
[gcmRef, tUnits] = NC_time_units(attTime);


if regexpbl(gcmRef,'unknown')
    gcmDateVec = nan(1,3);
    return
end
%Find calendar info:
strCal =  NC_cal(attTime);

if regexpbl(strCal, 'unknown')
    strCal = 'gregorian';
    warning('NC_avail_time:unknownCal','The calendar is unknown, but is being set to Gregorian.');
end
%Find time offset from the reference date:
if regexpbl(tUnits, 'days since')
    gcmDateVec = days_2_date(gcmTime, gcmRef, strCal);
elseif regexpbl(tUnits, 'Ka BP') && regexpbl(strCal,'noleap')
    gcmDateVec = time_kaBP_2_standard(gcmTime,gcmRef,strCal);
elseif regexpbl(tUnits, 'hours since')
    if numel(gcmRef) == 3
        gcmRef = [gcmRef, 0];
    end
    gcmDateVec = days_2_date(gcmTime/24, gcmRef, strCal);
else
    disp(['The GCM' char(39) 's time units are' tUntInfo '.']);
    error('NC_time:unitsUnknown','No case has been written to deal with the current time units.');
end

%Check ending date:
    %one reason to check is because different calendars can be used (e.g. 
    %Gregorian vs. no leap year)
%Calculate total days elapsed over time period included in the GCM:
daysCalcGcm = days_since(gcmDateVec(1,:), gcmDateVec(end,:), strCal);

%Calculate actual days in GCM:
if regexpbl(tUnits, 'days since')
    daysActGcm = gcmTime(end) - gcmTime(1);
elseif regexpbl(tUnits, 'hours since')
    daysActGcm = (gcmTime(end) - gcmTime(1))/24;
end

%Check correspondance between actual time and calculated time
    %Error if difference in calcualted and actual days is >= half a month. This should be changed for daily or subdaily data. 
if abs(daysCalcGcm - daysActGcm) > 14 && isempty(regexpi(tUnits, 'Ka BP'))
    disp(['The calculated number of days in the GCM is ' num2str(daysCalcGcm) '.']);
    disp(['The actual number of days in the GCM is ' num2str(daysActGcm)]);
    error('NC_time:gcmDayError','The number of days calculated to be in the GCM is not equal to the number actually in it.');
end


end