function [ind, closest] = geodata_time_ind(sGeo, currTime)
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

%Function takes geodata structure array and time interval and returns
%indice of data corresponding to indice of requested time.
%%INPUTS: 
%sGeo = structure array in Mosier geodata format
%currTime = time vector in format [year, month, day] (can also have hour
%and minute)
%%OUTPUTS:
%ind = indice of time element in sGeo
%closest = difference between time value and requested time (in days)

if ~isfield(sGeo,'attTime')
    ind = nan;
    closest = nan;
    return
else
    [tRef, units] = NC_time_units(sGeo.attTime);
    strCal =  NC_cal(sGeo.attTime);

    %Create date string [year, month, day]
    if length(tRef) == 2 %If no day of the month provided, add:
       tRef = [tRef, 0.5*eomday(tRef(1), tRef(2))]; 
    elseif isnan(tRef(3))
        tRef(3) = 0.5*eomday(tRef(1), tRef(2));
    end
    if length(currTime) == 2
       currTime = [currTime, 0.5*eomday(currTime(1), currTime(2))]; 
    end

    if regexpbl(units, {'days','since'},'and')
        nDays = days_since(tRef, currTime, strCal);
    elseif regexpbl(units, 'Unknown')
        nDays = NaN;
    else
        disp(['The geodata' char(39) 's time units are ' units '.']);
        error('geodata_curr_time:unitsUnknown','No case has been written to deal with the current time units.')
    end
    diffGeo = abs(sGeo.time - nDays);
    [closest, ind] = min(diffGeo);
    if abs(closest) > 31
        warning(['The difference between the current time and the ' ...
            'current time is ' num2str(closest)]);
    end
end