function [gcmTime, indGcmTime, attTime] = NC_monthly_time_use(ds, sMeta, varargin)
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


%Function compares time-series elements in NetCDf file with desired dates
%and, if necessary, adapts time-series elements to load.



%Find start and ending dates for NetCDF data:
[gcmTime, gcmTimeVec, attTime] = NC_avail_time(ds);
%Format of gcmTimeVec = [month, year, day]

%Choose which years to load based upon 'varargin'
[outTyp, yrsLd, mnthsLd] = output_type(varargin{1},sMeta);
if isempty(outTyp) || ~isnumeric(outTyp)
    yrsLd = (yrsLd(1):yrsLd(2));
else
    yrsLd = outTyp(:,1);
end

mnthsLd = unique(mnthsLd);

if numel(yrsLd(1,:)) > numel(yrsLd(:,1))
   yrsLd = yrsLd'; 
end

if numel(mnthsLd(1,:)) > numel(mnthsLd(:,1))
   mnthsLd = mnthsLd'; 
end

%Create date vector of monthly time-series elements to load:
timeLdVec = nan(numel(yrsLd)*numel(mnthsLd),3,'single');
if numel(timeLdVec(:,1)) ~= 0
    for ii = 1 : length(yrsLd)
        timeLdVec((ii-1)*numel(mnthsLd)+1:ii*numel(mnthsLd),:) = [yrsLd(ii)*ones(numel(mnthsLd),1), mnthsLd, round(0.5*eomday(yrsLd(ii),mnthsLd))];
        %day = 15 because assumes monthly data
    end
end

if isfield(sMeta,'days')
    for ii = numel(timeLdVec(:,1)) : -1 : 1
        currDays = (1 : eomday(timeLdVec(ii,1), timeLdVec(ii,2)))' + 0.5;
        if ii == numel(timeLdVec(:,1))
            timeLdVec = [timeLdVec(1:ii-1,:) ; [timeLdVec(ii,1)*ones(numel(currDays),1), timeLdVec(ii,2)*ones(numel(currDays),1), currDays]];
        elseif ii == 1
            timeLdVec = [[timeLdVec(ii,1)*ones(numel(currDays),1), timeLdVec(ii,2)*ones(numel(currDays),1), currDays]; timeLdVec(ii+1:end,:)]; 
        else
           timeLdVec = [timeLdVec(1:ii-1,:); [timeLdVec(ii,1)*ones(numel(currDays),1), timeLdVec(ii,2)*ones(numel(currDays),1), currDays]; timeLdVec(ii+1:end,:)]; 
        end
    end
end



%Assumes monthly data because only compares year and month:
if isfield(sMeta,'days')
    [dateBoth,indGcmTime,indLd] = intersect(gcmTimeVec, timeLdVec,'rows');
else
    [dateBoth,indGcmTime,indLd] = intersect(gcmTimeVec(:,1:2),timeLdVec(:,1:2),'rows');
end

%Display warning if some data not present
if length(indLd) < length(timeLdVec(:,1))
%     if isfield(sMeta,'days')
%         timeXAvail = setxor(dateBoth,timeLdVec(:,1:2),'rows');
%     else
%          timeXAvail = setxor(dateBoth,timeLdVec,'rows');
%     end
    
    warning('NC_monthly_time_use:notAvail',['GCM data not available '...
        'for some requested dates.' char(10) ...
        'This may cause program to crash.']);
end

%Keep only time entries corresponding to loaded data:
gcmTime = gcmTime(indGcmTime);

end