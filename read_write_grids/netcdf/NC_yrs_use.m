function [yrsLoad, nYrsUse, gcmTimeUseInd] = NC_yrs_use(sMeta, gcmStrtDate, gcmEndDate, varargin)
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

%Choose which years to load based upon 'varargin'
[~, yrsLoad, mnthsLoad] = output_type(varargin{1},sMeta);

%Find starting time indice to use:
if yrsLoad(1) >= gcmStrtDate(1) && gcmEndDate(1) >= yrsLoad(1) %Desired start year contained within GCM
    gcmTimeUseInd(1) = 12*(yrsLoad(1) - gcmStrtDate(1)) + 1;
elseif gcmStrtDate(1) >= yrsLoad(1) && gcmEndDate(1) >= yrsLoad(1) %GCM starts after start year but before end year
    if ~ischar(varargin{1}) || isempty(regexpi(varargin{1},'clim'))
        warning('gcm_load:startYr',['The downscaled time-series '...
            'will begin in ' num2str(gcmStrtDate(1)) ' instead of ' num2str(yrsLoad(1)) ' because '...
            'this is when the GCM data begins.']);
    else
        warning('gcm_load:startYr',['The GCM data used in the '...
            'climatology production begins in ' ...
            num2str(gcmStrtDate(1)) ' instead of ' ...
            num2str(yrsLoad(1)) '.  The historical climatology ' ...
            'used for quantile-based bias correction will ' ...
            'therefore be cropped to the available GCM data']);
    end
    yrsLoad(1) = gcmStrtDate(1);

    if gcmStrtDate(2) <= mnthsLoad(1)
        gcmTimeUseInd(1) = gcmStrtDate(2);
    else
        warning('gcm_load:startMnth',['The downscaled time-series '...
            'will begin in ' num2str(gcmStrtDate(1)) '-' ...
            num2str(gcmStrtDate(2)) ' instead of ' ...
            num2str(yrsLoad(1)) '/' num2str(mnthsLoad(1)) ...
            ' because this is when the GCM data begins.']);
    end

elseif yrsLoad(1) > gcmEndDate(1) %GCM ends before desired range of years begins
    error('gcm_load:noTimeOverlap',['The time-series data ends in ' ...
        num2str(gcmEndDate(1)) ' and the requested year to begin '...
        'downscaling is ' num2str(yrsLoad(1)) ...
        ', which is problematic...']);
else
    error('gcm_load:unknownTimeCase',['No logical structure has '...
        'been written to deal with this logically impossible case.']);
end

%Find ending time indice to use:
%These conditionals ensure that GCM time series includes all requested
%elements.  If not, it alters the time-series being downscaled 
if gcmEndDate(1) > yrsLoad(2) || (gcmEndDate(1) == yrsLoad(2) && gcmEndDate(2) >= mnthsLoad(end))
%Skip if statement in this case!

%If the CGM and requested downscaled time-series end in the same year & 
%not all mnthsLoad of that year are in GCM time-series:   
elseif gcmEndDate(1) == yrsLoad(2) && gcmEndDate(2) <= mnthsLoad(end) 
    yrsLoad(2) = yrsLoad(2) - 1;

    warning('gcm_load:endMnthTrunc',['The downscaled time-series '...
        'will end in ' num2str(yrsLoad(2)) ' because the GCM '...
        'time-series does not have all of the month data for the '...
        'final requested year.']);
%If the GCM ends before the requested year of the downscaled time-series
elseif gcmEndDate(1) < yrsLoad(2)
    yrsLoad(2) = gcmEndDate(1);

    if gcmEndDate(2) < mnthsLoad(end)
        yrsLoad(2) = yrsLoad(2) - 1;

        warning('gcm_load:EndYrMnth',['The downscaled time-series '...
            'will end in ' num2str(gcmEndDate(1)) '-' ...
            num2str(gcmEndDate(2)) ...
            ' instead of ' num2str(yrsLoad(2)) ' because the GCM does not extend through December.']);
    end

    warning('gcm_load:EndYr',['The output downscaled time-series ' ...
        'will end in the year ' num2str(gcmEndDate(1)) ...
        ' because this is when the GCM data ends.']);
else
    error('gcm_load:noYrSelected',['The GCM end date is ' ...
        num2str(gcmEndDate(2)) '/' num2str(gcmEndDate(1)) ' and the '...
        'final desired time-step to downscaled is ' ...
        num2str(mnthsLoad(end)) '/' num2str(yrsLoad(2)) ...
        '.  This case has not been coded for.']);
end

if yrsLoad(1) > yrsLoad(2) 
    warning('gcm_load:invalidEndYr', ['The requested ending year of ' ...
        'GCM data to load is less than the starting year.  Therefore, ' ...
        'the ending year will be changed from ' num2str(yrsLoad(2)) ...
        ' to ' num2str(yrsLoad(1)) ', which is the beginning year.']);
    yrsLoad(2) = yrsLoad(1);
end

%Define ending time indice:
gcmTimeUseInd(2) = 12*(yrsLoad(2) - yrsLoad(1) + 1) ...
    - 1 + gcmTimeUseInd(1);
%OLD VERSION:
% gcmTimeUseInd(2) = 12*(yrsLoad(2) - yrsLoad(1)) + mnthsLoad(end) ...
%     - 1 + gcmTimeUseInd(1);
%Total number of years selected:
nYrsUse = yrsLoad(2) - yrsLoad(1) + 1;

end