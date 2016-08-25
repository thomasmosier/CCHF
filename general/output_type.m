function [outTyp, yrs, mnths] = output_type(strTyp,sMeta)
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


%Input:
%'strTyp' = 'ts' (loads only current month but all years), 
    %'c' (loads current mnth and climatology years), 'all' (loads all months
    %and all years), 'none' indicating that the data has no temporal reference, 
    %or [year, month] numeric vector indicating one time-series element to load.
%'sMeta' = structure array with fields 'yrsOut' and 'currTime'.  These
    %fields not necessary is 'strTyp' is numeric.

%Output:
%'outTyp' = same as 'strTyp'
%'yrs' = years to load
%'mnths' = months to load


if ~isfield(sMeta,'mnths') && ~isfield(sMeta,'currTime')
    outTyp = strTyp;
    yrs = nan;
    mnths = nan;
    return
elseif isfield(sMeta,'currTime') && ischar(strTyp) && ~regexpbl(strTyp,'ts')
    outTyp = sMeta.currTime;
    yrs = sMeta.currTime(1);
    mnths = sMeta.currTime(2);
    return
elseif isfield(sMeta,'currTime') && isnumeric(strTyp)
    if all(strTyp > 12) && numel(strTyp) == 2
        outTyp = strTyp;
        yrs = strTyp;
        mnths = sMeta.currTime(2);
        return
    end
end


if iscell(strTyp) && ~isempty(strTyp)
    outTyp = strTyp{1};
elseif ischar(strTyp)
    outTyp = strTyp;
elseif isnumeric(strTyp)
    if sum(size(strTyp) > 1) > 1
        yrs = strTyp(:,1);
        mnths = strTyp(:,2);
        
        outTyp = strTyp;
        %error('output_type:multEntries','The first input argument appears to contain multiple time-series entries.  This function can only handle one at a time.');
    elseif sum(strTyp > 12) == numel(strTyp)
        if numel(strTyp) == 2
            yrs = (strTyp(1):strTyp(2))';
            mnths = sMeta.currTime(2);
        else
            yrs = strTyp;
            mnths = sMeta.currTime(2);
        end
        
        outTyp = [yrs,mnths*ones(numel(yrs),1)];
    elseif numel(strTyp) == 2
        yrs = strTyp(1,1);
        mnths = strTyp(1,2);
        outTyp = strTyp;
    else
        error('output_type:unknownNum','This case has not been coded for.');
    end
    
        
    return
else
    outTyp = '';
end


if isempty(outTyp)
    outTyp = 'ts';
end

if ~isnumeric(outTyp) && isempty(outTyp)
    yrs = NaN(2,1);
    mnths = NaN;
elseif ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'ts'))
    yrs = sMeta.yrsOut;
    if isfield(sMeta,'currTime')
        mnths = sMeta.currTime(2);
    else
        mnths = sMeta.mnths;
    end
elseif ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'c'))
    yrs = sMeta.yrsClim;
    if isfield(sMeta,'currTime')
        mnths = sMeta.currTime(2);
    else
        mnths = sMeta.mnths;
    end
elseif ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'all'))
    yrs = sMeta.yrsOut;
    mnths = sMeta.mnths;
elseif ~isnumeric(outTyp) && ~isempty(regexpi(outTyp,'none'))
    yrs = NaN(2,1);
    mnths = [];
elseif isnumeric(outTyp)
    yrs = outTyp;
    if isfield(sMeta,'currTime')
        mnths = sMeta.currTime(2);
    else
        mnths = sMeta.mnths;
    end
else
    error('read_geodata:vararg',['The optional argument is set to ' ...
        outTyp ', which is not a valid option.']);
end

%Warning and switch if yr(1) > yr(2)
if length(yrs(:)) > 1 && yrs(1) > yrs(2) 
    warning('gcm_load:invalidEndYr', ['The requested ending year of ' ...
        'GCM data to load is less than the starting year.  Therefore, ' ...
        'the ending year will be changed from ' num2str(yrs(2)) ...
        ' to ' num2str(yrs(1)) ', which is the beginning year.']);
    yrs = flipud(yrs);
end


end

