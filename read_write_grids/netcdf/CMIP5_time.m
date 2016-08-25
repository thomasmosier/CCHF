function varargout = CMIP5_time(fileData)
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


%Find yrs and mnths of CMIP5 GCM data and CRU data based upon name string.


%Initialize output:
yrsGcm = nan(2,1);
mnthsGcm = nan(2,1);
daysGcm = nan(2,1);


if exist(fileData,'file')
    [~, fileData, ~] = fileparts(fileData);
end
indDate = regexpi(fileData,'_');
indPer = regexpi(fileData,'\.');

if ~isempty(indDate) && isempty(indPer) %CMIP5 NetCDF format
    fileData = fileData(indDate(end)+1:end);
    indDateSep = regexpi(fileData,'-');
    yrsGcm = [str2double(fileData(1:4)); str2double(fileData(indDateSep +1 : indDateSep +4))];
    mnthsGcm = [str2double(fileData(5:6)); str2double(fileData(indDateSep +5 : indDateSep +6))];
    if ~isnan(str2double(fileData(7)))
        daysGcm = [str2double(fileData(7:8)); str2double(fileData(indDateSep +7 : indDateSep +8))];
    end
elseif ~isempty(indPer) %Maybe CRU NetCDF data or CFSR?
    indNum = regexpi(fileData,'(\d+)');
    dates = regexpi(fileData,'(\d+)','tokens');
    if ~isempty(dates)
        cntDates = 1;
        mnthsTemp = nan(2,1);
        for ii = 1 : length(dates(:))
            currNum = char(dates{ii});
            if numel(currNum) == 4
                yrsGcm(cntDates) = str2double(currNum);
                cntDates = cntDates + 1;
            elseif numel(currNum) == 6 %|| numel(currNum) == 7
                yrsGcm(cntDates) = str2double(currNum(1:4));
                mnthsGcm(cntDates) = str2double(currNum(5:end));
                cntDates = cntDates + 1;
            elseif numel(currNum) == 8
                yrsGcm(   cntDates) = str2double(currNum(1:4));
                mnthsGcm(cntDates) = str2double(currNum(5:6));
                daysGcm(  cntDates) = str2double(currNum(7:8));
                cntDates = cntDates + 1;
            end
        end
        if cntDates > 3
            yrsGcm = nan(2,1);
        end
        yrsGcm = sort(yrsGcm);
        if sum(~isnan(yrsGcm)) == 2 && yrsGcm(1) ~= yrsGcm(2) && isequal(mnthsTemp,[1;12])
            mnthsGcm = [1;12];
        end
    end
end

%If ending date nan, set to be same as beginning date:
if ~isnan(yrsGcm(1)) && isnan(yrsGcm(2))
    yrsGcm(2) = yrsGcm(1);
end
if ~isnan(mnthsGcm(1)) && isnan(mnthsGcm(2))
    mnthsGcm(2) = mnthsGcm(1);
end
if ~isnan(daysGcm(1)) && isnan(daysGcm(2))
    daysGcm(2) = daysGcm(1);
end

%write special code for GCM with units of 'BP' (before present)
if all(isnan(yrsGcm)) && regexpbl(fileData,'BP')
    yrsGcm = [-9999;-9999];
    mnthsGcm = [-9999;-9999];
end
if numel(yrsGcm) ~= 2
    warning('CMIP5_time:nYrs',['This function identified ' ...
        num2str(numel(yrsGcm)) ' that may be year references.  ' ...
        'Using normal CMIP5 naming convention, there should be a ' ...
        'starting year and ending year.']);
end

if nargout <= 1
    varargout{1} = [yrsGcm, mnthsGcm, daysGcm];
elseif nargout == 2
   varargout{1} = yrsGcm;
   varargout{2} = mnthsGcm;
elseif nargout == 3
   varargout{1} = yrsGcm;
   varargout{2} = mnthsGcm;
   varargout{3} = daysGcm; 
end
