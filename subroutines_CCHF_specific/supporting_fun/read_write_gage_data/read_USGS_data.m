% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


function [data, date] = read_USGS_data(path)





sepDate = {'\','/','-',':'};

[~, ~, extGage] = fileparts(path);

% rangeRead = [];
% if ~isempty(varargin(:))
%     rangeRead = varargin{1};
% end

if regexpbl(extGage,{'xls'})
    [data, txtRaw] = xlsread(path);

    data = data(:,end);
    
    %Find first row of dates in text output:
    rStart = [];
    for ii = 1 : numel(txtRaw(:,1))
       if strcmpi(txtRaw(ii,1),'usgs')
           rStart = ii;
           break
       end
    end

    if isempty(rStart)
       error('read_USGS_data:unknownTextFormat',['The streamgage data is '...
           'not is USGS format.  The current case has not been programmed for.']);
    end

    %Find column of dates in txt output:
    colDate = [];
    for jj = 1 : numel(txtRaw(rStart,:))
        if regexpbl(txtRaw(rStart,jj),sepDate)
            colDate = jj;
            break
        end
    end

    if isempty(colDate)
       error('read_USGS_data:unknownDateFormat',['The column containing '...
           'date information could not be found.']);
    else
        txtRaw = txtRaw(rStart:end, colDate);
    end
elseif regexpbl(extGage,{'txt'})  
    fid = fopen(path, 'rt');
    cellRaw = textscan(fid,'%s\t%s\t%s\t%s\t%s');
    fclose(fid);
    
    rStart = []; cntR = 1;
    while isempty(rStart)
        if regexpbl(cellRaw{1,1}(cntR),'usgs')
           rStart = cntR; 
        end
        cntR = cntR +1;
    end

    data = str2double(cellRaw{1,4}(rStart:end));
    txtRaw = cellRaw{1,3}(rStart:end);
else
    error('read_USGS_data:unknownFormat',['USGS data in ' char(39) ...
        extGage char(39) ' format has not been programmed for.']);
end

% %Convert from ft^3/s to m^3/s because USGS data typically in the former and
% %the SETI model uses standard metric units
% warning('read_USGS_data:convert',['The input gage data are being '...
%     'converted from ft^3/s to m^3/s.  This is supposed to happen as '...
%     'long as the input are in units of ft^3/s.']);
% data = 0.0283168*data;



%%Convert dates in text format to numeric array format:
%Identify seperator used:
for ll = 1 : numel(sepDate(:))
    if regexpbl(txtRaw{1},sepDate{ll});
        sepUse = sepDate{ll}; 
    end
end

colDateOut = numel(regexpi(txtRaw{1},sepUse)) + 1;
%Initialize output array:
date = nan(numel(txtRaw(:)), colDateOut);
%Convert each date:
for kk = 1 : numel(txtRaw(:))
    indSep = regexpi(txtRaw{kk},sepUse);
    
    if numel(date(1,:)) == 2
        date(kk,:) = [ ...
            str2double(txtRaw{kk}(1:indSep(1)-1)), ...
            str2double(txtRaw{kk}(indSep(1)+1:end))];
    elseif numel(date(1,:)) == 3
        date(kk,:) = [ ...
            str2double(txtRaw{kk}(1:indSep(1)-1)), ...
            str2double(txtRaw{kk}(indSep(1)+1:indSep(2)-1)), ...
            str2double(txtRaw{kk}(indSep(2)+1:end))];
    else
        error('read_USGS_data:DateNumber',['There are ' ...
            num2str(numel(date(1,:))) ' date values.  This has not '...
            'been programmed for.']);
    end
end

%Swap date columns:
if numel(date(1,:)) == 2
    date = [date(:,2), date(:,1)];
elseif numel(date(1,:)) == 3
    date = [date(:,3), date(:,1), date(:,2)];
else
    error('read_USGS_data:DateNumber',['There are ' ...
        num2str(numel(date(1,:))) ' date values.  This has not '...
        'been programmed for.']);
end



if numel(data(:,1)) ~= numel(date(:,1))
        error('read_USGS_data:diffLength',['There are ' ...
            num2str(numel(data(:,1))) ' entries in the data array and '...
            num2str(numel(date(:,1))) ' entries in the date array.  '...
            'The cause of this mismatch needs to be investigated.']);
end

