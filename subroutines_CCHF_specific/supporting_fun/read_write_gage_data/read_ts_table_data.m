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


function [data, date] = read_ts_table_data(path, inputUnits)





sepDate = {'\','/','-',':'};

[~, ~, extGage] = fileparts(path);

% rangeRead = [];
% if ~isempty(varargin(:))
%     rangeRead = varargin{1};
% end

if regexpbl(extGage,{'xls'})
    [dataTemp, txtRaw] = xlsread(path);

    %Find first row of dates in text output:
    rStart = [];
    for ii = 1 : numel(txtRaw(:,1))
       if strcmpi(txtRaw(ii,1),'usgs')
           rStart = ii;
           break
       end
    end
    
    if ~isempty(rStart) %USGS case
        [data, date] = read_USGS_data(path);
    else %Genereric case
        if numel(dataTemp(1,:)) ~= numel(txtRaw(1,:))
           error('readtstabledata:xlsformat',['The Excel data in ' path ...
               ' have unexpected dimensions. This is likely caused '...
               'by columns not formatted as numbers.']) 
        end
        
        rStart = 1;
        
        blYr = nan;
        blMnth = nan;
        blDay = nan;
        blData = nan;
        
        for jj = 1 : numel(txtRaw(rStart,:))
            if regexpbl(txtRaw(rStart,jj), {'year', 'yr'})
                blYr = jj;
            elseif regexpbl(txtRaw(rStart,jj), {'month'})
                blMnth = jj;
            elseif regexpbl(txtRaw(rStart,jj), {'day'})
                blDay = jj;
            elseif regexpbl(txtRaw(rStart,jj), {'discharge', 'flow', 'stream', 'river', 'volume'})
                blData = jj;
            end
        end
        
        data = dataTemp(:,blData);
        
        if ~isnan(blYr) && ~isnan(blMnth) && ~isnan(blDay)
            date = [dataTemp(:,blYr), dataTemp(:,blMnth), dataTemp(:,blDay)];
        elseif ~isnan(blYr) && ~isnan(blMnth)
            date = [dataTemp(:,blYr), dataTemp(:,blMnth)];
        else
            error('readtstabledata:missingDate1',['One of more date '...
                'indices in the Excel file (' path ')could not be found.']);
        end
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
    
    if ~isempty(rStart) %USGS case:
       [data, date] = read_USGS_data(path);
    else %generic case:
        error('readtstabledata:unreadableText',['.txt gage data in a '...
            'format other than USGS convention has not been programmed for.']);
    end


elseif regexpbl(extGage,{'csv'})
    dataTemp = csvread(path,1,0);
   
    blYr = nan;
    blMnth = nan;
    blDay = nan;

    for jj = 1 : numel(dataTemp(1,:))
        if all(dataTemp(:,jj) > 31) && max(diff(dataTemp(:,jj))) < 3
            blYr = jj;
        elseif all(dataTemp(:,jj) < 13) && max(diff(dataTemp(:,jj))) < 14
            blMnth = jj;
        elseif all(dataTemp(:,jj) < 35) && ~all(dataTemp(:,jj) < 13) && max(diff(dataTemp(:,jj))) < 35
            blDay = jj;
        end
    end
    
    if numel(dataTemp(1,:)) == 4
        if any([isnan(blYr), isnan(blMnth), isnan(blDay)])
            error('readtstabledata:missingDate2',['One of more date '...
                'indices in the Excel file (' path ')could not be found.']);
        end
        
        blData = setdiff((1:4), unique([blYr, blMnth, blDay]));
        
        date = [dataTemp(:,blYr), dataTemp(:,blMnth), dataTemp(:,blDay)];
    elseif numel(dataTemp(1,:)) == 3
        if any([isnan(blYr), isnan(blMnth)])
            error('readtstabledata:missingDate3',['One of more date '...
                'indices in the Excel file (' path ')could not be found.']);
        end
        
        blData = setdiff((1:4), unique([blYr, blMnth]));
        
        date = [dataTemp(:,blYr), dataTemp(:,blMnth)];
    else
        error('readtstabledata:csvFormat',['The CSV data in ' path ...
            ' have unexpected dimensions. This is likely caused '...
            'by columns not formatted as numbers.']) 
    end
    
    data = dataTemp(:,blData);
else
    error('read_ts_table_data:unknownFormat',['Time-series data in ' char(39) ...
        extGage char(39) ' format has not been programmed for.']);
end

%Convert units:
if regexpbl(inputUnits, 'ft3/s')
    %Convert from ft^3/s to m^3/s because USGS data typically in the former and
    %the SETI model uses standard metric units
    warning('read_ts_table_data:convert',['The input gage data are being '...
        'converted from ft^3/s to m^3/s.  This is supposed to happen as '...
        'long as the input are in units of ft^3/s.']);
    data = 0.0283168*data;
elseif ~regexpbl(inputUnits, 'm3/s')
    error('read_ts_tabled_ata:unknownUnits', ['The input units ' ...
        inputUnits  ' are not known.']);
end



% %%Convert dates in text format to numeric array format:
% %Identify seperator used:
% for ll = 1 : numel(sepDate(:))
%     if regexpbl(txtRaw{1},sepDate{ll});
%         sepUse = sepDate{ll}; 
%     end
% end

% colDateOut = numel(regexpi(txtRaw{1},sepUse)) + 1;
% %Initialize output array:
% date = nan(numel(txtRaw(:)), colDateOut);
% %Convert each date:
% for kk = 1 : numel(txtRaw(:))
%     indSep = regexpi(txtRaw{kk},sepUse);
%     
%     if numel(date(1,:)) == 2
%         date(kk,:) = [ ...
%             str2double(txtRaw{kk}(1:indSep(1)-1)), ...
%             str2double(txtRaw{kk}(indSep(1)+1:end))];
%     elseif numel(date(1,:)) == 3
%         date(kk,:) = [ ...
%             str2double(txtRaw{kk}(1:indSep(1)-1)), ...
%             str2double(txtRaw{kk}(indSep(1)+1:indSep(2)-1)), ...
%             str2double(txtRaw{kk}(indSep(2)+1:end))];
%     else
%         error('read_USGS_data:DateNumber',['There are ' ...
%             num2str(numel(date(1,:))) ' date values.  This has not '...
%             'been programmed for.']);
%     end
% end


if numel(data(:,1)) ~= numel(date(:,1))
        error('read_USGS_data:diffLength',['There are ' ...
            num2str(numel(data(:,1))) ' entries in the data array and '...
            num2str(numel(date(:,1))) ' entries in the date array.  '...
            'The cause of this mismatch needs to be investigated.']);
end

