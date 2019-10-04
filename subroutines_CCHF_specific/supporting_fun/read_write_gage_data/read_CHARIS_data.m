function [data, date, latitude, longitude] = ...
    read_CHARIS_data(path)
%read_CHARIS_data reads streamgage data format standardized by CHARIS
%
% Input
%   path - complete path and filename to CHARIS streamgage file
%
% Optional Input n/a
%   
% Output
%   data - numDays X 1 numeric array of daily flow data
%   date - numDays X 3 numberic array of year, month, day
%   latitude, longitude - lat/lon of gage, as parsed from file, or NaN
%
% Notes
% This function parses streamgage data obtained from various
% locations standardized/QC'd by the NSIDC CHARIS project to
% contain:
% 1) Any number of header lines beginning with '# '
% 2) Header latitude line: '# Latitude: dd.dddd'
% 3) Header longitude line: '# Longitude: dd.dddd'
% 4) Header units line: '# Units: m**3 s**-1'
% 5) Header COLUMN Names line: '# COLUMNS: year month day doy
% discharge'
% 6) Time series of discharge records as
% yyyy mm dd doy value
%
% CHARIS data is expected to contain on flow value per day.
% Some dates may be missing, they will be filled with NaNs on
% output

% Throws errors if:
%   Input file is not .dat
%   File does not contain expected column names from CHARIS
%   File metadata units are not m3/s
% 
% Example
% [data, date] =
%   read_CHARIS_data('/path/to/CHARIS_gage_day_discharge.dat');
%

% Copyright 2019 Thomas M. Mosier and the Regents of the University
% of Colorado
%
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
%

sepDate = {'\','/','-',':'};

[~, ~, extGage] = fileparts(path);

if regexpbl(extGage,{'dat'})  

    fprintf("%s: Reading dat file %s", mfilename(), path);
    % Read the comments to get lat/lon and column names
    % Can potentially also confirm units...
    latitude = NaN;
    longitude = NaN;
    expectedColumns = '# COLUMNS: year month day doy discharge';
    expectedUnits = '# Units: m**3 s**-1';
    foundExpectedColumns = false;
    inputUnits = 'Unknown';
    
    fid = fopen(path);
    tline = fgetl(fid);
    while strcmp(tline(1:2), '# ')
        if strcmp(tline, expectedColumns)
            foundExpectedColumns = true;
        elseif strcmp(tline, expectedUnits)
            inputUnits = 'm3/s';
        else
            tokenNames = regexp(tline, ...
                '^# (?<parameter>.*):(?<value>.*)', ...
                'names');
            if ~isempty(tokenNames)
                if strcmp(tokenNames.parameter, 'Latitude')
                    latitude = str2num(tokenNames.value);
                elseif strcmp(tokenNames.parameter, 'Longitude')
                    longitude = str2num(tokenNames.value);
                end
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    
    if ~foundExpectedColumns
        error('read_CHARIS_data:unknownFormat', ...
            ['Unknown format of time-series data in ' path]);
    end
    
    if ~strcmp(inputUnits, 'm3/s')
        error('read_CHARIS_data:unknownUnits', ...
            ['Unrecognized units in ' path]);
    end
    
    opts = detectImportOptions(path);
    T = readtable(path, opts);

    % Set the expected column names
    T.Properties.VariableNames = {'year' 'month' 'day' 'doy' 'discharge'};
    
    % Set discharge values of -9999.99 to NaNs
    T.discharge(T.discharge == -9999.99) = NaN;
    
    inDates = datetime(T.year, T.month, T.day);
    startDate = inDates(1);
    endDate = inDates(end);
    
    T = [T table(inDates)];
    
    % Initialize date and data with all NaNs
    % so that missing dates are populated with NaN
    date = (startDate:endDate)';
    discharge = NaN([numel(date) 1]);
    allData = [table(date), table(discharge)];
    
    
    for i=1:height(T)
        nextDate = T{i, "inDates"};
        nextDischarge = T{i, "discharge"};
        allData(allData.date == nextDate, :).discharge = nextDischarge;
    end

    yyyy = year(allData.date);
    mm = month(allData.date);
    dd = day(allData.date);
    
    % Extract date parts    
    allData = [allData, table(yyyy), table(mm), table(dd)];
    
    % And reformat to expected output numeric matrices
    data = table2array(allData(:,"discharge"));
    date = table2array(allData(:,{'yyyy' 'mm' 'dd'}));
    
else
    error("%s: ERROR: CHARIS gage file %s in unknownFormat: %s\n", ...
        mfilename(), path, extGage);
end

