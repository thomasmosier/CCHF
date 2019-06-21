function [data, date] = read_DHM_data(path)
%read_DHM_data reads streamgage data from Nepal DHM
%
% Input
%   path - complete path and filename to DHM streamgage file
%
% Optional Input n/a
%   
% Output
%   data - numDays X 1 numeric array of daily flow data, units = m^3/s
%   date = numDays X 3 numberic array of year, month, day
%
% Notes
% This function parses streamgage data obtained from the Nepal Government
% Department of Hydrology and Meteorology and reformats in a format
% expected by CCHF.
%
% DHM data includes confidence intervals and possibly multiple 
% measurements per day.  CCHF expects one flow value per day.
%
% Reads multi-column DHM data, ignoring lowerBound and upperBound
% measurements, and averaging meanDischarge values for any days with
% more than one measurement.
%
% TODO: Currently missing dates are passed through, do they need to
% be filled with NaNs?
%
% Throws errors if:
%   Input file is not .csv
%   File does not contain expected column names from DHM
% 
% Example
% [data, date] = read_DHM_data('/path/to/gage_from_DHM.csv');
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

if regexpbl(extGage,{'csv'})  

    fprintf("%s: Reading csv file %s", mfilename(), path);
    T = readtable(path);
    
    % Confirm that data columns are what we are expecting from DHM
    % N.B. this logical depends on readtable behavior that collapses
    % whitespace in column headers
    if ~all(ismember(T.Properties.VariableNames, ...
            {'Date', 'Time', 'lowerBound', ...
            'meanDischarge', 'upperBound'}))
        error("%s: ERROR: DHM file %s does not contain expected columns", ...
            mfilename(), path);
    end
    
    % Drop the confidence interval columns we will ignore
    % This feature is only available beginning with Matlab r2018a
    T.lowerBound = [];
    T.upperBound = [];
    
    % Find any dates with multiple values and
    % average them
    [G, TID] = findgroups(T(:,'Date'));
    dailyDischarge = splitapply(@mean, T(:, 'meanDischarge'), G);
    dailyData = table(TID, dailyDischarge);
    
    % Set any dates with missing discharge data to NaN
    startDate = dailyData{1, "TID"}.Date;
    endDate = dailyData{height(dailyData), "TID"}.Date;
    date = (startDate:endDate)';
    discharge = NaN([numel(date) 1]);
    allData = [table(date), table(discharge)];

    for i=1:height(dailyData)
        nextDate = dailyData{i, "TID"}.Date;
        nextDischarge = dailyData{i, "dailyDischarge"};
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
    error("%s: ERROR: DHM file %s in unknownFormat: %s\n", ...
        mfilename(), path, extGage);
end

