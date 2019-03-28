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

function dateAdj = adj_dates(strDataRes, dateCurr)

%Find time resolution of model run from 'dateCurr'
numModRes = numel(dateCurr(1,:));

if regexpbl(strDataRes,{'year','annual'})
    numDataRes = 1;
elseif regexpbl(strDataRes, 'month')
    numDataRes = 2;
elseif regexpbl(strDataRes,{'day','daily'})
    numDataRes = 3;
elseif regexpbl(strDataRes,{'hourly'})
    numDataRes = 4;
else
    error('adj_dates:dataRes',['The current data resolution is ' ...
        strDataRes '.  This data resolution has not been programmed for.']);
end


if numDataRes <= numModRes 
    if numDataRes < 3
        dateAdj = [dateCurr(1:numDataRes) + [zeros(1,numDataRes - 1), -1]; dateCurr(1:numDataRes); dateCurr(1:numDataRes) + [zeros(1,numDataRes - 1), 1]];
        
        if dateAdj(1,2) == 0
            dateAdj(1,:) = [dateAdj(1,1) - 1, 12];
        end
        if dateAdj(3,2) == 13
            dateAdj(3,:) = [dateAdj(3,1) + 1, 1];
        end
        
        
    elseif numDataRes == 3
        %MUST EDIT DAYS_2_DATE to make this work
%         dateAdj = days_2_date([-1;0;1], [1980,1,1], 'gregorian'); %For
%         testing
        dateAdj = days_2_date([-1;0;1], dateCurr, 'gregorian');
    else
        error('adj_dates:dataRes',['The data time resolution is '...
            'greater than monthly, which has not been programmed for.']);
    end
elseif numDataRes > numModRes 
    error('adj_dates:dataRes',['The time resolution of the data is '...
        'greater than the time-resolution of the model.  This has not '...
        'been programmed for.']);
end