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

function day = Julian_day(currTime, varargin)

%currTime = [year, month, day, ...]
%day = Julian day

%Test relative to http://quasar.as.utexas.edu/BillInfo/JulianDateCalc.html

% %Start of Julian Calendar in Julian is Jan 1, 4713 BC:
% day = days_since([-4713, 1, 1],currTime(:,1:3),'gregorian'); 



%If inputs have hour and/or minutes colum, convert to day decimal
if numel(currTime(1,:)) == 4
    currTime = [currTime(:,1:2), currTime(:,3) + currTime(:,4)/24];
elseif numel(currTime(1,:)) == 5
    currTime = [currTime(:,1:2), currTime(:,3) + currTime(:,4)/24 + currTime(:,5)/(60*24)];
end

day = nan(numel(currTime(:,1)),1);

for ii = 1 : numel(day)
    if ~isempty(varargin(:)) && regexpbl(varargin{1},'ofyear')
       day(ii) = ...
           datenum([num2str(currTime(ii,2)) '/' num2str(currTime(ii,3)) '/' num2str(currTime(ii,1))], 'mm/dd/yyyy') ...
           - datenum(['01/01/' num2str(currTime(ii,1))], 'mm/dd/yyyy') ...
           + 1;
    else
        %Difference between date number and Julian date is 1721058.5
        if numel(currTime(1,:)) == 3
            day(ii) = 1721058.5 + ...
                datenum([num2str(currTime(ii,2)) '/' num2str(currTime(ii,3)) '/' num2str(currTime(ii,1))], 'mm/dd/yyyy');
        else
           error('Julian_date:inputLength',['The indices in the current input is ' num2str(numel(currTime(ii,:))) ', which has not been programmed for.']); 
        end
    end
end

