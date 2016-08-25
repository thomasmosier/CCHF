function days = days_in_month(dateVec,strCal)
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


%Only want year and month:
if length(dateVec(1,:)) == 3
    dateVec = dateVec(:,1:2);
end

blCal = cal_class(strCal);
 
days = nan(length(dateVec(:,1)),1);

if blCal == 0 || blCal == 1 
    for ii = 1 : length(dateVec(:,1))
       days(ii) = eomday(dateVec(ii,1), dateVec(ii,2));
       %Subtract leap if calednar doesn't include it:
       if blCal == 1 && dateVec(ii,2) == 2
           days(ii) = days(ii) - 1;
       end
    end
elseif blCal == 2
    days(:) = 30;
else
    error('days_in_month:unknownCal',[char(39) num2str(blCal) char(39) ' is an unknown calander class from function ' char(39) 'cal_class' chaR(39)'.']);
end