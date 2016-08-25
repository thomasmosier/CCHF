function strMnth = mnth_str(numMnth)
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


%Convert number between 1 and 12 to corresponding month name

if numMnth > 12 || numMnth < 1 || isnan(numMnth)
    strMnth = 'Unknown'; 
elseif numMnth == 1
    strMnth = 'January';
elseif numMnth == 2
    strMnth = 'February';
elseif numMnth == 3
    strMnth = 'March';
elseif numMnth == 4
    strMnth = 'April';
elseif numMnth == 5
    strMnth = 'May';
elseif numMnth == 6
    strMnth = 'June';
elseif numMnth == 7
    strMnth = 'July';
elseif numMnth == 8
    strMnth = 'August';
elseif numMnth == 9
    strMnth = 'September';
elseif numMnth == 10
    strMnth = 'October';
elseif numMnth == 11
    strMnth = 'November';
elseif numMnth == 12
    strMnth = 'December';
end


end

