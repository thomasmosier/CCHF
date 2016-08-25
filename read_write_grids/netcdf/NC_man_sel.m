function ind = NC_man_sel(nType, varAvail, varFind)
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


%Display GUI for user to find variable if the automated system was unable to find it. 
if sum(nType) == 0    
    disp('');   %New line before prompt.
    for ii = 1 : length(varAvail)
        disp([num2str(ii) ') ' char(varAvail{ii})]);
    end
    indStr = input(['Enter the number of the GCM variable ' ...
        'corresponding to ' varFind '.' char(10) 'Enter a blank string if the '...
        'climate parameter is not present.' char(10)],'s');
    ind = str2double(indStr);
    if isnan(ind)
        ind = [];
    end
elseif sum(nType) == 1
    ind = find(nType ~= 0);
else
    error('gcm_load:varFind', ['Multiple variable indices selected.  '...
        'No clause has been written for this case.'])
end
if sum(nType) ~= 1 && isempty(ind)
    error('gcm_load:varFind', ['The script was unable to detect which '...
       'variable of the GCM time-series file represents ' ...
       varFind '.']); 
end