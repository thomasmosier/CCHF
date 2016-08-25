function gcmTs = NC_t_remove(gcmTs, mnthsUse, mnthsAvail, nYrsUse)
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


%This function presumes that 'gcmTs' is NetCDF data for continuous years.  
%It picks out which time-series elements to use based upon 

%Compare months to use with months available and index the former
%accordingly:
mnthsUseInd = zeros(length(mnthsAvail),1);
for ii = 1 : length(mnthsUse)
    indTemp = find(mnthsAvail == mnthsUse(ii), 1);
    
    if ~isempty(indTemp)
        mnthsUseInd(indTemp) = 1;
    end
end

if nYrsUse >= 1
    mnthIndKp = [];
    for ii = 1 : nYrsUse
        mnthIndKp = [mnthIndKp, mnthsUseInd];
    end
else
    error('gcmLoad:noYrs','No years of GCM data have been loaded.');
end

gcmTs = gcmTs(mnthIndKp == 1,:,:);