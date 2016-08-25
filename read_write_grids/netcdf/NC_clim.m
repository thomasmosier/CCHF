function [sClim] = NC_clim(path,sMeta)
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


[sClim] = NC_load(path,sMeta,1);

nYrsUse = abs(strucMeta.yrsClim(2) - strucMeta.yrsClim(1))+1;

if nYrsUse == length(sClim.Data(:,1,1))
    gcmClimTemp = zeros(size(sClim.Data(1,:,:)));
    for ii = 1 : nYrsUse
        gcmClimTemp = gcmClimTemp + sClim.Data(ii,:,:);
    end

    sClim.Data = gcmClimTemp / nYrsUse;
else
    error('gcm_load:climMnth',['The function is in climatology '...
        'mode, although there is an error with the number of '...
        'time-series elements present.']);
end

sClim.Data = squeeze(sClim.Data);
