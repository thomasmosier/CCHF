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

function varargout = flux_short(strType, varargin)

%varargin = sMeta required to evaluate incoming shortwave radiation
global sCryo sAtm sLand

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    return
end



%Find top of the atmosphere radiation:
indRSDT = find(ismember(sLand.rsdtdate(:,2:end),varargin{1}.currTime(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end

if regexpbl(strType,{'snow','swe'})
    sCryo.heatSR = (1-sCryo.albedoS).*sAtm.rstran.*squeeze(sLand.rsdt(indRSDT,:,:));
elseif regexpbl(strType,{'ice','glacier'})
    sCryo.heatSRI = (1-sCryo.albedoI).*sAtm.rstran.*squeeze(sLand.rsdt(indRSDT,:,:));
end
    
    
