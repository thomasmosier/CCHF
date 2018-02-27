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

function varargout = flux_precip(varargin)

global sCryo sAtm

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    return
else
    sMeta = varargin{1};
end
    
%Find seconds in current time-step:
dt = time2sec(1,sMeta.dt,sMeta.dateCurr);

lWat = find_att(sMeta.global, 'heat_cap_water')*find_att(sMeta.global, 'density_water');

%Sensible heat (units of Watts) = density of water (1000 kg/m^3) * specific heat capacity water (4179 J/kg/Deg C) * m of precipitation * (airTemp - snowTemp)
if isfield(sCryo,'tss')
    sCryo.hfcp = (lWat/dt)*squeeze(sAtm.pr(sAtm.indCurr,:,:)) ...
        .*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) - sCryo.tsn);
elseif isfield(sCryo,'tsn')
    sCryo.hfcp = (lWat/dt)*squeeze(sAtm.pr(sAtm.indCurr,:,:)) ...
        .*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) - sCryo.tsn);
else
    error('heat_precip:missingTmp',['Neither internal or surface '...
        'snow temperature are known.']);
end

%%LATENT HEAT IS DEALTH WITH IN ENERGY FUNCTION!!!
% %latent heat (units of Watts) = density of water (1000 kg/m^3) * latent heat of fusion of water (334000 J/kg) * m of rain
% %Add latent heat of rain impinging on snowpack
% indFreeze = find(sCryo.tsn < 0);
% sCryo.hfcp(indFreeze) = sCryo.hfcp(indFreeze) + (334000000/dt)*sAtm.rain(indFreeze);
%     sCryo.hfcp(isnan(sCryo.hfcp)) = 0;