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

function heat_precip_simple(varargin)

global sCryo sAtm


%VERSION WITH ONLY ONE FITTING PARAMETER
% if isempty(varargin(:))
% 	varargout{1} = cell(0,6);
%     varargout{1} = cat(1, varargout{1}, {'pre_heat_pwr', -2, 2, 0, 'heat_precip_simple', 'cryo'});
%     return
% else
%     pPwr = find_att(varargin{1}.coef,'pre_heat_pwr'); 
% end

%Find seconds in current time-step:
dt = time2sec(1,varargin{1}.dt,varargin{1}.currTime);

% 
% %Create field to store previous atmospheric temperature grid:
% if ~isfield(sCryo, 'tasp')
%     sCryo.tasp = squeeze(sAtm.tas(sAtm.indCurr,:,:));
% end

cSens   = 4.18*10^6; %= c_water*density_water: J / m^3 / Deg C
cLate   = 3.34*10^8; %= L_water*density_water: J / m^3

%Sensible heat (units of Watts) = density of water (1000 kg/m^3) * specific heat capacity water (4179 J/kg/Deg C) * m of precipitation * (airTemp - snowTemp)
sCryo.heatPC = (cSens/dt)*squeeze(sAtm.pr(sAtm.indCurr,:,:)).*squeeze(sAtm.tas(sAtm.indCurr,:,:));

%latent heat (units of Watts) = density of water (1000 kg/m^3) * latent heat of fusion of water (334000 J/kg) * m of rain
%Add latent heat of rain impinging on snowpack
indFreeze = find(sCryo.solid > 0 | sCryo.ice > 0);
sCryo.heatPC(indFreeze) = sCryo.heatPC(indFreeze) + (cLate/dt)*sAtm.rain(indFreeze);
    sCryo.heatPC(isnan(sCryo.heatPC)) = 0;
 
% sCryo.heatPC = 10^(pPwr)*sCryo.heatPC;
    