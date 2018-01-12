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

function dhfsnc = dflux_conduct_snow(sMeta)

%Positive values indicate heat flowing from internal snowpack to snow-air interface
%Negative values indicate heat flowing from snow-air interface to internal snowpack


global sCryo

c = find_att(sMeta.coef,'cond_pwr');

%Set maximum value for conduction:
condtMax = 50; %Units = Watts per meter squared

rhoDry = find_att(sMeta.global, 'density_snow_dry');
rhoWet = find_att(sMeta.global, 'density_snow_wet');
rhoWater = find_att(sMeta.global, 'density_water');
kDry = find_att(sMeta.global, 'thermal_conduct_snow_dry');
kWet = find_att(sMeta.global, 'thermal_conduct_snow_wet'); 

if isfield(sCryo, 'tsn') && isfield(sCryo, 'tsis') 
    dhfsnc = zeros(size(sCryo.snw),'single');
    
    indWet = find(sCryo.snlq > 0.01 * sCryo.snw); %this value somewhat arbitrary
    indDry = setdiff((1:numel(sCryo.snw)), indWet);
    dhfsnc(indWet) = -10^(c)*(2*kWet*rhoWet/rhoWater) ./ sCryo.snw(indWet);
    dhfsnc(indDry) = -10^(c)*(2*kDry*rhoDry/rhoWater) ./ sCryo.snw(indDry);
else
    error('heat_conduct_snow:missingTmp','Surface or internal snow temperature not defined.');
end

%Ensure no odd values (due to there not being any snow:
dhfsnc(isnan(dhfsnc)) = 0;
dhfsnc(sCryo.snw == 0) = 0;
dhfsnc(dhfsnc > condtMax) = condtMax;
dhfsnc(dhfsnc < -condtMax) = -condtMax;

%If temperatures equal, then no conduction
dhfsnc(sCryo.tsn == sCryo.tsis) = 0;