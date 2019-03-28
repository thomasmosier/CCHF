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

function dhfrl = dflux_long_simple(sMeta)
%R_long from atm is parameterized using Deacon (1970)

global sCryo sAtm

%VERSION WITH ONLY ONE FITTING PARAMETER



eSkyClear = find_att(sMeta.global, 'emm_clear');
eSkyCloud = find_att(sMeta.global, 'emm_cloud');

if ~isfield(sAtm,'emm')
   sAtm.emm = eSkyClear*ones(size(sAtm.rain)); 
end    

sAtm.emm(squeeze(sAtm.pr(sAtm.indCurr,:,:)) == 0) = eSkyClear;
sAtm.emm(squeeze(sAtm.pr(sAtm.indCurr,:,:)) > 0) = eSkyCloud;


% stefan = 5.67*10^(-8);

%Reference:
%Plüss, C., & Ohmura, A. (1997). Longwave radiation on snow-covered 
%mountainous surfaces. Journal of Applied Meteorology, 36(6), 818–824.

if isfield(sCryo,'tsis')
    dhfrl = -5.67*10^(-8)*(sCryo.tsis + 273.15).^3;
elseif ~isfield(sCryo,'tsis') && isfield(sCryo,'tsn')
    dhfrl =  -5.67*10^(-8)*(sCryo.tsn + 273.15).^3;
else
    error('heat_long_simple:missingTmp',['Neither internal or surface '...
        'snow temperature are known.']);
end

dhfrl(isnan(dhfrl)) = 0;