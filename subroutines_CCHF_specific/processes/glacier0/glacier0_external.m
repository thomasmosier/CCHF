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

function glacier0_external(sHydro, sMeta)
%output is ice thickness (in terms of water equivalent), assuming glacier-
%climate equilibrium.
global sCryo

%Empirical ice thickness based on
%Chen, J., & Ohmura, A. (1990). Estimation of Alpine glacier water 
%resources and their change since the 1870s. IAHS publ, 193, 127-135.

rhoI = find_att(sMeta.global,'density_ice'); 
rhoW = find_att(sMeta.global,'density_water');

if isfield(sCryo, 'icwe')
    sCryo.icbsdem = sHydro.dem - full(sCryo.icwe)*(rhoW/rhoI);
elseif isfield(sCryo, 'icdepth')
    sCryo.icwe = (rhoI/rhoW)*full(sCryo.icdepth);
    
    sCryo.icbsdem = sHydro.dem - full(sCryo.icdepth);
else
    error('glacier0External:noArea','Ice thickness cannot be estimated because there is no area approximation.');
end
