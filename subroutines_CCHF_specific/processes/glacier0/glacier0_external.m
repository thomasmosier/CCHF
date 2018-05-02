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


rhoI = find_att(sMeta.global,'density_ice'); 
rhoW = find_att(sMeta.global,'density_water');


if ~isfield(sCryo, 'igrddem')
    error('glacier0External:noArea', ['Glacier base elevation cannot be '...
        'estimated because no ice grid DEM is available (field name = igrddem)']);
end

if isfield(sCryo, 'icwe')
    sCryo.igrdbsdem = sCryo.igrddem - full(sCryo.icwe)*(rhoW/rhoI);
elseif isfield(sCryo, 'igrdwe')
    sCryo.igrdbsdem = sCryo.igrddem - full(sCryo.igrdwe)*(rhoW/rhoI);
elseif isfield(sCryo, 'icdepth')
    sCryo.igrdwe = (rhoI/rhoW)*full(sCryo.icdepth);
    
    sCryo.igrdbsdem = sCryo.igrddem - full(sCryo.icdepth);
elseif isfield(sCryo, 'igrddepth')
    sCryo.igrdwe = (rhoI/rhoW)*full(sCryo.igrddepth);
    
    sCryo.igrdbsdem = sHydro.dem - full(sCryo.igrddepth);
else
    error('glacier0External:noArea','Ice thickness and base DEM cannot be initialized because there is no depth or thickness approximation.');
end
