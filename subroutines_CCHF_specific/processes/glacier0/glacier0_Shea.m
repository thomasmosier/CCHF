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

function glacier0_Shea(sMeta)
%output is ice thickness (in terms of water equivalent), assuming glacier-
%climate equilibrium.
global sCryo

%Equilibrium glacier thickness (Eq. 8) of
%Shea, J. M., Immerzeel, W. W., Wagnon, P., Vincent, C., & Bajracharya, S. 
%(2015). Modelling glacier change in the Everest region, Nepal Himalaya. 
%The Cryosphere, 9(3), 1105–1128. http://doi.org/10.5194/tc-9-1105-2015

% disp('initializing ice initial conditions (e.g. water quivalent)')
tauNaught = find_att(sMeta.global,'equil_shear'); 
rhoI = find_att(sMeta.global,'density_ice'); 
rhoW = find_att(sMeta.global,'density_water');
g = find_att(sMeta.global,'grav_accel');

varAngle = 'igrdslopeangfdr';
varRiseRun = 'igrdslopefdr';
%Calculate angle of slope from rise/run
if ~isfield(sCryo, varAngle)
    sCryo.(varAngle) = real(atand(sCryo.(varRiseRun)));
    %Set minimal angle to 1.5 degrees to ensure no stagnant ice (based on Shea et al.)
    sCryo.(varAngle)(sCryo.(varAngle) > -1.5 & sCryo.(varAngle) < 0) = -1.5; 
    sCryo.(varAngle)(sCryo.(varAngle) >= 0 & sCryo.(varAngle) <= 1.5) = 1.5;
end

gThick = abs(tauNaught./(rhoI*g*sind(sCryo.(varAngle))));
gThick(gThick < 0) = 0;
% gThick(isnan(sCryo.icwe)) = nan;
% gThick(sCryo.icwe == 0) = 0;

if issparse(sCryo.igrdwe)
    gThick = sparse(double(gThick));
else
    gThick = single(gThick);
end

sCryo.igrdwe = (rhoI/rhoW)*gThick;

if issparse(sCryo.icwe)
    sCryo.igrdbsdem = sCryo.igrddem - full(gThick);
else
    sCryo.igrdbsdem = sCryo.igrddem - gThick;
end
