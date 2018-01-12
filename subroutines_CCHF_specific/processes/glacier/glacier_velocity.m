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

function glacier_velocity(sMeta)

%Weertman's sliding law (Eq. 7) of
%Shea, J. M., Immerzeel, W. W., Wagnon, P., Vincent, C., & Bajracharya, S. 
%(2015). Modelling glacier change in the Everest region, Nepal Himalaya. 
%The Cryosphere, 9(3), 1105–1128. http://doi.org/10.5194/tc-9-1105-2015

global sCryo

%Load constants:
tauNaught = find_att(sMeta.global,'equil_shear'); 
rho = find_att(sMeta.global,'density_ice'); 
g = find_att(sMeta.global,'grav_accel');

%Load assumptions:
nu = 0.1; %Bedrock roughness (unitless)
R = 1.80*10^9; %Material roughness coefficient (Pa m^{-2} s) - Shea et al. treat R as a fitting parameter (+-5.00X10^8)
n = 3; %creep cosntant of Glen's Flow law (unitless)

%Calculate the shear stress on the glaicer:
tau = sparse(max(0, rho*g*sCryo.iceWE.*abs(sCryo.iceSlope) - tauNaught)); 

sCryo.iceVel = (tau./(R*nu^2)).^((n+1)/2); %Units are m/s

%Note: Ice initialization uses similar formulation and assumes glaciers in
%equilibrium. Therefore, towards beginning of model runs, this formulation
%will likely return velocities close to 0.


%Consider errors in slope and impact this has:
%Vaze, J., Teng, J., & Spencer, G. (2010). Impact of DEM accuracy and 
%resolution on topographic indices. Environmental Modelling & Software, 
%25(10), 1086–1098.
