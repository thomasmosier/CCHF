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

function varargout = density_Cosima(varargin)

%See densification discussion in:
%Huintjes, E., Sauter, T., Schröter, B., Maussion, F., Yang, W., Kropá?ek, 
%J., ? Schneider, C. (2015). Evaluation of a Coupled Snow and Energy 
%Balance Model for Zhadang Glacier, Tibetan Plateau, Using Glaciological 
%Measurements and Time-Lapse Photography. Arctic, Antarctic, and Alpine 
%Research, 47(3), 573?590. htsCryo.tsn://doi.org/10.1657/AAAR0014-073


global sCryo

if isempty(varargin(:))
    varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'albedo_fresh', 0.5, 1, 0.71, 'density_Cosima','cryo'});

    
    return
else
%     aFresh = find_att(varargin{1}.coef,'albedo_fresh'); %Albedo of fresh, deep snow
end


rhoI = find_att(varargin{1}.global, 'density_ice');
rhoW = find_att(varargin{1}.global, 'density_water');
R = find_att(varargin{1}.global, 'gas_constant_univ');

K0   = 11;     % rate factors [-]
K1   = 575;
E0   = 10260;  % activation energy
E1   = 21400;
%Assume density of new snow is 250
rhoN = 250;


if ~isfield(sCryo,'rhosn')
    sCryo.rhosn = zeros(size(sCryo.snw), 'single');
end

if ~isfield(sCryo,'lwsnn')
    sCryo.lwsnn = zeros(size(sCryo.snw), 'single');
end



%DENSITY IMPACT OF MELT WATER:
%This is H_new/H_tot = SWE_new/rho_new / (SWE_new/rho_new +
%SWE_old/rho_old) = 1 / (1 + (rho_new*SWE_old)/(rho_old*SWE_new))
fracNM = 1./(1 + (rhoW*sCryo.snw)./(sCryo.lwsnn.*sCryo.rhosn));

%Calculate new density based on fresh melt and snowpack
sCryo.rhosn = fracNM*rhoW + (1 - fracNM).*sCryo.rhosn;

%Try to capture amount of new melt:
sCryo.lwsnn = abs(sCryo.lwsnl - sCryo.lwsnn - sCryo.snlr);
    

%DENSITY OF NEW SNOW:
%Calculate the fractional height of new snow versus old snow:
%This is H_new/H_tot = SWE_new/rho_new / (SWE_new/rho_new +
%SWE_old/rho_old) = 1 / (1 + (rho_new*SWE_old)/(rho_old*SWE_new))
fracN = 1./(1 + (rhoN*(sCryo.snw - sAtm.prsn))./(sAtm.prsn.*sCryo.rhosn));

%Calculate new density based on fresh snowfall and old snowpack
sCryo.rhosn = fracN*250 + (1 - fracN).*sCryo.rhosn;


%DENSIFICATION:
%Account for density changes corresponding to melt water and liquid drained
%from snowpack:
%Net liquid content change is (melt - refreeze) - snlr

%See Eq. 15 in reference:
indL550 = find(sCryo.rhosn < 550);
indG550 = setdiff((1:numel(sCryo.snw)),indL550);

sCryo.rhosn(indL550) = sCryo.rhosn(indL550) ...
    + K0*exp(-E0./(R.*sCryo.tsn(indL550))).*1000*sAtm.prsn(indL550).*((rhoI-sCryo.rhosn(indL550))./rhoI);
sCryo.rhosn(indG550) = sCryo.rhosn(indG550) ...
    + K1.*exp(-E1./(R.*sCryo.tsn(indG550))).*1000*sAtm.prsn(indG550).*((rhoI-sCryo.rhosn(indG550))./rhoI);

sCryo.rhosn(sCryo.rhosn > rhoI) = rhoI;