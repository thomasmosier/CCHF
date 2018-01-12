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

function varargout = ETI_LST(sHydro, varargin)
%See reference:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm sLand


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {  'Watt_per_deg',  0, 35, 20, 'ETI_LST','cryo'});
    varargout{1} = cat(1,varargout{1}, heat_long_out());
    varargout{1} = cat(1,varargout{1}, heat_conduct_c());
    
%     varargout{1} = cat(1,varargout{1}, {'shortwave_melt', 0, 0.2, 0.002, 'ETI_Pellicciotti','cryo'});
        
    return
else
    sMeta = varargin{1};
    wattperdeg = find_att(varargin{1}.coef,'Watt_per_deg'); 
%     srf = find_att(varargin{1}.coef,'shortwave_melt'); 
end

% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
%Units of J-m^{-2}-hr-mm^{-1}-s^{-1} := W-m^{-2}-hr-mm^{-1} (because [dis] = mm/hr/deg 


%Find top of the atmosphere radiation:
indRSDT = find(ismember(sLand.rsdtdate(:,2:end),sMeta.currTime(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end
    

%Shortwave melt energy:
sCryo.heatSR = (1-sCryo.albedoS).*sAtm.rstran.*squeeze(sLand.rsdt(indRSDT,:,:));
sCryo.heatSRI = (1-sCryo.albedoI).*sAtm.rstran.*squeeze(sLand.rsdt(indRSDT,:,:));

%Longwave radiation:
heat_long_simple(varargin{1});

%%Temperature melt energy:
sCryo.heatT = wattperdeg*squeeze(sAtm.tas(sAtm.indCurr,:,:));

% heat_conduct_c(varargin{1});

%Calculate melt potential using Pellicciotti's formulation: 
%Because of conversion factor, each term has units of w/m^2
sCryo.heat  = sCryo.heatT +  sCryo.heatSR + sCryo.heatLR;
sCryo.heatI = sCryo.heatT + sCryo.heatSRI + sCryo.heatLR;
% sCryo.heat  = sCryo.heatT +  sCryo.heatSR + 10^(lwPwr)*sCryo.heatLR + sCryo.heatPC;
% sCryo.heatI = sCryo.heatT + sCryo.heatSRI + 10^(lwPwr)*sCryo.heatLR + sCryo.heatPC;

%VERSION WITH RAAD SCALING (REDUNDENT):
% sCryo.heat = dis*fConvert*squeeze(sAtm.tas(sAtm.indCurr,:,:)) ...
%     + srf*(1-sCryo.albedoS).*sAtm.rstran.*squeeze(sLand.rsdt(indRSDT,:,:));