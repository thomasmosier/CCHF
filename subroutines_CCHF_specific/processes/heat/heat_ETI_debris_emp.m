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

function varargout = heat_ETI_debris_emp(varargin)
%Empricial debris formulation inspired by:
%Kraaijenbrink, P. D. A., Bierkens, M. F. P., Lutz, A. F., &
%Immerzeel, W. W. (2017). Impact of a global temperature rise of 1.5
%degrees Celsius on Asia?s glaciers. Nature, 549(7671), 257?260.
%https://doi.org/10.1038/nature23878

%Enhanced Temperature Index formation inspired by Pellicciotti:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'Watt_per_deg_snow', 5, 40, 10, 'heat_ETI_debris','cryo'});
    varargout{1} = cat(1,varargout{1}, {'Watt_per_deg_ice',  5, 50, 20, 'heat_ETI_debris','cryo'});
     
    %For simple degree index (i.e. no solar radiation)
    %The degree index factor here is approx 3.9 times degree-day factors
    %with units of mm / C / day. Shea et al., 2015 find snow DD factor of
    %about 5 mm / C / day, which is approximately 19 - 20, and clean ice DD
    %factor of 9.7, which is 30-40 (range comes from approx standard
    %deviation cites in their paper)
    
    return
else
    sMeta = varargin{1};
    wattperdegS = find_att(varargin{1}.coef, 'Watt_per_deg_snow');  
    wattperdegI = find_att(varargin{1}.coef, 'Watt_per_deg_ice'); 
end

% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
%Units of J-m^{-2}-hr-mm^{-1}-s^{-1} := W-m^{-2}-hr-mm^{-1} (because [dis] = mm/hr/deg 


%Find top of the atmosphere radiation:
indRSDT = find(ismember(sAtm.datersdt(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end


%Calculate shortwave radiation for snow and ice:
sCryo.hfrs  = (1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
sCryo.hfrsi = (1-sCryo.icalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));


%Temperature sensible energy (parameterization for convenction):
sCryo.hft  = wattperdegS*squeeze(sAtm.tas(sAtm.indtas,:,:));
sCryo.hfti = wattperdegI*squeeze(sAtm.tas(sAtm.indtas,:,:));
    sCryo.hft(isnan(sCryo.icx)) = nan;
    sCryo.hfti(isnan(sCryo.icx)) = nan;

%Calculate snow melt potential using Pellicciotti's formulation: 
%Because of conversion factor, each term has units of w/m^2
sCryo.hfnet = sCryo.hft + sCryo.hfrs;
%Calculate ice heatflux (based on clean ice; debris is modified below):
sCryo.hfneti = sCryo.hft + sCryo.hfrsi; 

sCryo.hfnet(isnan(sCryo.icx)) = nan;
sCryo.hfneti(isnan(sCryo.icx)) = nan;


%%Modify heat flux based on debris cover / clean ice
%Inspired by Kraaijenbrink
if isfield(sCryo, 'icdbr') %If debris cover information available
    if ~isfield(sCryo, 'icdbrmelt')
        sCryo.icdbrmelt = debris_melt_empirical(sCryo.icdbr, 'm');
    end

    sCryo.hfneti = sCryo.icdbrmelt.*sCryo.hfneti;
end


%Melt is enhanced by factor of 10 at locations that are entirely ponds.
%Use fractional relation to account for partially ponded grid cells
%(when fraction is 0, there is no impact; when fraction is 1, the
%factor is 10)
%Inspired by Kraaijenbrink
fldLake = '';
if isfield(sCryo, 'icpndx')
    fldLake = 'icpndx';
elseif isfield(sCryo, 'iclk')
    fldLake = 'iclk';
end

if ~isempty(fldLake)
    sCryo.hfneti = (1+9*sCryo.(fldLake)).*sCryo.hfneti;
end
