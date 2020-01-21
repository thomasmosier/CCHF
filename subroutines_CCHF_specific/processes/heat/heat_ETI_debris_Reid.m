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

function varargout = heat_ETI_debris_Reid(varargin)
%Debris portion inspired by:
%Reid, T. D., & Brock, B. W. (2010). An energy-balance model for 
%debris-covered glaciers including heat conduction through the debris 
%layer. Journal of Glaciology, 56(199), 903?916.

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
    wattperdegI  = find_att(varargin{1}.coef, 'Watt_per_deg_ice'); 
end

% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
%Units of J-m^{-2}-hr-mm^{-1}-s^{-1} := W-m^{-2}-hr-mm^{-1} (because [dis] = mm/hr/deg 


%Find top of the atmosphere radiation:
indRSDT = find(ismember(sAtm.datersdt(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end
   
%Initialization temperature for ice:
tmpIce = 0;

%Calculate shortwave radiation for snow and ice:
sCryo.hfrs  = (1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
sCryo.hfrsi = (1-sCryo.icalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));

%Temperature sensible energy (parameterization for convenction):
sCryo.hft  = wattperdegS*squeeze(sAtm.tas(sAtm.indtas,:,:));
sCryo.hfti = wattperdegI*squeeze(sAtm.tas(sAtm.indtas,:,:));

%Calculate snow melt potential using Pellicciotti's formulation: 
%Because of conversion factor, each term has units of w/m^2
sCryo.hfnet = sCryo.hft + sCryo.hfrs;

%Initialize ice heatflux:
sCryo.hfneti = nan(size(sCryo.hfnet), 'single');

%%Modify heat flux based on debris cover / clean ice
%Find indices of debris covered and clean ice:
if isfield(sCryo, 'icdbr') %If debris cover information available
    %Find debris-free and clean ice indices:
    indCln = find(sCryo.icdbr == 0);
    indDeb = find(sCryo.icdbr > 0);
else
    indCln = (1:numel(sCryo.hfnet));
    indDeb = [];
    
    warning('heatEtiDebris:noDebrisField','No debris field is present. This heat representation is designed for use with a debris field.')
end


%Find heat propogation through debris
%(necessary for determining heat into ice at bottom of debris because 
%this depends on the local-gradient)
if ~isempty(indDeb)    
    %Define values used in debris conductivity. Values come from: 
    %Rounce, D. R., & McKinney, D. C. (2014). Debris thickness of glaciers in 
    %the Everest area (Nepal Himalaya) derived from satellite imagery using a 
    %nonlinear energy balance model.
    rhoD = find_att(sMeta.global, 'density_debris'); 
    cpD  = find_att(sMeta.global, 'heat_cap_debris'); 
    kD   = find_att(sMeta.global, 'thermal_conduct_debris'); 
    
    %Define constant used in heat propogation:
    %C = k*deltaT/(2*rho*cp*deltaZ^2);
    Cf = kD*time2sec(1,sMeta.dt,sMeta.dateCurr)/(2*rhoD*cpD);
    
    %Create grid that stores the number of layers used to model debris heat
    %transfer at each grid cell
    if ~isfield(sCryo,'icdbrlyr')
        sCryo.icdbrlyr = sparse(numel(sCryo.icx), 1); %Uses linear indexing
        sCryo.icdbrdz  = sparse(numel(sCryo.icx), 1); 
        for ii = 1 : numel(indDeb)
            thickLayer = 0.01; %Each layer should be approximately one 1 cm thick
            
            %Use N+1 layers because they are actually "skin" temperatures
            %at boundaries between layers
            sCryo.icdbrlyr(indDeb(ii)) = max(round(sCryo.icdbr(indDeb(ii))/thickLayer), 1) + 1;
            
            sCryo.icdbrdz(indDeb(ii)) = sCryo.icdbr(indDeb(ii))/(sCryo.icdbrlyr(indDeb(ii))-1);
        end
    end
    
    %Initialize grid to hold temperature at boundary of each layer (=N
    %because bottom temperature is constant for all grid cells)
    if ~isfield(sCryo, 'icdbrtmp')
        %Use 2-D sparse array (using class from Matlab FileExchange)
        sCryo.icdbrtmp = sparse(numel(sCryo.icx), full(max(sCryo.icdbrlyr(:)))); %1st dimension is linear index of geographic domain, 2nd index is number of layers
        
        %Initialize debris temperatures to ice tmp (requires spin-up):
        for ii = 1 : numel(indDeb)
        	sCryo.icdbrtmp(indDeb(ii), 1 : sCryo.icdbrlyr(indDeb(ii))) = tmpIce;
        end
    end
    
    %Initialize field to store previous downward heatflux over debris
    if ~isfield(sCryo, 'icdbrq')
        sCryo.icdbrq = sparse(numel(sCryo.icx), 1);
    end
    
    %Enforce boundary conditions (surface heat flux = 'Neumann' and temperature
    %at base = 'dirchelet'):
    %Upper boundary condition:
    icdbrbc1 = (Cf./(sCryo.icdbr(indDeb).^2)).*full(2*sCryo.icdbrdz(indDeb)/kD)...
        .*(sCryo.hfrsi(indDeb) + sCryo.hfti(indDeb) + full(sCryo.icdbrq(indDeb)));
    %lower boundary condition (commented out becuase it is enforces in diffusion eqn):
%     for ii = 1 : numel(indDeb)
%         
%         sCryo.icdbrtmp(indDeb(ii), sCryo.icdbrlyr(indDeb(ii))) = tmpIce;
%     end
    
    %Update previous downward heatflux over debris for use in next
    %iteration:
    sCryo.icdbrq(indDeb) = sCryo.hfrsi(indDeb) + sCryo.hfti(indDeb);

    %Propogate heat through debris:
    for ii = 1 : numel(indDeb)
        sCryo.icdbrtmp(indDeb(ii),1:sCryo.icdbrlyr(indDeb(ii))) ...
            = heat_diffusion_CN(full(sCryo.icdbrtmp(indDeb(ii),1:sCryo.icdbrlyr(indDeb(ii)))), 'neumann', icdbrbc1(ii), 'dirichlet', tmpIce, Cf/(sCryo.icdbr(indDeb(ii))^2));
    end
    
    %For testing
%     sCryo.icdbrtmp(indDeb(2),1:sCryo.icdbrlyr(indDeb(2)))

    %Find heat at ice surface (under debris):
    sCryo.hfneti(indDeb) = kD*(full(sCryo.icdbrtmp(sub2ind(size(sCryo.icdbrtmp), indDeb, sCryo.icdbrlyr(indDeb)-1))) - tmpIce*ones(numel(indDeb), 1)) ...
        ./sCryo.icdbrdz(indDeb);
end

%Calculate melt for clean ice (no conduction)
if ~isempty(indCln)
    sCryo.hfneti(indCln) = sCryo.hft(indCln) + sCryo.hfrsi(indCln); 
end


%Melt is enhanced by factor of 10 at locations that are entirely ponds.
%Use fractional relation to account for partially ponded grid cells
%(when fraction is 0, there is no impact; when fraction is 1, the
%factor is 10)
fldLake = '';
if isfield(sCryo, 'icpndx')
    fldLake = 'icpndx';
elseif isfield(sCryo, 'iclk')
    fldLake = 'iclk';
end

if ~isempty(fldLake)
    sCryo.hfneti = (1+9*sCryo.(fldLake)).*sCryo.hfneti;
end


sCryo.hfnet(isnan(sCryo.icx)) = nan;
sCryo.hfneti(isnan(sCryo.icx)) = nan;