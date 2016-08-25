function varargout = heat_ETI_Pellicciotti(varargin)
%See reference:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'Watt_per_deg',  0, 30, 13.57, 'heat_ETI_Pellicciotti','cryo'});
    varargout{1} = cat(1,varargout{1}, {'SW_pwr',       -1, 1, -0.43, 'heat_ETI_Pellicciotti','cryo'}); 
        
    return
else
    sMeta = varargin{1};
    wattperdeg = find_att(varargin{1}.coef,'Watt_per_deg');  
    srf = find_att(varargin{1}.coef,'SW_pwr'); 
end

% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
%Units of J-m^{-2}-hr-mm^{-1}-s^{-1} := W-m^{-2}-hr-mm^{-1} (because [dis] = mm/hr/deg 


%Find top of the atmosphere radiation:
indRSDT = find(ismember(sAtm.datersdt(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end
    


%Shortwave melt energy:
% if isfield(sLand,'rsdf') %Direct and diffuse shortwave:
%     sCryo.hfrs  = (1-sCryo.snalb).*((1-sAtm.rstran).*squeeze(sLand.rsdf(indRSDT,:,:)) + ...
%         sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));
%     sCryo.hfrsi = (1-sCryo.icalb).*((1-sAtm.rstran).*squeeze(sLand.rsdf(indRSDT,:,:)) + ...
%         sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));
% else %Direct radiation only
%     sCryo.hfrs  = (1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
%     sCryo.hfrsi = (1-sCryo.icalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
% end
sCryo.hfrs  = 10^(srf)*(1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
sCryo.hfrsi = 10^(srf)*(1-sCryo.icalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));

%Temperature melt energy:
sCryo.hft = wattperdeg*squeeze(sAtm.tas(sAtm.indCurr,:,:));

%Calculate melt potential using Pellicciotti's formulation: 
%Because of conversion factor, each term has units of w/m^2
sCryo.hfnet = sCryo.hft + sCryo.hfrs;
sCryo.hfneti = sCryo.hft + sCryo.hfrsi;