function varargout = heat_LST_Mosier(varargin)
%See reference:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {  'Watt_per_deg',  0, 35, 10, 'LSD_Mosier','cryo'});
    varargout{1} = cat(1,varargout{1}, flux_long_simple());
%     varargout{1} = cat(1,varargout{1}, {'sw_pwr', -2, 2, -0.43, 'LSD_Mosier','cryo'}); 
%     varargout{1} = cat(1,varargout{1}, {'lw_pwr', -2, 2, 0, 'LSD_Mosier', 'cryo'});
    
    return
else
    sMeta = varargin{1};
    wattperdeg = find_att(varargin{1}.coef,'Watt_per_deg'); 
%     srf = find_att(varargin{1}.coef,'SW_pwr'); 
%     lrf = find_att(varargin{1}.coef,'lw_pwr'); 
end

% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
%Units of J-m^{-2}-hr-mm^{-1}-s^{-1} := W-m^{-2}-hr-mm^{-1} (because [dis] = mm/hr/deg 

if ~isfield(sCryo, 'tsn')
   sCryo.tsn = zeros(size(sCryo.snw), 'single'); 
end  

%Find top of the atmosphere radiation:
indRSDT = find(ismember(sAtm.datersdt(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end
    

%Shortwave melt energy:
if isfield(sAtm,'rsdf') %Diffuse shortwave present ('rsdf'):
    sCryo.hfrs  = (1-sCryo.snalb).*((1-sAtm.rstran).*squeeze(sAtm.rsdf(indRSDT,:,:)) + ...
        sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));
    sCryo.hfrsi = (1-sCryo.icalb).*((1-sAtm.rstran).*squeeze(sAtm.rsdf(indRSDT,:,:)) + ...
        sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));
else %Direct radiation only
    sCryo.hfrs  = (1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
    sCryo.hfrsi = (1-sCryo.icalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));
end
    

%Longwave radiation (dependent on surface temperature):
if isfield(sCryo,'tssn') && isfield(sCryo,'tsic')
    sCryo.hfrl = flux_long_simple(sCryo.tssn, varargin{1});
    sCryo.hfrli = flux_long_simple(sCryo.tsic, varargin{1});
elseif isfield(sCryo,'tsn')
    sCryo.hfrl = flux_long_simple(sCryo.tsn, varargin{1});
    sCryo.hfrli = sCryo.hfrl;
else
    error('heat_LST_Mosier:noSnTmp',['A snow temperature field is '...
        'needed for the longwave readiation function but is not available.']);
end
% sCryo.hfrl = 10^(lrf)*sCryo.hfrl;

%%Temperature melt energy:
sCryo.hft = wattperdeg*squeeze(sAtm.tas(sAtm.indCurr,:,:));

%Calculate melt potential using Pellicciotti's formulation: 
%Because of conversion factor, each term has units of w/m^2
sCryo.hfnet  = sCryo.hfrs  + sCryo.hft + sCryo.hfrl;
sCryo.hfneti = sCryo.hfrsi + sCryo.hft + sCryo.hfrli;
% sCryo.hfnet  = sCryo.hft +  sCryo.hfrs + 10^(lwPwr)*sCryo.hfrl + sCryo.hfnetPC;
% sCryo.hfneti = sCryo.hft + sCryo.hfrsi + 10^(lwPwr)*sCryo.hfrl + sCryo.hfnetPC;

%VERSION WITH RAAD SCALING (REDUNDENT):
% sCryo.hfnet = dis*fConvert*squeeze(sAtm.tas(sAtm.indCurr,:,:)) ...
%     + srf*(1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:));