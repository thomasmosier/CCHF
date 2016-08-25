function varargout = heat_ETI_Hock(varargin)
%See reference:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'tas_watt_const',  0, 50, 1, 'heat_ETI_Hock','cryo'}); %W/m^2
    varargout{1} = cat(1,varargout{1}, { 'raad_snow_pwr', -2,   2,  0, 'heat_ETI_Hock','cryo'}); %Unitless
    varargout{1} = cat(1,varargout{1}, {  'raad_ice_pwr', -2,   2,  0, 'heat_ETI_Hock','cryo'}); %Unitless
%     varargout{1} = cat(1,varargout{1}, {'shortwave_melt', 0, 0.2, 0.002, 'ETI_Hock','cryo'});
        
    return
else
    sMeta  = varargin{1};
    wCon   = find_att(varargin{1}.coef, 'tas_watt_const');  
    mSnPwr = find_att(varargin{1}.coef,  'raad_snow_pwr'); 
    mIcPwr = find_att(varargin{1}.coef,   'raad_ice_pwr'); 
%     srf = find_att(varargin{1}.coef,'shortwave_melt'); 
end

%Units of J-m^{-2}-hr-mm^{-1}-s^{-1} := W-m^{-2}-hr-mm^{-1}

% %Determine number of time steps per day
% if regexpbl(varargin{1}.dt,'month')
%     nStep = 1/eomday(varargin{1}.dateCurr(1), varargin{1}.dateCurr(2));
% elseif regexpbl(varargin{1}.dt,{'day','daily'})
%     nStep = 1;
% elseif regexpbl(varargin{1}.dt,'hour')
%     nStep = 24;
% else
%     error('groundwater_bucket:dtUnknown',[varargin{1}.dt ' is an unknown time-step.'])
% end


%Find top of the atmosphere radiation:
indRSDT = find(ismember(sAtm.datersdt(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end

%Calculate melt potential using Hock's formulation: 
%Because of conversion factor, LHS has units of mm/hr
sCryo.hfnet = squeeze(sAtm.tas(sAtm.indCurr,:,:)).* ...
	(wCon + (10^mSnPwr)*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));
sCryo.hfneti = squeeze(sAtm.tas(sAtm.indCurr,:,:)).* ...
	(wCon + (10^mIcPwr)*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));


% %Convert from units of mm/hr to W/m2
% densW = find_att(sMeta.global,'density_water'); 
% cLate = densW*find_att(sMeta.global,'latent_water');  
% 
% convFact = cLate*(1/1000)*(1/3600); %Latent fusion * length conver * time conver; units = J-m-hr-m^-3-mm^-1-s^-1
% sCryo.hfnet = convFact*sCryo.hfnet;


% sCryo.hfnet = squeeze(sAtm.tas(sAtm.indCurr,:,:)).* ...
% 	(wConst + 10^(srf)*(1-sCryo.snalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));
% sCryo.hfneti = squeeze(sAtm.tas(sAtm.indCurr,:,:)).* ...
% 	(wConst + 10^(srf)*(1-sCryo.icalb).*sAtm.rstran.*squeeze(sAtm.rsdt(indRSDT,:,:)));
