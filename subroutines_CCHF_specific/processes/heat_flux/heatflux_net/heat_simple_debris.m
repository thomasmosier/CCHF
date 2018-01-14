function varargout = heat_simple_debris(varargin)
%Notes:
%Allows ice melt if 'ice' field present
%Same melt rate for snow and ice 
%No albedo factor

%See reference on formulation and various degre index factors:
%Singh, P., Kumar, N., & Arora, M. (2000). Degree–day factors for snow and 
%ice for Dokriani Glacier, Garhwal Himalayas. Journal of Hydrology, 235(1), 
%1-11.

global sCryo sAtm


if isempty(varargin(:))
    varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {    'watt_per_deg',   0, 60, 42.35, 'heat_simple_degree','cryo'});
    
    return
else
    wattperdegS = find_att(varargin{1}.coef,'watt_per_deg');
end


% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
% %Units of J-hr-m^{-2}-s^{-1}

%Calculate equivalent of 'heatflux' using simple degree index model
sCryo.hfnet = wattperdegS*squeeze(sAtm.tas(sAtm.indtas,:,:)); %units to Watts per m^2

%Heat flux for ice
if isfield(sCryo, 'icedbr') %If debris cover information available
    relMelt = debris_melt_emp(sCryo.icedbr);
    
    sCryo.hfneti = relMelt.*sCryo.hfnet; 
else
    sCryo.hfneti = sCryo.hfnet;
    warning('heat_simple_debris:noDebris','This function is designed to work with a debris cover layer, which is not present.');
end

