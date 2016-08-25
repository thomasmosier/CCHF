function varargout = flux_precip(varargin)

global sCryo sAtm

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    return
else
    sMeta = varargin{1};
end
    
%Find seconds in current time-step:
dt = time2sec(1,sMeta.dt,sMeta.dateCurr);

lWat = find_att(sMeta.global, 'heat_cap_water')*find_att(sMeta.global, 'density_water');

%Sensible heat (units of Watts) = density of water (1000 kg/m^3) * specific heat capacity water (4179 J/kg/Deg C) * m of precipitation * (airTemp - snowTemp)
if isfield(sCryo,'tss')
    sCryo.hfcp = (lWat/dt)*squeeze(sAtm.pr(sAtm.indCurr,:,:)) ...
        .*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) - sCryo.tsn);
elseif isfield(sCryo,'tsn')
    sCryo.hfcp = (lWat/dt)*squeeze(sAtm.pr(sAtm.indCurr,:,:)) ...
        .*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) - sCryo.tsn);
else
    error('heat_precip:missingTmp',['Neither internal or surface '...
        'snow temperature are known.']);
end

%%LATENT HEAT IS DEALTH WITH IN ENERGY FUNCTION!!!
% %latent heat (units of Watts) = density of water (1000 kg/m^3) * latent heat of fusion of water (334000 J/kg) * m of rain
% %Add latent heat of rain impinging on snowpack
% indFreeze = find(sCryo.tsn < 0);
% sCryo.hfcp(indFreeze) = sCryo.hfcp(indFreeze) + (334000000/dt)*sAtm.rain(indFreeze);
%     sCryo.hfcp(isnan(sCryo.hfcp)) = 0;