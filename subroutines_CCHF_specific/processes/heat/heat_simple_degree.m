function varargout = heat_simple_degree(varargin)
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
    varargout{1} = cat(1,varargout{1}, {'watt_per_deg_snow', 10, 40, 20, 'heat_simple_degree','cryo'});
    varargout{1} = cat(1,varargout{1}, {'watt_per_deg_ice',  10, 60, 20, 'heat_simple_degree','cryo'});
    varargout{1} = cat(1,varargout{1}, {  'tas_offset',  0,   4,     1, 'heat_simple_degree','cryo'}); %Units of depth melt
    
    %The degree index factor here is approx 3.9 times degree-day factors
    %with units of mm / C / day. Shea et al., 2015 find snow DD factor of
    %about 5 mm / C / day, which is approximately 19 - 20, and clean ice DD
    %factor of 9.7, which is 30-40 (range comes from approx standard
    %deviation cites in their paper)

    
%     if isfield(sCryo, 'icedbr')
%         varargout{1} = cat(1,varargout{1}, {    'watt_per_deg_debris',   0, 60, 20, 'heat_simple_degree','cryo'});
%     end
    return
else
    wattperdegS = find_att(varargin{1}.coef,'watt_per_deg_snow'); 
    wattperdegI = find_att(varargin{1}.coef,'watt_per_deg_ice'); 
    tasOffset = find_att(varargin{1}.coef, 'tas_offset');
%     if isfield(sCryo, 'icedbr')
%         wattperdegDeb = find_att(varargin{1}.coef,'watt_per_deg_debris'); 
%     end
end


% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
% %Units of J-hr-m^{-2}-s^{-1}

%Calculate equivalent of 'heatflux' using simple degree index model
sCryo.hfnet = wattperdegS*(squeeze(sAtm.tas(sAtm.indtas,:,:)) - tasOffset); %units to Watts per m^2
sCryo.hfneti = wattperdegI*(squeeze(sAtm.tas(sAtm.indtas,:,:)) - tasOffset); %units to Watts per m^2

% %Heat flux for ice
% if isfield(sCryo, 'icedbr') %If debris cover information available
%     sCryo.hfneti = wattperdegDeb*(squeeze(sAtm.tas(sAtm.indtas,:,:)) - tasOffset); 
% 
%     if ~isfield(sCryo, 'indicecln')
%         debThresh = find_att(varargin{1}.global,'debris_threshold');
%         sCryo.indicecln = find(sCryo.icedbr < debThresh);
%     end
%     
%     if ~isempty(sCryo.indicecln)
%         sCryo.hfneti(sCryo.indicecln) = sCryo.hfnet(sCryo.indicecln);
%     end
% else
%     sCryo.hfneti = sCryo.hfnet;
%     %sCryo.hfneti = wattperdegI*squeeze(sAtm.tas(sAtm.indCurr,:,:)); %units to Watts per m^2
% end

