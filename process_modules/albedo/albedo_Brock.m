function varargout = albedo_Brock(varargin)

%Primary reference:
%Brock, B. W., Willis, I. C., & Sharp, M. J. (2000). Measurement and 
%parameterization of albedo variations at Haut Glacier d’Arolla, 
%Switzerland. Journal of Glaciology, 46(155), 675–688.

%See albedo discussion in:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm



if isempty(varargin(:))
    varargout{1} = snow_albedo_Brock(0.1, 0.8, 0.5); %Doesn't matter because this is only for parameter mode
    return
end

%NO FITTING PARAMETERS
aUnder = find_att(varargin{1}.global,'albedo_ground'); %Underlying or debris albedo (Unitless)
aIce = find_att(varargin{1}.global,'albedo_ice');

aFresh = find_att(varargin{1}.global,'albedo_snow_fresh'); %Albedo of fresh snow
aOld   = find_att(varargin{1}.global,'albedo_snow_old');  %Albedo of old/melting snow

%FITTING PARAMETERS:
% if isempty(varargin(:))
%     argout = snow_albedo_Brock(sSnowpack, aUnder); 
%     argout(end+1,:) = {'albedo_ice', 0.2, 0.7, 0.45, 'albedo_Brock'}; 
%     return
% else  
%     aIce = find_att(varargin{1}.coef,'albedo_ice');
% end

if ~isfield(sAtm,'tasmaxc')
    sAtm.tasmaxc = zeros(size(sCryo.snw),'single');
end

%Update accumulated daily maximum air temperature since last snow event (used in albedo calculation):
sAtm.tasmaxc = sAtm.tasmaxc + squeeze(sAtm.tasmax(sAtm.indCurr,:,:));
sAtm.tasmaxc( sAtm.prsn > 0 ) = 0;
sAtm.tasmaxc(sAtm.tasmaxc < 0) = 0;

snow_albedo_Brock(aUnder, aFresh, aOld, varargin{1});

%For now, set constant albedo of ice:
sCryo.icalb = aIce*ones(size(sCryo.snalb));