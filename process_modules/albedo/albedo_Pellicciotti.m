function varargout = albedo_Pellicciotti(varargin)

%See albedo discussion in:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm

if isempty(varargin(:))
    varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'albedo_fresh', 0.5, 1, 0.71, 'albedo_Pellicciotti','cryo'});
%     varargout{1} = cat(1,varargout{1}, {'albedo_decay', 0.0, 0.5, 0.11, 'albedo_Pellicciotti','cryo'});
%     varargout{1} = cat(1,varargout{1}, {  'albedo_ice', 0.2, 0.7, 0.45, 'albedo_Pellicciotti','cryo'});
    
    return
else
%     aFresh = find_att(varargin{1}.coef,'albedo_fresh'); %Albedo of fresh, deep snow
%     aDecay = find_att(varargin{1}.coef,'albedo_decay');  
%     aIce = find_att(varargin{1}.coef,'albedo_ice');
end

aIce = find_att(varargin{1}.global,'albedo_ice');

aFresh = find_att(varargin{1}.global,'albedo_snow_fresh'); %Albedo of fresh snow
aOld   = find_att(varargin{1}.global,'albedo_snow_old');  %Albedo of old/melting snow


if ~isfield(sAtm,'tasmaxc')
    sAtm.tasmaxc = zeros(size(sCryo.snw),'single');
end

%Update accumulated daily maximum air temperature since last snow event (used in albedo calculation):
sAtm.tasmaxc = sAtm.tasmaxc + squeeze(sAtm.tasmax(sAtm.indCurr,:,:));
sAtm.tasmaxc( sAtm.prsn > 0 ) = 0;
sAtm.tasmaxc(sAtm.tasmaxc < 0) = 0;

%Albedo calculation:
sCryo.snalb = real(aFresh - aOld*log10(sAtm.tasmaxc));

%Set Maximum albedo to that of fresh snow
sCryo.snalb(sCryo.snalb > 1) = aFresh;
sCryo.snalb( isnan(sCryo.snalb) | sCryo.snalb < aIce) = aIce; 
    %This value shouldn't actually matter, because seperate albedo used for ice
    
sCryo.icalb = aIce*ones(size(sCryo.snalb));