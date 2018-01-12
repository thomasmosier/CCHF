function varargout = snow_albedo_Brock( aUnder, aFresh, aOld, varargin)
global sCryo
%Formulation from:
%Brock, B. W., Willis, I. C., & Sharp, M. J. (2000). Measurement and 
%parameterization of albedo variations at Haut Glacier d'Arolla, 
%Switzerland. Journal of Glaciology, 46(155), 675-688.

%a_s = (1-exp[d/d*])*a_ds + exp[-d/d*]*a_ss
    %a_ds = 0.713 - 0.112*Log10[T_a]
    %a_ss = a_u + 0.442*exp[-0.058*T_a]
    
    %a_u = underlying ice or debris albedo
    

%WITHOUT FITTING PARAMETERS:    
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    
    return
end

    %From Brock paper:
%     aFresh = 0.71; %Albedo of fresh, deep snow
%     aOld = 0.11;  
    
%WITH FITTING PARAMETERS:
% if isempty(varargin(:))
% 	argout = cell(2,5);
%     
%     argout(1,:) = {   'albedo_fresh', 0.5, 1, 0.71, 'snow_albedo_Brock'}; 
%     argout(2,:) = {   'albedo_decay', 0.0, 0.5, 0.11, 'snow_albedo_Brock'}; 
%     
%     return
% else
%     aFresh = find_att(varargin{1}.coef,'albedo_fresh'); %Albedo of fresh, deep snow
%     aDecay = find_att(varargin{1}.coef,'albedo_decay');  
% end

dChar = 0.024; %Characteristic depth of snow, where all albedo property of snow (Units of meters)

%Albedo calculation:
sCryo.snalb = real((1-exp(-sCryo.snw/dChar)).*(aFresh-aOld*log10(sCryo.tasmaxCum))...
    +exp(-sCryo.snw/dChar).*(aUnder + 0.442*exp(-0.058*sCryo.tasmaxCum)));

%Set Maximum albedo to that of fresh snow
sCryo.snalb(sCryo.snalb > 1) = aFresh;
sCryo.snalb( isnan(sCryo.snalb) | sCryo.snalb < aUnder) = aUnder;    