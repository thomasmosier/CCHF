function varargout = heat_sensible_simple(varargin)

%Formulation taken from pages 161-164 of 
%DeWalle, D. R., & Rango, A. (2008). Principles of snow hydrology. 
%Cambridge University Press.

global sCryo sAtm sLand

%TWO FITTING PARAMETERS:
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {  'stable_wind', 0, 10,   1, 'heat_sensible_regimes','cryo'}); 
    varargout{1} = cat(1,varargout{1}, {'unstable_wind', 0, 10, 2.5, 'heat_sensible_regimes','cryo'}); 
    varargout{1} = cat(1,varargout{1}, ['tsn',cell(1,5)]);
    return
else
    uStab = find_att(varargin{1}.coef,  'stable_wind'); 
    uUnst = find_att(varargin{1}.coef,'unstable_wind'); 
end

riStab = 0.8; 
riUnst = -1; 


%MANY FITTING PARAMETERS:
% if isempty(varargin(:))
% 	argout = cell(4,5);
%     argout(1,:) = {      'ri_stable', 0.3,  0.1,  0.8, 'heat_sensible_regimes'};
%     argout(2,:) = {    'ri_unstable', -10, -0.1,   -1, 'heat_sensible_regimes'}; 
%     argout(3,:) = {  'scalar_stable',   0,    10,    1, 'heat_sensible_regimes'}; 
%     argout(4,:) = {'scalar_unstable',   0,    10,    1, 'heat_sensible_regimes'}; 
%     return
% else
%     riStab = find_att(varargin{1}.coef,      'ri_stable'); 
%     riUnst = find_att(varargin{1}.coef,    'ri_unstable'); 
%     scStab = find_att(varargin{1}.coef,  'scalar_stable'); 
%     scUnst = find_att(varargin{1}.coef,'scalar_unstable'); 
% end

%T_air - T_snow
%  Stable when  T_air >> T_snow (i.e. difference positive)
%Unstable when T_snow >>  T_air (i.e. difference negative)
tempDiff = squeeze(sAtm.tas(sAtm.indCurr,:,:)) - sCryo.tsn;
tempAvg = 0.5*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) + sCryo.tsn) + 273.15; %Units of Kelvin

%Richardson Number:
%g*Z_air*(T_air-T_snow)/(T_avg * u_air^2)
indPos = find(tempDiff >= 0);
indNeg = find(tempDiff <= 0);

%Define Richardson Index (using fitting parameters for wind):
if ~isfield(sLand, 'ri')
   sLand.ri = zeros(size(sCryo.solid)); 
end
sLand.ri(:) = 0; 
sLand.ri(indPos) = 2*9.81*tempDiff(indPos)./(tempAvg(indPos)*uStab) ; %Dimensionless
sLand.ri(indNeg) = 2*9.81*tempDiff(indNeg)./(tempAvg(indNeg)*uUnst) ; %Dimensionless

%C_h = 0.4^2/Log(Z_air/Z_0)^2 (see page 161-162)
    %Za = 2 m (measurement height of air temp)
    %Z_0 = 0.0025 m (aerodynamic roughness height for snow)
%chnS = (0.4^2)/log(2/0.0025)^2;
    %For glaciers, use Z_0 = 0.02 m
%chnI = (0.4^2)/log(2/0.02)^2;
chnS = 0.0036; %Dimensionless
chnI = 0.0075;


%Find indices of unstable, stable, and neutral atmospheric conditions:
indUnst = find(sLand.ri  < riUnst               );
indNeut = find(sLand.ri >= riUnst & sLand.ri <= riStab);
indStab = find(                sLand.ri > riStab);

%Define bulk transfer coeffient for neutral conditions:
if ~isfield(sCryo,'chn')
   sCryo.chn = nan(size(sCryo.solid)); 
end
sCryo.chn(:) = nan;
sCryo.chn(sCryo.solid > 0                 ) = chnS;
sCryo.chn(sCryo.solid <= 0 & sCryo.ice > 0) = chnI;

%Initialize heat argument:
if ~isfield(sCryo,'heatSC')
    sCryo.heatSC = zeros(size(sCryo.solid));
end
sCryo.heatSC(:) = 0;

%Populate heat array:
%Sensible heat transfer for three different regimes (see DeWalle and Rango, pg. 161-164):
%Density of air (1.3 kg/m^3) * specific heat capacity (1005 J/kg/Deg C) * bulk transfer coefficient (related to atm stability; unitless) * u_air (NOT INCLUDED; m/s) * (Tair - Tsnow)
sCryo.heatSC(indUnst) = 1306.5*uUnst.*sCryo.chn(indUnst).*((1-16*sLand.ri(indUnst)).^(0.75)).*tempDiff(indUnst);
sCryo.heatSC(indStab) = 1306.5*uStab.*sCryo.chn(indStab).*((1- 5*sLand.ri(indStab)).^2).*tempDiff(indStab);
sCryo.heatSC(indNeut) = 1306.5*mean([uUnst,uStab]).*sCryo.chn(indNeut).*tempDiff(indNeut);

% %Loop over three regimes (stable, neutral, and unstable)
% %Calculate C_hn, which is c_h*(C_h/C_hn)
% for ii = 1 : 3
%     switch ii
%         case 1 %Unstable (Ch/Chn -> 0)
%             ind = indUnst;
%                 %Ch/Chn = (1 - 16*Ri)^0.75
%             scalar = sCryo.chn(ind).*(1-16*ri(ind)).^(0.75);
%         case 2 %Stable (Ch/Chn >> 1)
%             ind = indStab;
%                 %Ch/Chn = (1 - 5*Ri)^2
%             scalar = sCryo.chn(ind).*(1-5*ri(ind)).^2;
%         case 3 %Neutral (Ch/Chn ~ 1)
%             ind = indNeut;
%             scalar = sCryo.chn(ind);
%     end
%     %Density of air (1.3 kg/m^3) * specific heat capacity (1005 J/kg/Deg C) * bulk transfer coefficient (related to atm stability; unitless) * u_air (NOT INCLUDED; m/s) * (Tair - Tsnow)
%     sCryo.heatSC(ind) = 1306.5*scalar.*tempDiff(ind); %Units are W/m^2
% end
