function varargout = density_Liston(varargin)

%See albedo discussion in:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm

if isempty(varargin(:))
    varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'albedo_fresh', 0.5, 1, 0.71, 'density_Liston','cryo'});
    
    return
else
    sMeta = varargin{1};
%     aFresh = find_att(varargin{1}.coef,'albedo_fresh'); %Albedo of fresh, deep snow
end


rhoIce = find_att(sMeta.global, 'density_ice');
rhoW = find_att(sMeta.global, 'density_water');
% 
% aFresh = find_att(varargin{1}.global,'albedo_snow_fresh'); %Albedo of fresh snow
% aOld   = find_att(varargin{1}.global,'albedo_snow_old');  %Albedo of old/melting snow
A1 = 0.0013; %m^-1 s^-1
A2 = 0.021; %m^3 kg^-1


if ~isfield(sCryo,'rhosn')
    sCryo.rhosn = zeros(size(sCryo.snw), 'single');
end
if ~isfield(sCryo,'lhsnme')
    sCryo.lhsnme = zeros(size(sCryo.snw), 'single');
end


tsnow = 0.5*(sCryo.tsn - 1);
dt = time2sec(1,sMeta.dt,sMeta.dateCurr);

%DENSITY IMPACT OF MELT WATER:
%This is H_new/H_tot = SWE_new/rho_new / (SWE_new/rho_new +
%SWE_old/rho_old) = 1 / (1 + (rho_new*SWE_old)/(rho_old*SWE_new))
fracNM = 1./(1 + (rhoW*sCryo.snw)./(sCryo.lhsnme.*sCryo.rhosn));
fracNM(isnan(fracNM)) = 0;

%Calculate new density based on fresh melt and snowpack
sCryo.rhosn = fracNM*rhoW + (1 - fracNM).*sCryo.rhosn;


%DENSITY OF NEW SNOW:
%Model density of new snow to be function of wetbulb temperature:
%Assume relative humidity is 8% (this is pretty arbitrary but result not very sensitive):
RH = 10;
%Tdew = squeeze(sAtm.tas) - ((100 - RH)/5);
%Dewpoint depression = tas - tDew = ((100 - RH)/5)
dewPtDep = (100 - RH)/5;
tasWB = squeeze(sAtm.tas) - dewPtDep/3;
%rhoNS = 50 + 1.7*(tasWB + 273.15 - 258.16).^1.5;
rhoNS = 50 + 1.7*(tasWB + 14.99).^1.5;
%Calculate the fractional height of new snow versus old snow:
%This is H_new/H_tot = SWE_new/rho_new / (SWE_new/rho_new +
%SWE_old/rho_old) = 1 / (1 + (rho_new*SWE_old)/(rho_old*SWE_new))
fracN = 1./(1 + (rhoNS.*(sCryo.snw - sAtm.prsn))./(sAtm.prsn.*sCryo.rhosn));
fracN(isnan(fracN)) = 0;


%Calculate new density based on fresh snowfall and old snowpack
sCryo.rhosn = fracN.*rhoNS + (1 - fracN).*sCryo.rhosn;

%DENSIFICATION:
sCryo.rhosn = sCryo.rhosn ...
    + dt*A1*0.5*sCryo.snw.*sCryo.rhosn.*exp(0.08*tsnow).*exp(-A2*sCryo.rhosn);

%Ensure density of snowpack layer less than or equal to density of ice:
sCryo.rhosn(sCryo.rhosn > rhoIce) = rhoIce;

