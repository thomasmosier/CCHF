% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function varargout = mass_enbal_v2(varargin)

global sCryo


%To evaluate snowpack energy balance, varargin{1} = sMeta
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'sens_pwr', -3, 2, 0, 'snowpack_enbal','cryo'}); %decimal (0->1)
    varargout{1} = cat(1,varargout{1}, {'cond_pwr', -3, 2, -1, 'snowpack_enbal','cryo'}); %decimal (0->1)
    varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr' , -2, 2, -0.2348, 'melt_cc','cryo'});
    return
else
	sMeta = varargin{1};
    SensPwr = find_att(varargin{1}.coef, 'sens_pwr');
    condPwr = find_att(varargin{1}.coef, 'cond');
    dii = find_att(varargin{1}.coef, 'ice_melt_pwr');
end


%Available snow fields:
%sSnowpack.solid
%sSnowpack.tsn
%sSnowpack.liquid

%Calculate heat required to raise snowpack to 0 Deg C:
%DeltaE = m*c_ice*DeltaT
%dt*power = density_H2O*depth_w.e.*c_ice*DeltaT


%Constants used here:
%c_ice   =   2100 J / kg / Deg C
%L_water = 334000 J / kg
%density_water = 1000 kg / m^3
cSensW = 4.18*10^6; %= c_water*density_water: J / m^3 / Deg C
cSensI = 2.027*10^6; %= c_ice*density_water: J / m^3 / Deg C (specific heat capacity of ice/snow ranges from about 2030-2090 J/kg/K)
cLate  = 3.34*10^8; %= L_water*density_water: J / m^3

% 
% if isequal(sMeta.currTime(1:2),[2002,8])
%     keyboard
% end


%Initialize fields on first iteration:
if ~isfield(sCryo,'solid')
    sCryo.solid = zeros(size(sCryo.solid),'single');
end
if ~isfield(sCryo,'liquid')
    sCryo.snlq = zeros(size(sCryo.solid),'single');
end
if ~isfield(sCryo,'tsn')
    sCryo.tsn = zeros(size(sCryo.solid),'single');
end

%%Reset fields every iteration: ('dSnow' is reset in 'snow_accum')
sCryo.sensibleSnow = zeros(size(sCryo.solid),'single');
sCryo.latentSnow   = zeros(size(sCryo.solid),'single');
sCryo.latentIce    = zeros(size(sCryo.solid),'single');
sCryo.dIce       = zeros(size(sCryo.solid),'single');
sCryo.ilr        = zeros(size(sCryo.solid),'single'); %Ice liquid release


%Only execute this function once at end of this function:
% %%ADJUST SOLID AND LIQUID CONTENT DUE TO SNOWPACK TEMPERATURE
% sSnowpack = snowpack_equilibriate(sSnowpack,cSens,cLate);


%%CALCULATE ENERGY (from heat) AND ADJUST SNOWPACK PROPERTIES
%Energy has units of Joules/m^2
energS = 10^(condPwr)*sCryo.heat*time2sec(1,sMeta.dt,sMeta.currTime);



%Find indices with liquid water and those with snow
indLiq = find(sCryo.snlq > 0);
indSnow = find(sCryo.solid > 0);

%Find indices with negative energy and those with positive (treated
%seperately below)
indNeg = find(energS < 0);
indPos = find(energS > 0);


%%FOR LOCATIONS WITH NEGATIVE ENERGY BAL:
%Freeze water content:
indReFreeze = intersect(indNeg,indLiq);
if ~isempty(indReFreeze)
%     sCryo.latentSnow(indReFreeze) = -10^(phasePwr)*cLate*sCryo.snlq(indReFreeze);
    sCryo.latentSnow(indReFreeze) = -cLate*sCryo.snlq(indReFreeze);

    indMaxRf = intersect(indReFreeze, find(sCryo.latentSnow < energS));
    sCryo.latentSnow(indMaxRf) = energS(indMaxRf);

    dSnowCurr = -sCryo.latentSnow(indReFreeze)/(cLate);

    %Refrozen snow has temperature of 0 deg, therefore need to calculate new
    %temperature that is weighting of old and new:
    %ORDER MATTERS: Must be before updating .solid
    dTempCurr = sCryo.tsn(indReFreeze).*(sCryo.solid(indReFreeze)./(sCryo.solid(indReFreeze)+dSnowCurr) - 1);
        dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
    sCryo.tsn(indReFreeze) = sCryo.tsn(indReFreeze) + dTempCurr; 
    
    sCryo.solid(indReFreeze) = sCryo.solid(indReFreeze) + dSnowCurr;
    sCryo.snlq(indReFreeze) = sCryo.snlq(indReFreeze) - dSnowCurr;
    energS(indReFreeze) = energS(indReFreeze) - sCryo.latentSnow(indReFreeze);
end

%Set liquid that is approximately 0 to zero:
epsLq = 10^(-5);
sCryo.snlq(sCryo.snlq > -eps & sCryo.snlq < eps) = 0;

%Lower snow temperature using energy that didn't go into refreezing snow
indLowTemp = intersect(intersect(indNeg,indSnow), find(sCryo.snlq == 0));
if ~isempty(indLowTemp)
    sCryo.sensibleSnow(indLowTemp) = energS(indLowTemp);
    dTempCurr = sCryo.sensibleSnow(indLowTemp) ./ (cSensI*sCryo.solid(indLowTemp));
    % dTempCurr = 10^(sensPwr)*sCryo.sensibleSnow(indLowTemp) ./ (cSens*sCryo.solid(indLowTemp).*sHydro.area(indLowTemp));
        dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
    % sCryo.dTmp(indLowTemp) = sCryo.dTmp(indLowTemp) + dTempCurr;
    sCryo.tsn(indLowTemp) = sCryo.tsn(indLowTemp) + 10^(SensPwr)*dTempCurr;
    energS(indLowTemp) = energS(indLowTemp) - sCryo.sensibleSnow(indLowTemp);
end



%%FOR LOCATIONS WITH POSITIVE ENERGY BAL:
%Raise temperature of snowpack:
indRaisTemp = intersect(indPos,indSnow);
if ~isempty(indRaisTemp)
    %Sensible heat necessary to raise snowpack temp:
    sCryo.sensibleSnow(indRaisTemp) = -cSensI*sCryo.solid(indRaisTemp).*sCryo.tsn(indRaisTemp);
    % sCryo.sensibleSnow(indRaisTemp) = -cSens*sCryo.solid(indRaisTemp).*sCryo.tsn(indRaisTemp).*sHydro.area(indRaisTemp);
    %Max sensible energy is that available:
    indMaxRaise = intersect(indRaisTemp, find(sCryo.sensibleSnow > energS));
    sCryo.sensibleSnow(indMaxRaise) = energS(indMaxRaise);

    dTempCurr = sCryo.sensibleSnow(indRaisTemp) ./ (cSensI*sCryo.solid(indRaisTemp));
    % dTempCurr = 10^(sensPwr)*sCryo.sensibleSnow(indRaisTemp) ./ (cSens*sCryo.solid(indRaisTemp).*sHydro.area(indRaisTemp));
        dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
    % sCryo.dTmp(indRaisTemp) = sCryo.dTmp(indRaisTemp) + dTempCurr;
    sCryo.tsn(indRaisTemp) = sCryo.tsn(indRaisTemp) + dTempCurr;
    energS(indRaisTemp) = energS(indRaisTemp) - sCryo.sensibleSnow(indRaisTemp);
end

%Make sure temperature not increaased too much:
epsTsn = 10^(-1);
sCryo.tsn(sCryo.tsn > -epsTsn) = 0;

%Melt snow:
indMeltS = intersect(find(sCryo.solid > 0), intersect(indPos, find(sCryo.tsn == 0)));
if ~isempty(indMeltS)
    sCryo.latentSnow(indMeltS) = cLate*sCryo.solid(indMeltS);
    % sCryo.latentSnow(indMelt) = cLate*sCryo.solid(indMelt).*sHydro.area(indMelt);
    indMaxMelt = intersect(indMeltS, find(sCryo.latentSnow > energS));
    sCryo.latentSnow(indMaxMelt) = energS(indMaxMelt);

    dSnowCurr = -sCryo.latentSnow(indMeltS)/(cLate);
    % dSnowCurr = -10^(phasePwr)*sCryo.latentSnow(indMelt)./(cLate*sHydro.area(indMelt));


    sCryo.solid(indMeltS) = sCryo.solid(indMeltS) + dSnowCurr;
    sCryo.snlq(indMeltS) = sCryo.snlq(indMeltS) - dSnowCurr;
    % sCryo.dSnow(indMelt) = sCryo.dSnow(indMelt) + dSnowCurr;
    energS(indMeltS) = energS(indMeltS) - sCryo.latentSnow(indMeltS);
end


%Melt ice:
indIceMlt = intersect(intersect(indPos, find(sCryo.ice ~= 0 & ~isnan(sCryo.ice))), energS > 0);
%Calculate changes using residual energy at locations with ice:
if ~isempty(indIceMlt)
    %Calculate ice melt, which involves different degree-index factor
    sCryo.latentIce(indIceMlt) = 10^(dii)*energS(indIceMlt).*sCryo.heatI(indIceMlt)./sCryo.heat(indIceMlt);

    sCryo.latentIce( isnan(sCryo.latentIce)) = 0;
    sCryo.ilr(indIceMlt) = sCryo.latentIce(indIceMlt)/cLate;
    sCryo.dIce(indIceMlt) = -sCryo.latentIce(indIceMlt)/cLate;
end


% sCryo.dIce(indMeltIce) = -10^(dii)*sCryo.latentIce(indMeltIce) ./ (cLate*sHydro.area(indMeltIce)); 
% sCryo.ilr(indMeltIce)  = -10^(dii)*sCryo.latentIce(indMeltIce) ./ (cLate*sHydro.area(indMeltIce)); 



%%MAKE CORRECTIONS TO VARIOUS PROPERTIES:
%Set temperature to 0 where there is liquid water:
sCryo.tsn(sCryo.snlq > 0) = 0;
%Set temp to zero if solid snow exists and temp positive:
sCryo.tsn(sCryo.tsn > 0 & sCryo.solid > 0) = 0;
%Set snow temperature to nan where there is no snow
sCryo.tsn(sCryo.solid == 0) = nan;
%Set netagive snow values to 0:
sCryo.solid(sCryo.solid < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.snlq(sCryo.snlq < 0 ) = 0;

if any(any(sCryo.tsn < sCryo.tsnmin & sCryo.solid > 0.01))
    if ~regexpbl(sMeta.mode, 'calib')
        warning('snowpack_enbal:negSnowTmp',['Snowpack temperature is less than the minimum, which is ' ...
            num2str(sCryo.tsnmin) ' for at least one cell. It is being set to the minimum.']);
    end
    sCryo.tsn(  sCryo.tsn < sCryo.tsnmin) = sCryo.tsnmin; 
end




%%%%%%%%%%%%%%%%%%%%OLD METHOD%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%SENSIBLE HEAT TRANSFER (i.e. raise or lower snow temp):
% indDTmp = find(sSnowpack.solid > 0);
% sSnowpack.dTmp(indDTmp) = energy(indDTmp) ./ (cSensSnow*1000*sSnowpack.solid(indDTmp)); %Energy / (c*m)
% 
% %Ensure maximum snow temperature is 0 Deg C by changing delta Tmp accordingly:
% indSensGZero = find( sSnowpack.tsn + sSnowpack.dTmp > 0 & sSnowpack.solid > 0);
% if ~isempty(indSensGZero)
%     sSnowpack.dTmp(indSensGZero) = -sSnowpack.tsn(indSensGZero);
% end
% %Calculate new snowpack temperature:
% sSnowpack.tsn = sSnowpack.tsn + sSnowpack.dTmp;
% 
% %Calculate sensible snowpack heat transfer:
% sSnowpack.sensibleSnow(indDTmp) = cSensSnow*1000*sSnowpack.solid(indDTmp).*sSnowpack.dTmp(indDTmp);
% %Subtract sensible heat from energy
% energy(indDTmp) = energy(indDTmp) - sSnowpack.sensibleSnow(indDTmp);
% 
% 
% %%LATENT HEAT TRANSFER (i.e. melt or refreeze snow)
% 
% %Refreeze liquid water at locations where snowpack temp less than 0 C:
% indRf = find(sSnowpack.tsn < 0 & sSnowpack.liquid > 0 & sSnowpack.solid > 0);
%     %Cells where snow is below 0 Deg C, there is liquid water, and solid snowpack:
% if ~isempty(indRf)
%     %Calculate refreeze potential (depth of w.e.):
%     %depth_freeze = (specific heat capacity / latent heat of fusion) * depth_SWE * (0 - Temp)
%     rfPot(indRf) = -(cSensSnow/cLateFuse)*sSnowpack.solid(indRf).*sSnowpack.tsn(indRf); 
%     
%     %Set maximum refreeze potential to be liquid water content:
%     rfPot(rfPot > sSnowpack.liquid) = sSnowpack.liquid(rfPot > sSnowpack.liquid);
% 
%     %Add modified refreeze potential to solid snowpack
%     sSnowpack.solid(indRf) = sSnowpack.solid(indRf) + rfPot(indRf); 
%     %Subtract modified refreeze potential from liquid water in snowpack
%     sSnowpack.liquid(indRf) = sSnowpack.liquid(indRf) - rfPot(indRf); 
%     
%     %Update snowchange:
%     sSnowpack.dSnow(indRf) = sSnowpack.dSnow(indRf) + rfPot(indRf);
% 
%     %Adjust snow temperature based on amount of water refrozen:
%     %DeltaTemp = (Latent / sensible) * (depth_refreeze / depth_SWE)
%     sSnowpack.dTmpRefreeze(indRf) = (cLateFuse/cSensSnow)*(rfPot(indRf)./sSnowpack.solid(indRf));
%     sSnowpack.tsn(indRf) = sSnowpack.tsn(indRf) + sSnowpack.dTmpRefreeze(indRf); 
% end
% 
% %Excess energy melts snow:
% indDSnow = find(sSnowpack.solid > 0 & energy > 0);
% if ~isempty(indDSnow)
%     sSnowpack.dSnow(indDSnow) = - energy(indDSnow) ./ cLateFuse;
%     %Set Maximum melt depth to be snow water equivalent:
%     indLateLZero = find( sSnowpack.solid <= -sSnowpack.dSnow & sSnowpack.dSnow ~= 0);
%     if ~isempty(indLateLZero)
%         sSnowpack.dSnow(indLateLZero) = -sSnowpack.solid(indLateLZero);
%     end
% 
%     %Calculate new snowpack water equivalent and liquid water content:
%     sSnowpack.solid(indDSnow) = sSnowpack.solid(indDSnow) + sSnowpack.dSnow(indDSnow);
%     sSnowpack.liquid(indDSnow) = sSnowpack.liquid(indDSnow) - sSnowpack.dSnow(indDSnow);
%     %Update snowchange:
%     sSnowpack.dSnow(indDSnow) = sSnowpack.dSnow(indDSnow) + sSnowpack.dSnow(indDSnow);
% 
%     %Calculate latent heat:
%     sSnowpack.latentSnow(indDSnow) = cLateFuse*1000*sSnowpack.dSnow(indDSnow);
%     %Subtract sensible heat from energy
%     energy(indDSnow) = energy(indDSnow) - sSnowpack.latentSnow(indDSnow);
% end
% 
% %Excess energy melts glaciers:
%     %Assumes glaciers are ripe (i.e. no sensible heating or other
%     %glacier thermodynamics)
% indIceMelt = find(sSnowpack.ice == 1 & energy > 0);
% if ~isempty(indIceMelt)
%     %Because of difference in shortwave radiation absorption by ice (diff. 
%     %albedo), need to scale energy available for melting ice. 
%     sSnowpack.latentIce(indIceMelt) = (sSnowpack.heatI(indIceMelt) ./ sSnowpack.heat(indIceMelt)) .* energy(indIceMelt);
%     sSnowpack.dIce(indIceMelt) = -sSnowpack.latentIce(indIceMelt) / cLateFuse; %Negative is ice shrinking, positive if growing (units of depth)
%     sSnowpack.liquid(indIceMelt)    = sSnowpack.liquid(indIceMelt) - sSnowpack.dIce(indIceMelt); %Subtract ice change, since negative ice change makes positive contribution to water availability
%     energy(indIceMelt) = 0;
% end
% 
% %If somehow snowpack set to nan, correct this:
% sSnowpack.solid(isnan(sSnowpack.solid)) = 0;
% 
% %Set minimum snowpack temp (hopefully only matters when there are small amounts
% %of snow)
% sSnowpack.tsn(sSnowpack.tsn < minSnowTmp) = minSnowTmp;
% 
% %%WRITE OUTPUT
% argout = sSnowpack;
% 
% 
% 
