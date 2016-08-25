function varargout = mass_enbal_noconduct(varargin)

global sCryo


%To evaluate snowpack energy balance, varargin{1} = sMeta
if isempty(varargin(:))
	varargout{1} = cell(2,6);
    varargout{1}(1,:) = {'fuse_pwr', -3, 2, -1, 'snowpack_enbal','cryo'}; %decimal (0->1)
    varargout{1}(2,:) = {'sens_pwr', -3, 2, -1.5, 'snowpack_enbal','cryo'}; %decimal (0->1)
    return
else
	sMeta = varargin{1};
    fusePwr = find_att(varargin{1}.coef, 'fuse_pwr');
    sensPwr = find_att(varargin{1}.coef, 'sens_pwr');
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
cSens   = (10^sensPwr)*2.027*10^6; %= c_ice*density_water: J / m^3 / Deg C
cLate   = (10^fusePwr)*3.34*10^8; %= L_water*density_water: J / m^3
% sCryo.tsnmin = -50; %Minimum temperature snow can be (Deg C)



%Initialize fields on first iteration:
if ~isfield(sCryo,'solid')
    sCryo.solid = zeros(size(sCryo.solid),'single');
end
if ~isfield(sCryo,'liquid')
    sCryo.liquid = zeros(size(sCryo.solid),'single');
end
if ~isfield(sCryo,'tsn')
    sCryo.tsn = zeros(size(sCryo.solid),'single');
end

%%Reset fields every iteration: ('dSnow' is reset in 'snow_accum')
sCryo.dTmp         = zeros(size(sCryo.solid),'single');
sCryo.sensibleSnow = zeros(size(sCryo.solid),'single');
sCryo.latentSnow   = zeros(size(sCryo.solid),'single');
sCryo.latentIce    = zeros(size(sCryo.solid),'single');
sCryo.dIce         = zeros(size(sCryo.solid),'single');

%Only execute this function once at end of this function:
% %%ADJUST SOLID AND LIQUID CONTENT DUE TO SNOWPACK TEMPERATURE
% sSnowpack = snowpack_equilibriate(sSnowpack,cSens,cLate);



%%CALCULATE ENERGY (from heat) AND ADJUST SNOWPACK PROPERTIES
%Energy has units of Joules/m^2
energyS = sCryo.heat*time2sec(1,sMeta.dt,sMeta.currTime);



%Find indices with liquid water and those with snow
indLiq = find(sCryo.liquid > 0);
indSnow = find(sCryo.solid > 0);

%Find indices with negative energy and those with positive (treated
%seperately below)
indNeg = find(energyS < 0);
indPos = find(energyS > 0);

%%FOR LOCATIONS WITH NEGATIVE ENERGY BAL:
%Freeze water content:
indReFreeze = intersect(indNeg,indLiq);

sCryo.latentSnow(indReFreeze) = -cLate*sCryo.liquid(indReFreeze);

indMaxRf = intersect(indReFreeze, find(sCryo.latentSnow < energyS));
sCryo.latentSnow(indMaxRf) = energyS(indMaxRf);

dSnowCurr = -sCryo.latentSnow(indReFreeze) / cLate;

sCryo.solid(indReFreeze) = sCryo.solid(indReFreeze) + dSnowCurr;
sCryo.liquid(indReFreeze) = sCryo.liquid(indReFreeze) - dSnowCurr;
sCryo.dSnow(indReFreeze) = sCryo.dSnow(indReFreeze) + dSnowCurr;
energyS(indReFreeze) = energyS(indReFreeze) - sCryo.latentSnow(indReFreeze);
%Refrozen snow has temperature of 0 deg, therefore need to calculature new
%temperature that is weighting of old and new:
dTempCurr = sCryo.tsn(indReFreeze).*(dSnowCurr./sCryo.solid(indReFreeze) - 1);
    dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
sCryo.dTmp(indReFreeze) = sCryo.dTmp(indReFreeze) + dTempCurr;
sCryo.tsn(indReFreeze) = sCryo.tsn(indReFreeze) + dTempCurr; 


%Lower snow temperature using energy that didn't go into refreezing snow
indLowTemp = intersect(indNeg,indSnow);
sCryo.sensibleSnow(indLowTemp) = energyS(indLowTemp);
dTempCurr = sCryo.sensibleSnow(indLowTemp) ./ (cSens*sCryo.solid(indLowTemp));
    dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
sCryo.dTmp(indLowTemp) = sCryo.dTmp(indLowTemp) + dTempCurr;
sCryo.tsn(indLowTemp) = sCryo.tsn(indLowTemp) + dTempCurr;
energyS(indLowTemp) = energyS(indLowTemp) - sCryo.sensibleSnow(indLowTemp);


%%FOR LOCATIONS WITH POSITIVE ENERGY BAL:
%Raise temperature of snowpack:
indRaisTemp = intersect(indPos,indSnow);


%Sensible heat necessary to raise snowpack temp:
sCryo.sensibleSnow(indRaisTemp) = -cSens*sCryo.solid(indRaisTemp).*sCryo.tsn(indRaisTemp);
%Max sensible energy is that available:
indMaxRaise = intersect(indRaisTemp, find(sCryo.sensibleSnow > energyS));
sCryo.sensibleSnow(indMaxRaise) = energyS(indMaxRaise);

dTempCurr = sCryo.sensibleSnow(indRaisTemp) ./ (cSens*sCryo.solid(indRaisTemp));
    dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
sCryo.dTmp(indRaisTemp) = sCryo.dTmp(indRaisTemp) + dTempCurr;
sCryo.tsn(indRaisTemp) = sCryo.tsn(indRaisTemp) + dTempCurr;
energyS(indRaisTemp) = energyS(indRaisTemp) - sCryo.sensibleSnow(indRaisTemp);

%Melt snow:
indMelt = intersect(indPos, find(sCryo.tsn == 0));
sCryo.latentSnow(indMelt) = cLate*sCryo.solid(indMelt); %Joules is would take to melt all snow in cell
indMaxMelt = intersect(indMelt, find(sCryo.latentSnow > energyS));
sCryo.latentSnow(indMaxMelt) = energyS(indMaxMelt);

dSnowCurr = -sCryo.latentSnow(indMelt) ./ (cLate);

sCryo.solid(indMelt) = sCryo.solid(indMelt) + dSnowCurr;
sCryo.liquid(indMelt) = sCryo.liquid(indMelt) - dSnowCurr;
sCryo.dSnow(indMelt) = sCryo.dSnow(indMelt) + dSnowCurr;
energyS(indMelt) = energyS(indMelt) - sCryo.latentSnow(indMelt);

%Melt ice:
indMeltIce = intersect(intersect(indPos, find(sCryo.ice > 0)), energyS > 0);
sCryo.latentIce(indMeltIce) = energyS(indMeltIce) .* sCryo.heatI(indMeltIce) ./ sCryo.heat(indMeltIce);
sCryo.dIce(indMeltIce) = -sCryo.latentIce(indMeltIce) / cLate; 
sCryo.ice(indMeltIce) = sCryo.ice(indMeltIce) - sCryo.dIce(indMeltIce);

%%ADJUST SOLID AND LIQUID CONTENT DUE TO SNOWPACK TEMPERATURE
sCryo = snowpack_equilibriate(sCryo,cSens,cLate);


%If snow temp less than minimum allowed, reset:
% sCryo.tsn(sCryo.tsn < sCryo.tsnmin) = sCryo.tsnmin;
% sCryo.dTmp(sCryo.dTmp < sCryo.tsnmin) = sCryo.tsnmin;
sCryo.dTmp(sCryo.dTmp > abs(sCryo.tsnmin)) = abs(sCryo.tsnmin);

%If liquid slightly negative, set to 0:
sCryo.liquid(-0.01 < sCryo.liquid & sCryo.liquid < 0) = 0;


%Ensure no snowpack temperature values greater than zero
if any(any(sCryo.tsn > 0.001 & sCryo.solid > 0.01))
    warning('snowpack_enbal:posSnowTmp',['Snowpack temperature is '...
        'positive (mean pos. = ' num2str(mean2d(sCryo.tsn(sCryo.tsn > 0 & sCryo.solid > 0))) ...
        '), which should not happen.']);
end

if any(any(sCryo.tsn < sCryo.tsnmin & sCryo.solid > 0.01))
    if ~regexpbl(sMeta.mode, 'calib')
        warning('snowpack_enbal:negSnowTmp',['Snowpack temperature is less than the minimum, which is ' ...
            num2str(sCryo.tsnmin) ' for at least one cell. It is being set to the minimum.']);
    end
    sCryo.tsn(  sCryo.tsn < sCryo.tsnmin) = sCryo.tsnmin; 
    sCryo.dTmp(sCryo.dTmp < sCryo.tsnmin) = sCryo.tsnmin;
end

%Ensure no negative snowpack values:
if any(any(sCryo.solid < -0.0001))
    warning('snowpack_enbal:negSnow','Snowpack is negative, which should not happen.')
end

%Ensure liquid water not negative
if any(any(sCryo.liquid < 0))
    warning('snowpack_enbal:negLiq','Liquid water content of the snowpack is negative, which should not happen.')
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
% sSnowpack.tsn(sSnowpack.tsn < sCryo.tsnmin) = sCryo.tsnmin;
% 
% %%WRITE OUTPUT
% argout = sSnowpack;
% 
% 
% 
