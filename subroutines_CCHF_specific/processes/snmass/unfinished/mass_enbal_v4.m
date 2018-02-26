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

function varargout = mass_enbal_v4(varargin)

global sCryo sAtm


%To evaluate snowpack energy balance, varargin{1} = sMeta
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'surf_cond_tsn', 0, 1, 0.8, 'snowpack_enbal','cryo'}); %decimal (0->1)
%     varargout{1} = cat(1,varargout{1}, {'late_pwr', -2, 2, -1, 'snowpack_enbal','cryo'}); %decimal (0->1)
    varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr' , -2, 2, 0.11, 'melt_cc','cryo'});
    return
else
	sMeta = varargin{1};
    dii = find_att(varargin{1}.coef, 'ice_melt_pwr');
    condSur = find_att(varargin{1}.coef, 'surf_cond_tsn');
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
densW = find_att(sMeta.global,'density_water');

cSensS   = find_att(sMeta.global,'heat_cap_snow')*densW;%= c_ice*density_water: J / m^3 / Deg C
cLate   = find_att(sMeta.global,'latent_water')*densW; %= L_water*density_water: J / m^3

sdl = cSensS/cLate;
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
sCryo.lhme = sCryo.heat*time2sec(1,sMeta.dt,sMeta.currTime);


%Update snow temperature (then equilibriate snowpack based on positive or 
%negative temperature):
%Assume constant conduction factor for dry snow and full energy absorption
%for wet snow
indLiq = find(sCryo.snlq > 0 & sCryo.solid > 0); 
indNoLiq = intersect(find(sCryo.solid > 0), setdiff((1:numel(sCryo.solid)), indLiq));
sCryo.tsn(indLiq) = sCryo.tsn(indLiq) + sCryo.lhme(indLiq)./(cSensS.*sCryo.solid(indLiq));
sCryo.tsn(indNoLiq) = sCryo.tsn(indNoLiq) + condSur*sCryo.lhme(indNoLiq)./(cSensS.*sCryo.solid(indNoLiq));





%Find liquid water that should be refrozen (equilibriate snowpack
%temperature and liquid content)
indLiqFrz = find(sCryo.tsn < 0 & sCryo.snlq > 0 & sCryo.solid > 0); %Indices where there is positive cold content and there is liquid water in snow
if ~isempty(indLiqFrz)
    %Indices where there's more liquid water than energy to freeze
    indFrzMax = find(sCryo.snlq(indLiqFrz) > sdl*sCryo.tsn.*sCryo.solid );
    if ~isempty(indFrzMax)
        %Order matters:
        sCryo.snlq(indLiqFrz(indFrzMax))  = sCryo.snlq(indLiqFrz(indFrzMax))  - sdl*sCryo.tsn(indLiqFrz(indFrzMax)).*sCryo.solid(indLiqFrz(indFrzMax));
        sCryo.solid(indLiqFrz(indFrzMax)) = sCryo.solid(indLiqFrz(indFrzMax)) + sdl*sCryo.tsn(indLiqFrz(indFrzMax)).*sCryo.solid(indLiqFrz(indFrzMax));
        sCryo.tsn(indLiqFrz(indFrzMax)) = 0;
    end
    
    %Indices where there's sufficient energy to freeze all liquid:
    if isempty(indFrzMax)
        indFrzNMax = indLiqFrz;
    else
        indFrzNMax = setdiff(indLiqFrz, indLiqFrz(indFrzMax));
    end
    if ~isempty(indFrzNMax)
        %*Order of below lines matters*:
        sCryo.tsn(indFrzNMax) = sCryo.tsn(indFrzNMax) - (sCryo.snlq(indFrzNMax)./sCryo.solid(indFrzNMax))/sdl;
        sCryo.solid(indFrzNMax) = sCryo.solid(indFrzNMax) + sCryo.snlq(indFrzNMax);
        sCryo.snlq(indFrzNMax) = 0;
    end
end
 

% 
% %%FREEZE RAIN AND ADD TO ENERGY:
% indRainFrz = find(sAtm.rain > 0 & sCryo.tsn < 0 & sCryo.solid > 0);
% if ~isempty(indRainFrz)
%     frzPot = -sCryo.solid(indRainFrz).*sCryo.tsn(indRainFrz)*cSensS/cLate;
%     indMazFrz = find(frzPot < sAtm.rain(indRainFrz));
% end



%%MELT SNOW
%Melt potential all goes into melting snow but is limited by amount of snow
%available
indMelt = find( sCryo.tsn > 0 & sCryo.solid > 0);
if ~isempty(indMelt)
    sCryo.latentSnow(indMelt) = sdl*sCryo.tsn(indMelt).*sCryo.solid(indMelt);
        indMaxSnow = find(sCryo.latentSnow(indMelt) > sCryo.solid(indMelt));
        sCryo.latentSnow(indMelt(indMaxSnow)) = sCryo.solid(indMelt(indMaxSnow));
    sCryo.tsn(indMelt) = sCryo.tsn(indMelt) - (sCryo.latentSnow(indMelt)./sCryo.solid(indMelt))/sdl;
    
    %Add melted snow to 'release' field and remove from 'solid':
    sCryo.snlq(indMelt) = sCryo.snlq(indMelt) + sCryo.latentSnow(indMelt);
    sCryo.solid(indMelt) = sCryo.solid(indMelt) - sCryo.latentSnow(indMelt);
    sCryo.lhme(indMelt) = sCryo.lhme(indMelt) - cLate*sCryo.latentSnow(indMelt);
end


%%MELT ICE:
%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'ice')
    indIceMlt = find( sCryo.ice ~= 0 & ~isnan(sCryo.ice) & sCryo.lhme > 0);
else
    indIceMlt = [];
end

%Calculate changes using residual energy at locations with ice:
if ~isempty(indIceMlt)
    %Calculate ice melt, which involves different degree-index factor
    sCryo.latentIce(indIceMlt) = 10^(dii)*sCryo.lhme(indIceMlt).*sCryo.heatI(indIceMlt)./sCryo.heat(indIceMlt);
        sCryo.lhme(indIceMlt) = 0;
        
    %Scale melt for debris covered glacier cells:
    if isfield(sCryo,'debris')
        indDebrisMlt = intersect(indIceMlt, find(~isnan(sCryo.debris))); 
    	sCryo.latentIce(indDebrisMlt) = debScl*sCryo.latentIce(indDebrisMlt);
    end
    sCryo.latentIce( isnan(sCryo.latentIce)) = 0;
    sCryo.ilr(indIceMlt) = sCryo.latentIce(indIceMlt);
    sCryo.dIce(indIceMlt) = -sCryo.latentIce(indIceMlt);
end




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
%     keyboard
%     if ~regexpbl(sMeta.mode, 'calib')
%         warning('snowpack_enbal:negSnowTmp',['Snowpack temperature is less than the minimum, which is ' ...
%             num2str(sCryo.tsnmin) ' for at least one cell. It is being set to the minimum.']);
%     end
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
