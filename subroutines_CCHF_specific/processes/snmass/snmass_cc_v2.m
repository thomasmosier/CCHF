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

function varargout = snmass_cc_v2(varargin)

    
global sCryo

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'cc_pwr', -2,   1,    -1, 'snmass_cc_v2','cryo'}); %Unitless; impacts accumulation of cold-content (negative energy)
    
    %Snow temperature only impacts heat flux (e.g. longwave out)
    if isfield(sCryo,'tsn')
        varargout{1} = cat(1,varargout{1}, {'tsn_pwr', -2,   2,    -1, 'snmass_cc_v2','cryo'});
    end
    
    return
else
    ccPwr = find_att(varargin{1}.coef,'cc_pwr'); 
    
    if isfield(sCryo,'tsn')
        tsnPwr = find_att(varargin{1}.coef,'tsn_pwr');
    end
    
    sMeta = varargin{1};
end

%This conversion complements a similar one in 'ETI_Pellicciotti' and
%'simple_degree'. Makes units consistent with SETI
densW = find_att(sMeta.global,'density_water'); 
cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3

%Initialize release field:
sCryo.lhsnme = zeros(size(sCryo.snw),'single');

if ~isfield(sCryo, 'snlw')
    sCryo.snlw = zeros(size(sCryo.snw),'single');
end
if ~isfield(sCryo, 'sncc')
    sCryo.sncc = zeros(size(sCryo.snw),'single');
end

if isfield(sCryo,'hfsnc')
   error('mass_cc:conduct', ['The CC conservation of energy and mass '...
       'formulation cannot be used with a heat formulation that '...
       'includes snowpack heat conduction because it assumes the snow '...
       'is isothermal.']); 
end

% if isequal(sMeta.dateCurr, [2003,5]) || isequal(sMeta.dateCurr, [2003,6]) 
%    keyboard 
% end

%Calculate melt potential using simple degree indec formulation (units of m): 
sCryo.lhpme = time2sec(1,sMeta.dt,sMeta.dateCurr)*sCryo.hfnet/cLate;


%Update cold-conent:
%Add thermal inertia to the cold-content parameter:
indNegEn = find(sCryo.lhpme < 0);
if ~isempty(indNegEn)
   sCryo.sncc(indNegEn) = sCryo.sncc(indNegEn) - 10^(ccPwr)*sCryo.lhpme(indNegEn);
end

%Remove cold-content and melt energy:
indPosEn = find(sCryo.lhpme > 0 & sCryo.sncc > 0);
if ~isempty(indPosEn)
    tempDcc = sCryo.sncc(indPosEn) - sCryo.lhpme(indPosEn);
    indMaxCC = find(tempDcc < 0);
    if ~isempty(indMaxCC)
        sCryo.lhpme(indPosEn(indMaxCC)) = -tempDcc(indMaxCC);
        sCryo.sncc(indPosEn(indMaxCC)) = 0;
    end
    indNMaxCC = setdiff((1:numel(indPosEn)),indMaxCC);
    if ~isempty(indNMaxCC)
        sCryo.sncc(indPosEn(indNMaxCC)) = tempDcc(indNMaxCC);
        sCryo.lhpme(indPosEn(indNMaxCC)) = 0;
    end
end


%Reset cold content if no snow:
sCryo.sncc(sCryo.sncc ~= 0 & sCryo.snw == 0) = 0;



%
%Refreeze liquid in snow; if snow contains liquid and there is
%rain, release rain, otherwise hold in snowpack
indLiqFrz = find(sCryo.sncc > 0 & sCryo.snlw > 0 & sCryo.snw > 0); %Indices where there is positive cold content and there is liquid water in snow
if ~isempty(indLiqFrz)
    %Indices where there's more liquid water than energy to freeze
    indFrzMax = find(sCryo.snlw(indLiqFrz) > sCryo.sncc(indLiqFrz) );
    if ~isempty(indFrzMax)
        sCryo.snw(indLiqFrz(indFrzMax)) = sCryo.snw(indLiqFrz(indFrzMax)) + sCryo.sncc(indLiqFrz(indFrzMax));
        sCryo.snlw(indLiqFrz(indFrzMax)) = sCryo.snlw(indLiqFrz(indFrzMax)) - sCryo.sncc(indLiqFrz(indFrzMax));
        sCryo.sncc(indLiqFrz(indFrzMax)) = 0;
    end
    
    %Indices where there's sufficient energy to freeze all liquid:
    if isempty(indFrzMax)
        indFrzNMax = indLiqFrz;
    else
        indFrzNMax = setdiff(indLiqFrz, indLiqFrz(indFrzMax));
    end
    if ~isempty(indFrzNMax)
        sCryo.snw(indFrzNMax) = sCryo.snw(indFrzNMax) + sCryo.snlw(indFrzNMax);
        sCryo.sncc(indFrzNMax) = sCryo.sncc(indFrzNMax) - sCryo.snlw(indFrzNMax);
        sCryo.snlw(indFrzNMax) = 0;
    end
end


% %Freeze Rain if cold content > 0
% indRainFrz = find(sAtm.rain > 0 & sCryo.sncc > 0);
% if ~isempty(indRainFrz)
%     indFrzMax = find(sAtm.rain(indRainFrz) > sCryo.sncc(indRainFrz));
%     if ~isempty(indFrzMax)
%         sCryo.snw(indRainFrz(indFrzMax)) = sCryo.snw(indRainFrz(indFrzMax)) + sCryo.sncc(indRainFrz(indFrzMax));
%         sCryo.rFrz(indRainFrz(indFrzMax)) = sCryo.sncc(indRainFrz(indFrzMax));
%         sAtm.rain(indRainFrz(indFrzMax)) = sAtm.rain(indRainFrz(indFrzMax)) - sCryo.rFrz(indRainFrz(indFrzMax));
%         %Must occur last:
%         sCryo.sncc(indRainFrz(indFrzMax)) = 0;
%     end
%     
%     %Indices where there's sufficient energy to freeze all liquid:
%     if isempty(indFrzMax)
%         indFrzNMax = indRainFrz;
%     else
%         indFrzNMax = setdiff(indRainFrz, indRainFrz(indFrzMax));
%     end
%     if ~isempty(indFrzNMax)
%         sCryo.rFrz(indFrzNMax) = sAtm.rain(indFrzNMax);
%         sCryo.snw(indFrzNMax) = sCryo.snw(indFrzNMax) + sCryo.rFrz(indFrzNMax);
%         sCryo.sncc(indFrzNMax) = sCryo.sncc(indFrzNMax) - sCryo.rFrz(indFrzNMax);
%         sAtm.rain(indFrzNMax) = 0;
%     end
% end


%Melt potential all goes into melting snow but is limited by amount of snow
%available
indMelt = find( sCryo.lhpme > 0 & sCryo.sncc <= 0);
if ~isempty(indMelt)
    sCryo.lhsnme(indMelt) = sCryo.lhpme(indMelt);
    indMaxSnow = find(sCryo.lhsnme(indMelt) > sCryo.snw(indMelt));
    sCryo.lhsnme(indMelt(indMaxSnow)) = sCryo.snw(indMelt(indMaxSnow));
    sCryo.lhpme(indMelt) = sCryo.lhpme(indMelt) - sCryo.lhsnme(indMelt);
    
    %Add melted snow to 'release' field and remove from 'solid':
    sCryo.snlw(indMelt) = sCryo.snlw(indMelt) + sCryo.lhsnme(indMelt);
    sCryo.snw(indMelt) = sCryo.snw(indMelt) - sCryo.lhsnme(indMelt);
end



%Set negative cold content values to 0:
sCryo.sncc(sCryo.sncc < 0) = 0;
%Set netagive solid snow values to 0:
sCryo.snw(sCryo.snw < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.snlw(sCryo.snlw < 0 ) = 0;


%ADJUST SNOW TEMPERATURE BASED ON COLD-CONTENT
if isfield(sCryo,'tsn') 
    minTemp = find_att(sMeta.global,'snow_temp_min'); 
    
    cSensW = densW*find_att(sMeta.global, 'heat_cap_water'); %= c_water*density_water: J / m^3 / Deg C
    cSensI = densW*find_att(sMeta.global, 'heat_cap_ice'); %= c_ice*density_water: J / m^3 / Deg C (specific heat capacity of ice/snow ranges from about 2030-2090 J/kg/K)
%     cSensI = 2.027*10^6; %= c_ice*density_water: J / m^3 / Deg C (specific heat capacity of ice/snow ranges from about 2030-2090 J/kg/K)

    indSnow = find(sCryo.snw ~= 0);
    %From def. of cold-content, T = (cc*rho_w*L)/(rho_s*c*SWE)
    sCryo.tsn(indSnow) = -10^(tsnPwr)...
        *(cLate*sCryo.sncc(indSnow)...
        ./(cSensI*sCryo.snw(indSnow) + cSensW*sCryo.snlw(indSnow)));
        %In the above, there shouldn't actually be any liquid water at
        %locations where cc is positive

    
    %%MAKE CORRECTIONS TO TEMPERATURE BASED ON LIMITS AND ASSUMPTIONS:
    %ADJUST SNOW INTERNAL TEMP AND OTHER PARAMETERS:
    %Set netagive solid snow values to 0:
    sCryo.snw(sCryo.snw < 0 ) = 0;
    %Set negative snow liquid values to 0:
    sCryo.snlw(sCryo.snlw < 0 ) = 0;

    sCryo.tsn(sCryo.tsn > 0) = 0;
    sCryo.tsn(sCryo.tsn < minTemp) = minTemp; 
    sCryo.tsn(sCryo.snlw > 0.005*sCryo.snw) = 0;
    sCryo.tsn(sCryo.snw == 0) = nan; 
    
    %Set surface skin temperature equal to average internal temperature (i.e. isothermal assumption):
    sCryo.tssn = sCryo.tsn;
end
