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

function varargout = melt_cc(varargin)

    
global sCryo sAtm

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'inertia_pwr',  -4,   1,    -1.5, 'melt_inertia','cryo-routing'}); %Units of C
    varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr'   , -2, 2, -0.2348, 'melt_inertia','cryo'});  %Units of mm/hr/Deg C
        
    return
else
    inertia = find_att(varargin{1}.coef,'inertia_pwr'); 
    dii    = find_att(varargin{1}.coef,'ice_melt_pwr');  %Units of mm/hr/Deg C 
    sMeta  = varargin{1};
end

%This conversion complements a similar one in 'ETI_Pellicciotti' and
%'simple_degree'. Makes units consistent with SETI
cLate = 3.34*10^8; %Density water * Latent heat fusion ice; %Units = J/m^3



%Initialize fields on first iteration:
if ~isfield(sCryo,'liquid')
    sCryo.liquid = zeros(size(sCryo.solid),'single'); 
end
if ~isfield(sCryo,'sncc')
    sCryo.sncc       = zeros(size(sCryo.solid),'single'); %Coldcontent array
end

%Reset fields every iteration:
sCryo.latentSnow = zeros(size(sCryo.solid),'single');
sCryo.latentIce  = zeros(size(sCryo.solid),'single');
sCryo.dIce       = zeros(size(sCryo.solid),'single');




%%Calculate coldceont:
%Properties of Cold Content:
    %When positive, snow temp is negative
    %When negative, snow melts

%     if isequal([2000,8],sMeta.currTime)
%         keyboard
%     end

%Calculate potential melt energy (units of m): 
meltDeptPot = time2sec(1,sMeta.dt,sMeta.currTime)*sCryo.heat/cLate;

%Original cold-content definition:
% sCryo.sncc = sCryo.sncc - dRatio*(2100/(3.34*10^4))*sAtm.prsn*squeeze(sAtm.tas(sAtm.indCurr,:,:)); %See DeWalle and Rano pg. 63-64

indNeg = find(meltDeptPot < 0);
indPos = find(meltDeptPot > 0);
%Calculate sensible heat changes to cold-content:
sCryo.sncc(indNeg) = sCryo.sncc(indNeg) - (10^inertia)*meltDeptPot(indNeg);
sCryo.sncc(indPos) = sCryo.sncc(indPos) - meltDeptPot(indPos);

%Calculate latent heat changes to cold-content:
%Two cases are:
    %(1) if rain falls and cc is pos, freeze rain (add to snowpack) and
    %decrease cc (moves towards 0)
    %(2) if snow falls and cc is neg, melt snow (add to release) and
    %increase cc (moves towards 0)
indRainFrz = find(sAtm.rain > 0 & sCryo.sncc > 0);
if ~isempty(indRainFrz)
    rainFrz = sAtm.rain(indRainFrz);
    
    indMax = find(rainFrz > sCryo.sncc(indRainFrz));
    rainFrz(indMax) = sCryo.sncc(indRainFrz(indMax));
    
%     sCryo.sncc(indRainFrz)  = sCryo.sncc(indRainFrz)  - rainFrz;
    sCryo.solid(indRainFrz) = sCryo.solid(indRainFrz) + rainFrz;
    sAtm.rain(indRainFrz)   = sAtm.rain(indRainFrz)   - rainFrz;
end
indSnowMlt = find(sAtm.prsn > 0 & sCryo.sncc < 0);
if ~isempty(indSnowMlt)
    snowMlt = sAtm.prsn(indSnowMlt);
    
    indMax = find(snowMlt > -sCryo.sncc(indSnowMlt));
    snowMlt(indMax) = -sCryo.sncc(indSnowMlt(indMax));
    
%     sCryo.sncc(indSnowMlt)   = sCryo.sncc(indSnowMlt)   + snowMlt;
    sCryo.solid(indSnowMlt)  = sCryo.solid(indSnowMlt)  - snowMlt;
    sCryo.liquid(indSnowMlt) = sCryo.liquid(indSnowMlt) - snowMlt;
    sAtm.prsn(indSnowMlt)    = sAtm.prsn(indSnowMlt)    - snowMlt;
end

%Reset cold content when no snow
sCryo.sncc(sCryo.solid == 0) = 0;

% %Set melt potential to zero at cells where temperature less than threshold:
meltDeptPot(sCryo.sncc > 0) = 0;


%%Remove melt from SWE and add to release 
%(doesn't account for snowpack liquid holding capacity):

%Melt potential all goes into melting snow but is limited by amount of snow
%available
sCryo.latentSnow = meltDeptPot;
indMaxSnow = find(sCryo.latentSnow > sCryo.solid);
sCryo.latentSnow(indMaxSnow) = sCryo.solid(indMaxSnow);
meltDeptPot = meltDeptPot - sCryo.latentSnow;

%Add melted snow to 'release' field and remove from 'solid':
sCryo.liquid = sCryo.liquid + sCryo.latentSnow;
sCryo.solid = sCryo.solid - sCryo.latentSnow;
%Update snowchange field:
sCryo.dSnow = sCryo.dSnow - sCryo.latentSnow;
% sCryo.sncc  = sCryo.sncc  - sCryo.latentSnow;

%%MELT ICE:
%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'ice')
    indIce = find( sCryo.ice ~= 0 & ~isnan(sCryo.ice));
else
    indIce = [];
end

%Calculate changes using residual energy at locations with ice:
if ~isempty(indIce) && any(meltDeptPot(indIce) > 0)
    %Calculate ice melt, which involves different degree-index factor
    sCryo.latentIce(indIce) = 10^(dii)*meltDeptPot(indIce).*sCryo.heatI(indIce)./sCryo.heat(indIce);
    sCryo.latentIce( isnan(sCryo.latentIce)) = 0;
    sCryo.liquid(indIce) = sCryo.liquid(indIce) + sCryo.latentIce(indIce);
    sCryo.dIce(indIce) = -sCryo.latentIce(indIce);
    sCryo.ice(indIce) = sCryo.ice(indIce) - sCryo.latentIce(indIce);
%     sCryo.sncc(indIce)  = sCryo.sncc(indIce)  - sCryo.latentIce(indIce);
end