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

function snowpack_equilibriate(cSens,cLate)

global sCryo
%%ADJUST SOLID AND LIQUID CONTENT DUE TO SNOWPACK TEMPERATURE

%If temperature greater than 0, melt snow and lower temperature:
indMelt = find(sCryo.tsn > 0 & sCryo.solid > 0);
dSnowCurr = -(cSens/cLate)*sCryo.solid(indMelt).*sCryo.tsn(indMelt);
%Set max melt to amount of available snow:
indMaxMelt = find(dSnowCurr < -sCryo.solid(indMelt));
dSnowCurr(indMaxMelt) = -sCryo.solid(indMelt(indMaxMelt));
%Determine temperature effect:
dTempCurr = (cLate/cSens)*dSnowCurr./sCryo.solid(indMelt);
    dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
%Adjust snowpack properties:
sCryo.solid( indMelt) = sCryo.solid( indMelt) + dSnowCurr;
sCryo.liquid(indMelt) = sCryo.liquid(indMelt) - dSnowCurr;
% sCryo.dSnow( indMelt) = sCryo.dSnow( indMelt) + dSnowCurr;
% sCryo.dTmp(indMelt) = sCryo.dTmp(indMelt) + dTempCurr; 
sCryo.tsn( indMelt) = sCryo.tsn( indMelt) + dTempCurr;

%If temperature less than 0, freeze liquid and raise temperature:
indFreez = find(sCryo.tsn < 0 & sCryo.liquid > 0);
dSnowCurr = -(cSens/cLate)*sCryo.liquid(indFreez).*sCryo.tsn(indFreez);
%Set max freeze to amount of available water:
indMxFreez = find(dSnowCurr > sCryo.liquid(indFreez));
dSnowCurr(indMxFreez) = sCryo.liquid(indFreez(indMxFreez));
%Determine temperature effect:
dTempCurr = (cLate/cSens)*dSnowCurr./sCryo.liquid(indFreez);
    dTempCurr(isnan(dTempCurr)) = 0; %Occurs when numerator and denominator are NaN
%Adjust snowpack properties:
sCryo.solid( indFreez) = sCryo.solid( indFreez) + dSnowCurr;
sCryo.liquid(indFreez) = sCryo.liquid(indFreez) - dSnowCurr;
% sCryo.dSnow( indFreez) = sCryo.dSnow( indFreez) + dSnowCurr;
% sCryo.dTmp(indFreez) = sCryo.dTmp(indFreez) + dTempCurr; 
sCryo.tsn( indFreez) = sCryo.tsn( indFreez) + dTempCurr;

