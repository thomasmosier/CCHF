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

function snow_accum

global sCryo sAtm sLand

if ~isfield(sAtm,'prsn')
    warning('snow_accum:missingPRSN','prsn field missing from sAtm.')
    return
end

if ~isfield(sCryo,'snw')
    sCryo.snw = zeros(size(squeeze(sAtm.pr(1,:,:))),'single');
end
if ~isfield(sCryo,'lwsnl')
    sCryo.lwsnl = zeros(size(squeeze(sAtm.pr(1,:,:))),'single');
end

sLand.rnrf = zeros(size(squeeze(sAtm.pr(1,:,:))),'single');

%Add new snowfall to snowpack:
sCryo.snw = sCryo.snw + sAtm.prsn;
sCryo.sndwe = sCryo.sndwe + sAtm.prsn;

%Where there's snow, add rain to snow liquid, else add to land
indSnow = find(sCryo.snw > 0);
if~isempty(indSnow)
    sCryo.lwsnl(indSnow)  = sCryo.lwsnl(indSnow) + sAtm.rain(indSnow);
    sCryo.sndwe(indSnow) = sCryo.sndwe(indSnow) + sAtm.rain(indSnow);
end
indNoSnow = find(sCryo.snw <= 0);
if~isempty(indNoSnow)
    sLand.rnrf(indNoSnow) =  sAtm.rain(indNoSnow);
end

%%Set internal average snow temp to nan if no snow present
if isfield(sCryo,'tsn')
    sCryo.tsn(sCryo.snw == 0 & sCryo.lwsnl == 0 & sCryo.icx == 0) = nan;
    
    %If snow temp = nan and there is snow, set to atm temp
    indNew = find(isnan(sCryo.tsn) & (sCryo.snw ~= 0 | sCryo.lwsnl ~= 0 | sCryo.icx ~= 0));
    sCryo.tsn(indNew) = 0;
    indNew3 = ind2_3d(size(sAtm.tas),indNew,sAtm.indtas);
    indTas3 = find(sAtm.tas(indNew3) < 0);
    sCryo.tsn(indNew(indTas3)) = sAtm.tas(indNew3(indTas3));
end

% %Estimate density of new snow based on air temperature and formulation in
% %Pomeroy, J. W., & Brun, E. (2001). Physical properties of snow. Snow Ecology, 45–126.
% sCryo.densityF = nan(size(sCryo.snw));
% indSnow = find(sAtm.prsn > 0);
% 
% sCryo.densityF(indSnow) = 67.9 + 51.25*exp(sAtm.tas(ind2_3d(size(sAtm.tas),indSnow,sAtm.indtas))/2.59);


