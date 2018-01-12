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

function snow_tsn2tas
%If no snow at cell, set snow temperature to air temperature 
%(may matter in energy balance):

global sCryo


if ~isfield(sCryo,'tsn')
    sCryo.tsn = zeros(size(sCryo.solid),'single');
end
% if ~isfield(sCryo,'dTmp')
%     sCryo.dTmp = zeros(size(sCryo.solid),'single');
% end

% ind2NoSnow = find(sCryo.solid == 0 & sCryo.liquid == 0);
% [snowRow, snowCol] = ind2sub(size(sCryo.solid), ind2NoSnow);
% szTmp = size(sAtm.tas);
% ind3NoSnow = sAtm.indCurr + (snowRow-1)*szTmp(1) + (snowCol-1)*szTmp(2)*szTmp(1);
% sCryo.tsn(ind2NoSnow) = sAtm.tas(ind3NoSnow);

% sCryo.dTmp(ind2NoSnow) = 0;

sCryo.tsn(sCryo.solid == 0 & sCryo.liquid == 0) = nan;

%If snow temp = nan and there is snow, set to 0
sCryo.tsn(isnan(sCryo.tsn) & (sCryo.solid ~= 0 | sCryo.liquid ~= 0)) = 0;