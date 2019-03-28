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

function varargout = tLag_Liston(sHydro, varargin)

global sCryo sLand
%%Calculate velocity through each cell (from Liston's 'timecoefs.f')

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    varargout{1} = cat(1, varargout{1}, {'time_cf_slow', 0.2, 5, 2.778, 'tLag_Liston','routing'}); %Unitless scalar
    varargout{1} = cat(1, varargout{1}, {'time_cf_fast', 0.5, 5, 3.838, 'tLag_Liston','routing'}); %Unitless scalar
    
    return
else
    sMeta = varargin{1};
    alphaS = find_att(sMeta.coef,'time_cf_slow'); %Unitless scalar
    alphaF = find_att(sMeta.coef,'time_cf_fast'); %Unitless scalar
end



vSlowScIce = 0.12;
vSlowSfIce = 0.20;
vSlowScLand = 0.10;
vSlowSfLand = 0.08;


tFast = sparse(sHydro.fdr);
    tFast(sHydro.dl ~= 0) = nan;
tSlow = sparse(sHydro.fdr);
    tSlow(sHydro.dl ~= 0) = nan;



%Slope as degree:
if ~isfield(sHydro,'slopeLis') || ~isfield(sHydro,'sHydro.slopeCorr')
    sHydro.slopeLis = atand(sHydro.slope); %Slope is in rise/run, but Liston uses degrees
    sHydro.slopeLis(sHydro.slope < 5) = 5; 
    sHydro.slopeLis(sHydro.slope > 45) = 45;
    sHydro.slopeLis = sHydro.slopeLis.^0.8;
    sHydro.slopeLis(sHydro.dl == 0) = 0;
    %
    sHydro.slopeCorr = sHydro.slopeLis/(15^0.8);
end


%Calculate velocities at each point (m/s):
indDl = find(sHydro.dl ~= 0);
[rDl, ~] = ind2sub(size(sHydro.dl), indDl);
for ii = 1 : numel(indDl)
    [rC, cC] = ind2sub(size(sHydro.dem), rDl(ii));
    if sCryo.snw(rC,cC) > 0 && sCryo.icx(rC,cC) > 0 %Snow covered ice
        vFast = max(vSlowScIce, vSlowScIce*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowScIce;
    elseif sCryo.snw(rC,cC) == 0 && sCryo.icx(rC,cC) > 0 %Snow free ice
        vFast = max(vSlowSfIce, vSlowSfIce*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowSfIce;
    elseif sCryo.snw(rC,cC) > 0 && sCryo.icx(rC,cC) == 0 %Snow covered land
        vFast = max(vSlowScLand, vSlowScLand*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowScLand;
    else %Snow free land
        vFast = max(vSlowSfLand, vSlowSfLand*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowSfLand;
    end
    
    %Populate residence time, 'tLag', array:
    tFast(indDl(ii)) = sHydro.dl(indDl(ii)) ./ vFast; %units of seconds
    tSlow(indDl(ii)) = sHydro.dl(indDl(ii)) ./ vSlow; %units of seconds
end
    
%Set cap for slow flow at one month: 
dtMonth = 24*3600*eomday(sMeta.datecurr(1), sMeta.dateCurr(2));
tSlow(tSlow > dtMonth) = dtMonth - 1;



%Adjust slow and fast flow times with fitting parameters: 
sLand.tlags = alphaS*tSlow; 
sLand.tLag = alphaF*tFast;