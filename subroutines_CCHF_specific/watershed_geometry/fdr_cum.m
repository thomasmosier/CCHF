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

function [totFdr, fdrSteps, cumTime] = fdr_cum(sHydro, tLag, fdrSteps, cumTime, tMax)


totFdr = sparse(numel(sHydro.dem),numel(sHydro.dem));
% totFdr = zeros(numel(sHydro.dem),numel(sHydro.dem),'single');


cumTime = double(cumTime); %Must be double to add to a sparse matrix.

nGrid = numel(sHydro.fdr(1,:));
indRun = (1:nGrid)';

%First index is that of cell upstream of runoff routing from past
    %time-step:
ii = double(min2d(fdrSteps) + 1);  

flagFl = 1;
while flagFl == 1 && ii < nGrid %Create a fdr grid for the current set of flow-from -> flow-to points
    currFdr = sHydro.fdr^ii;  
        %In FDR: row is cell of origin and col is downstream cell

    %Find indices that have been routed fewer steps than others:
    indFdrSkip = find(fdrSteps + 1 > ii);
    indCurrFdr = setdiff(indRun, indFdrSkip);
    
    %Route Water (do here so that all water routed at least one step):    
    totFdr = totFdr + currFdr;
    fdrSteps(indCurrFdr) = fdrSteps(indCurrFdr) + 1;
    
    %Add time lag to cumulative time array
    tTravCurr = diag(currFdr*tLag,0);

    cumTime(indCurrFdr) = cumTime(indCurrFdr) + tTravCurr(indCurrFdr); %(units = seconds)

    currFdr(cumTime > tMax,:) = 0;
    currFdr(indFdrSkip,:) = 0;
    ii = ii + 1;

    %Exit when no more water to route or tiem has been exceeded at all
        %points:
    if ~any(any(currFdr ~= 0)) || all(all(cumTime > tMax)) %Logical tests must be in this form because currFDR is much larger than cumTime
       flagFl = 0; 
    end
end


% while flagRf == 1 && ii < nGrid
%     currFdr = full(sHydro.fdr^ii);   
%     
%     %Add time lag to cumulative time array
%     cumTime(indRun) = cumTime(indRun) + diag(currFdr(indRun,:)*sLand.tLag(:,indRun),0); %(units = seconds)
%     fdrSteps(indRun) = ii;
%     
%     indRem = find(cumTime(indRun) >= tMax);
%     indRun(indRem) = [];
%     currFdr(indRem,:) = 0;
%     ii = ii + 1;
%     
%     totFdr = totFdr + currFdr;
%     
%     %Exit when no more water to route or tiem has been exceeded at all
%     %points:
%     if all(all(currFdr == 0)) || isempty(indRun)
%        flagRf = 0; 
%     end
% end