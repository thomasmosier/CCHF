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

function [area, fdrSteps, totFdr] = area_upslope(sHydro)

nCells = numel(sHydro.lat)*numel(sHydro.lon);
szGrid = [numel(sHydro.lat), numel(sHydro.lon)];
totFdr = zeros(nCells, nCells,'single');
% totFdr = zeros(size(sHydro.fdr));

nGrid = numel(sHydro.fdr(1,:));
indRun = (1:nGrid)';

fdrSteps = zeros(szGrid, 'single');

%First index is that of cell upstream of runoff routing from past
    %time-step:
ii = 0;  


flagFl = 1;
while flagFl == 1 && ii < nGrid %Create a fdr grid for the current set of flow-from -> flow-to points
    currFdr = full(sHydro.fdr^ii);   

    %Find indices that have been routed fewer steps than others:
    indCurrFdr = setdiff(indRun, find(sum(currFdr) == 0));
    
    %Route Water (do here so that all water routed at least one step):    
    totFdr = totFdr + currFdr;
    fdrSteps(indCurrFdr) = fdrSteps(indCurrFdr) + 1;
    
    ii = ii + 1;

    %Exit when no more water to route or tiem has been exceeded at all
        %points:
    if all(all(currFdr == 0)) || ii > nGrid
       flagFl = 0; 
    end
end

area = reshape(totFdr'*reshape(sHydro.area,[],1), szGrid);