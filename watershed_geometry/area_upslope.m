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