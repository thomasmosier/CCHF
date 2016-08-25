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