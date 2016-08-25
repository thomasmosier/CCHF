function newChrom = chrom_roulette(oldChrom, error, sOpt)

perFit = find_att(sOpt.att,'fitness requirement');
nElite = find_att(sOpt.att,'elite','no_warning');

if isempty(nElite)
   nElite = 0; 
end

%Implement a piecewise function resembling exponential decay.  This
%determines what percentage of the population are eligible for reproduction
%based on their fitness score.  The idea if there is large deviation
%amongst fitness of members, a smaller portion of the population used;
%conversely, if fitnesses are somewhat close, large percentage of
%population reproduces.
if ischar(perFit) && regexpbl(perFit,'exp_decay')
    keyboard
    if std(error) < 0.5
       perFit = 1; 
    elseif std(error) < 1
        perFit = 0.85;
    elseif std(error) < 10
        perFit = 0.5;
    elseif std(error) < 50
        perFit = 0.1;
    elseif std(error) < 200
         perFit = 0.05;
    else
        perFit = 0.005;
    end
end
    

%Initialize 'newChrom':
newChrom = cell(size(oldChrom));

if ~isempty(nElite)
    nSelect = numel(error) - nElite;
else
    nSelect = numel(error);
end


%%USE ONLY CHROMOSOMES MEETING FITNESS REQUIREMENT:
%Take top n individuals, according to fitness requirement
nFit = ceil(perFit)*numel(error);
%Adjust number of entries selected (min of two and not more than max available):
if nFit == 1 && numel(error) > 1
    nFit = 2;
elseif nFit > numel(error)
    nFit = numel(error);
end

[error, ordError ] = sort(error);
ordError = ordError(1:nFit);
error = error(1:nFit);
oldChrom = oldChrom(ordError);

%%USE ROULETTE SELECTION TO DETERMINE PARENTS FOR NEXT GENERATION
%Calculate weighting factors for roullette selection based on errors.  Use 
%inverse of errors so that members with lowest errors are the most likely
%to be selected.
accumError = cumsum(error.^(-1)) / sum(error.^(-1));

%Calculate random numbers used for selection:
vecRand = rand(nSelect,1);
%Implement fitness proportionate selection
for ii = 1 : nSelect
    indSel = find(vecRand(ii) <= accumError, 1, 'first');
    newChrom{ii} = oldChrom{indSel};
end
 
if ~isempty(nElite) && nElite > 0
    newChrom{nSelect+1:nSelect + nElite} = oldChrom{1:nElite};
end
