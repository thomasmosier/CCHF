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
