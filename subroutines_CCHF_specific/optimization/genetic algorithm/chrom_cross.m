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

function chromosomes = chrom_cross(chromosomes, sOpt)


nCrossPts = find_att(sOpt.att,'crossover sites');
perCross = find_att(sOpt.att,'crossover rate');

%Randomly select indices for chromosomes to cross
xTemp = rand(numel(chromosomes), 1);

%Using indices of the sort function ensures no repreated indices
[~, indCross] = sort(xTemp);
%Keep only the percent that are supposed to be crossed:
indCross = indCross(1:ceil(perCross*numel(chromosomes)));
%Ensure pairs:
if mod(numel(indCross),2) == 1    
    indCross(end) = []; 
end


%Loop over chromosomes, crossing each consecutive pair:
for ii = 1 : length(indCross)/2
    %Select sites where crossing occurs:
    indsites = randi([1, numel(chromosomes{1})], nCrossPts, 1);
    %Perform crossover for current pair:
    for jj = 1 : nCrossPts
        chromosomes{indCross(2*(ii-1)+1), : } = [chromosomes{indCross(2*(ii-1)+1)}(1:indsites(jj)), chromosomes{indCross(2*(ii-1)+2)}(indsites(jj)+1:end)];
        chromosomes{indCross(2*(ii-1)+2), : } = [chromosomes{indCross(2*(ii-1)+2)}(1:indsites(jj)), chromosomes{indCross(2*(ii-1)+1)}(indsites(jj)+1:end)];
    end
end