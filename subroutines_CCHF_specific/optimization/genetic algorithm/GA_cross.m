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
% <http://www.gnu.org/licenses/>.s


function coef = GA_cross(coef, sOpt)


nCrossPts = find_att(sOpt.att,'crossover sites');
nMemCross = find_att(sOpt.att,'crossover pop');
nElite = find_att(sOpt.att,'elite pop','no_warning');

if isempty(nElite)
   nElite = 0; 
end


% %Randomly select indices for chromosomes to cross
% xTemp = rand(numel(coef(:,1)), 1);

% %Using indices of the sort function ensures no repreated indices
% [~, indCross] = sort(xTemp);
% %Keep only the percent that are supposed to be crossed:
% indCross = indCross(1:nMemCross);
% %Ensure pairs:
% if mod(numel(indCross),2) == 1    
%     indCross(end) = []; 
% end

% [~, ordX] = sort(rand(numel(coef(:,1)),1));

if mod(nMemCross,2)
    nMemCross = nMemCross + 1;
end
%Loop over parameter sets, crossing each consecutive pair (assumes they've been randomized):
for ii = 1 : nMemCross/2
    indCurr = nElite + [2*(ii-1)+1, 2*(ii-1)+2];
%     memX = ordX(indCurr);
    %Select sites where crossing occurs:
    indsites = sort(randi([2, numel(coef(1,:))-1], nCrossPts, 1));
    %Perform crossover for current pair:
    for jj = 1 : nCrossPts
        coef(indCurr(1), : ) = [coef(indCurr(1),1:indsites(jj)), coef(indCurr(2),indsites(jj)+1:end)];
        coef(indCurr(2), : ) = [coef(indCurr(2),1:indsites(jj)), coef(indCurr(1),indsites(jj)+1:end)];
    end
end