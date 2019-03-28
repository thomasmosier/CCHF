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

function coef = coef_randomize(coef, score, sOpt)
%Puts elite performing coefficient sets first and then randomized the rest:

nElite = find_att(sOpt.att,'elite pop','no_warning');

if isempty(nElite)
   nElite = 0; 
end

[~, coefOrd] = sort(rand( numel(coef(:,1)), 1));

for ii = 1 : nElite
    coefOrd(coefOrd == ii) = [];
end

if nElite > 0
    [~, sortOrd] = sort(score);
    coef = coef(sortOrd,:);
    
    if nElite < numel(coef(:,1))
        coef = [coef(1:nElite,:) ; coef(coefOrd,:)];
    end
else
    coef = coef(coefOrd,:);
end