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

function [coef, nRunPr] = uniform_sample(prmBnds, nRun)

nCoef = numel(prmBnds(:,1));

intDelta = round(nRun^(1/nCoef));

%Round down number of divisions in parameter space if runs required is
%significantly greater than number expected by user.
if intDelta^nCoef ~= nRun
    if intDelta^nCoef + 100 > nRun
        intDelta = intDelta - 1;
    end
    nRunPr = intDelta^nCoef;
    warning('uniform_sample_SETI:runsChange',['The total number of '...
        'model runs specified was ' num2str(nRun) ', which is being '...
        'changed to ' num2str(nRunPr) ' in order to implement a '...
        'uniform parameter space sampling.' char(10) ...
        'Adding an additional sampling point would result in ' ...
        num2str((intDelta+1)^nCoef) ' model runs.']);
else
    nRunPr = nRun;
end

%Initialize coef output matrix:
% coefValues = nan(nCoef,intDelta);
coef = nan(numel(prmBnds(:,1)),nRunPr);

%Determine all parameter values for uniform parameter space sampling:
coefValues = repmat(1:intDelta,nCoef,1)/(intDelta+1).* ...
    repmat(prmBnds(:,2) - prmBnds(:,1),1,intDelta) + repmat(prmBnds(:,1),1,intDelta); 

%Find all permutations for each coeficient value:
for ii = 1 : nCoef
    nConsec = intDelta^(nCoef-ii);
    
    for jj = 1 : nRunPr / nConsec
        coef(ii, nConsec*(jj-1) + 1 : nConsec*jj) = coefValues(ii,mod(jj-1,intDelta)+1);
    end
end

coef = coef';