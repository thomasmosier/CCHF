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

function upslope = upstream(fdr, lon, lat, pt)


%fdr = sparse matrix where row is flowFrom cell and column is the cell that
%flow goes to
%pt = [lat, lon]

%upslope = boolean with all upstream cells marked as 1 and others nan.
szUp = [numel(lat) numel(lon)];
upslope = nan(szUp);

colPt = min(abs(pt(1)-lon));
rowPt = min(abs(pt(2)-lat));

indPt = sub2ind(szUp, rowPt, colPt);

indFlowF = find(fdr(:,indPt) == 1);

upslope(indFlowF) = 1;

iiLoop = 1;
while ~isempty(indFlowF)
    iiLoop = iiLoop + 1;
    
    flowTemp = fdr^ii;
    indFlowF = find(flowTemp(:,indFlowF) == 1);

    upslope(indFlowF) = 1;
end