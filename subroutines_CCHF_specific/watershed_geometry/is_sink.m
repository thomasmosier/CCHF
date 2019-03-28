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

function [sinkBl, dem] = is_sink(dem)


dimDem = size(dem);
sinkBl = nan(dimDem); 

[~, indDem] = sort(dem(:),'ascend');
%Remove all nan elements:
indDem(isnan(dem(indDem))) = [];

for ii = 1 : numel(indDem) %Loop over columns
    [r,c] = ind2sub(dimDem,indDem(ii));
    
    if r == 1 || r == dimDem(1) || c == 1 || c == dimDem(2)
       continue
    end
    
    demCurr = reshape(dem(r-1:r+1, c-1:c+1),[],1);
    [val,~] = min(demCurr);

    ind = find(demCurr == val);
    if numel(ind) == 1 && ind == 5 %If center is lowest, it may be a sink
        sinkBl(r,c) = 1;
        demCurr(5) = [];
        dem(r,c) = min(demCurr);
    end
end