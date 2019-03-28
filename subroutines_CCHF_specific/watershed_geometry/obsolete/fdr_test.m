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

function fdrO = fdr_test(dem)
%FROM ARCGIS WEBSITE:
%If a cell is lower than its eight neighbors, that cell is given the value 
%of its lowest neighbor, and flow is defined toward this cell. If multiple 
%neighbors have the lowest value, the cell is still given this value, but 
%flow is defined with one of the two methods explained below. This is used 
%to filter out one-cell sinks, which are considered noise.

%If a cell has the same change in z-value in multiple directions and that 
%cell is part of a sink, the flow direction is referred to as undefined. 
%In such cases, the value for that cell in the output flow direction raster 
%will be the sum of those directions. For example, if the change in z-value 
%is the same both to the right (flow direction = 1) and down (flow 
%direction = 4), the flow direction for that cell is 1 + 4 = 5. Cells with 
%undefined flow direction can be flagged as sinks using the Sink function.

%If a cell has the same change in z-value in multiple directions and is not 
%part of a sink, the flow direction is assigned with a lookup table 
%defining the most likely direction. See Greenlee (1987).


%Initalize:
fdrO = nan(size(dem));

sinkVal = 10^3;
for ii = 2 : length(dem(1,:))-1 %Loop over columns
    for jj = 2 : length(dem(:,1))-1
        %Find steepest descent
        [val,ind] = min(reshape(dem(jj-1:jj+1,ii-1:ii+1),1,[]));
        
        if ind == 5 %Local sink (flag)
            fdrO(jj,ii) = sinkVal;
        else
            indMin = find(dem(jj-1:jj+1,ii-1:ii+1) == val);
            
            if numel(indMin) == 1
                fdrO(jj,ii) = fdr_val(indMin);
            else
                fdrO(jj,ii) = sinkVal;
            end
        end
    end
end

%Deal with sinks and multiple descent paths:
indSpec = find(fdrO == sinkVal);
while numel(indSpec) > 0
    disp(num2str(numel(indSpec)))
    
    for ii = 1 : numel(indSpec)
        [r,c] = ind2sub(size(fdrO), indSpec(ii)); 
        fdrO(r,c) = mode2d(fdrO(r-1:r+1,c-1:c+1));
    end
indSpec = find(fdrO == sinkVal);
end
