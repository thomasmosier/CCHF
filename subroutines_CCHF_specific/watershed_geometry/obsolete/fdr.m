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

function fdrO = fdr(dem,varargin)
%RULES FROM http://spatial-analyst.net/ILWIS/htm/ilwisapp/flow_direction_algorithm.htm
%STEEPEST SLOPE METHOD:
%For each block of 3x3 input pixels, the operation calculates the height 
    %difference between the central pixel (CP) and each of the 8 neighbour 
    %pixels.
%If, for a neighbour, the height difference is positive (i.e. central pixel 
    %has larger value than the specific neighbour), then:
        %for corner neighbours, height differences are divided by (distance) 1.4
        %for horizontal neighbours, height differences are divided by (distance) 1
%This determines the steepness between the central pixel and its neigbours.
%Then, the (position of the) neighbour with the largest 'steepness' value 
    %is the output flow direction for the current central pixel.

%Additional rules that are applied:
%Pixels along the edges of the input map (margins and corners) will always return the undefined value in the output map.
%If all neighbouring pixels of a central pixel have a larger value than the central pixel itself (i.e. a sink or pit), the undefined value will be returned in the output map. In principle, this situation should not occur; it is advised to use a sink-free input DEM.
%When a central pixel has the undefined value, undefined will be returned in the output map.
%Neighbour pixels that have the undefined value are ignored during the calculation.
%If from eight neighbour pixels considered, three adjacent neighbour pixels in a single row or column are found to have the same steepest slope or the same smallest height value, then the position of the neighbour pixel that is in the middle of those three neighbours is used.
%If from eight neighbour pixels considered, two neighbour pixels are found to have the same steepest slope or smallest height value, then the position of one of these two neighbour pixels is used arbitrarily.



dimDem = size(dem);
fdrO = nan(dimDem); 

%Determine how to treat frame:
strFrame = 'nan';
if ~isempty(varargin(:))
    if regexpbl(varargin{1},{'in','out'})
        strFrame = varargin{1};
    else
        warning('fdr:borderUnknown','An unknown option was selected for defining flow at the border.  Frame will be left as nan.');
    end
end


%Find sinks and use sink-filled DEM in processing
[~, dem] = is_sink(dem);

%Sort DEM highest to lowest
[~, indDem] = sort(dem(:),'ascend');
%Remove all nan elements:
indDem(isnan(dem(indDem))) = [];


%IndDiag used to adjust corner slopes.
indDiag = [1,3,7,9];

% indVer = [4,6];
% indHor = [2,8];
% jj = round(length(lat(:,1))/2);
% ii = round(length(lon(1,:))/2);
% %Differential length factors:
% facDiag = acosd(sind(lat(jj,ii))*sind(lat(jj-1,ii-1))+cosd(lat(jj,ii))*cosd(lat(jj-1,ii-1)).*cosd(lon(jj,ii)-lon(jj-1,ii-1)));
% facHor = acosd(sind(lat(jj,ii))*sind(lat(jj,ii-1))+cosd(lat(jj,ii))*cosd(lat(jj,ii-1)).*cosd(lon(jj,ii)-lon(jj,ii-1)));
% facVer = acosd(sind(lat(jj,ii))*sind(lat(jj-1,ii))+cosd(lat(jj,ii))*cosd(lat(jj-1,ii)).*cosd(lon(jj,ii)-lon(jj-1,ii)));
% 
% facHor = facHor / facDiag;
% facVer = facVer / facDiag;

%Actual value:
facDiag = 1.414;

for ii = 1 : numel(indDem)
    [r,c] = ind2sub(dimDem,indDem(ii));
    
    if r == 1 || r == dimDem(1) || c == 1 || c == dimDem(2)
       continue
    end
    
    %Find change in Elevation for from current cell to surrounding frame
    deltaHCurr = reshape(dem(r-1:r+1,c-1:c+1),[],1) - dem(r,c);
%     deltaHCurr(indHor) = deltaHCurr(indHor) / facHor;
%     deltaHCurr(indVer) = deltaHCurr(indVer) / facVer;
    deltaHCurr(indDiag) = deltaHCurr(indDiag) / facDiag;

    [val,~] = min(deltaHCurr);

    %Find if multiple cells have same minimum value:
    ind = find(deltaHCurr == val);
    
    %Use sink-filled DEM and don't let water route to same cell
    ind(ind == 5) = [];

    %Multiple indices:
    if numel(ind) == 3
        ind = ind(2);
    elseif numel(ind) == 2
        ind = ind(1);
    end
    
    %Write FDR output:
    fdrO(r,c) = fdr_val_arc(ind);
end

% indSink = find(~isnan(sinkBl));
% for ii = 1 : numel(indSink)
%     
% end




if regexpbl(strFrame,'in')
    %First and last rows:
    for jj = 1 : 2
        if jj == 1
            ri = 1;
        else
            ri = length(dem(:,1));
        end

        for ii = 2 : length(dem(1,:))-1 %Loop over columns
            if jj == 1
                deltaHCurr = reshape([nan(1,3); dem(1:2,ii-1:ii+1)],[],1) - dem(1,ii);
            else
                deltaHCurr = reshape([dem(end-1:end,ii-1:ii+1); nan(1,3)],[],1) - dem(end,ii);
            end
            
            deltaHCurr(indDiag) = deltaHCurr(indDiag) / facDiag;
            
            %Find steepest descent
            [val,~] = min(deltaHCurr);

            %Find if multiple cells have same minimum value:
            ind = find(deltaHCurr == val);

            %Write FDR output:
            fdrO(ri,ii) = fdr_val_arc(ind);
        end
    end

    %First and last columns:
    for ii = 1 : 2
        if ii == 1
            ci = 1;
        else
            ci = length(dem(1,:));
        end

        for jj = 2 : length(dem(:,1))-1 %Loop over columns
            if ii == 1
                deltaHCurr = reshape([nan(3,1), dem(jj-1:jj+1,1:2)],[],1) - dem(jj,1);
            else
                deltaHCurr = reshape([dem(jj-1:jj+1,end-1:end), nan(3,1)],[],1) - dem(jj,1);
            end
            
            deltaHCurr(indDiag) = deltaHCurr(indDiag) / facDiag;
            
            %Find steepest descent
            [val,~] = min(deltaHCurr);

            %Find if multiple cells have same minimum value:
            ind = find(deltaHCurr == val);

            %Write FDR output:
            fdrO(jj,ci) = fdr_val_arc(ind);
        end
    end
elseif regexpbl(strFrame,'out')
    for jj = 2 : length(dem(:,1))-1
        fdrO(jj,1) = fdr_val_arc(2); 
        fdrO(jj,end) = fdr_val_arc(8);     
    end
    for ii = 2 : length(dem(1,:))-1
        fdrO(1,ii) = fdr_val_arc(4); 
        fdrO(end,ii) = fdr_val_arc(6);     
    end
    fdrO(1,1) = fdr_val_arc(1);
    fdrO(end,1) = fdr_val_arc(3);
    fdrO(1,end) = fdr_val_arc(7);
    fdrO(end,end) = fdr_val_arc(9);
end
