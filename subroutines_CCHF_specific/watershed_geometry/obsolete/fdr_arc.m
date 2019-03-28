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

function fdrO = fdr_arc(dem,lat,lon,varargin)
%Important info on ArcMap algorithm at:
%http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//009z00000063000000.htm

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
facDiag = 1.4142;

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
            
            deltaHCurr(indHor) = deltaHCurr(indHor) / facHor;
            deltaHCurr(indVer) = deltaHCurr(indVer) / facVer;
            
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
            
            deltaHCurr(indHor) = deltaHCurr(indHor) / facHor;
            deltaHCurr(indVer) = deltaHCurr(indVer) / facVer;
            
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
