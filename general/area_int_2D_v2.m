function dataR = area_int_2D_v2(lonD,latD,data,lonR,latR, varargin)
% Copyright 2013, 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the GlobalClimateData Downscaling Package.
% 
% The GlobalClimateData Downscaling Package is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The Downscaling Package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%Fields used in input data:
%data
%latD
%lonD
%areaD
%latR
%lonR
%areaR
%varargin{1} = 'nansum' (nan's do not affect processing), 'sum' (output is
%nan if any input is nan

%Output:
%'reData' = sData area-weight integrated to 'sRef' grid.



%Varargin specifies 
if ~isempty(varargin(:)) && regexpbl(varargin{1}, 'nansum')
    type = 'nansum';
else
    type = 'sum';
end

%Make sure data does not have superfuluous index:
data = squeeze(data);
if length(size(data)) == 3
   error('area_int:dataD','Function only designed to work with 2D arrays'); 
end

%%Reorient data and reference vectors:
%Transpose reference longitude from row to column vector:
[rLonR, cLonR] = size(lonR);
if rLonR ~= 1 && cLonR == 1
    lonR = lonR';  
elseif rLonR ~= 1 && cLonR ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Transpose reference longitude from row to column vector:
[rLonD, cLonD] = size(lonD);
if rLonD ~= 1 && cLonD == 1
    lonD = lonD';  
elseif rLonD ~= 1 && cLonD ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Transpose reference latitude from colum to row vector:
[rLatR, cLatR] = size(latR);
if rLatR == 1 && cLatR ~= 1
    latR = latR';  
elseif rLatR ~= 1 && cLatR ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Transpose data latitude from colum to row vector:
[rLatD, cLatD] = size(latD);
if rLatD == 1 && cLatD ~= 1
    latD = latD';  
elseif rLatD ~= 1 && cLatD ~= 1
    error('area_int_1D:refMat',['The reference longitudes are in matrix '...
        'form but a vector is required.']);
end

%Flip and sort data so that it is in map orientation:
if ~issorted(lonD)
    [lonD, indLonD] = sort(lonD);
    data = data(:,indLonD);
end
if ~issorted(flipud(latD))
    [latD, indLatD] = sort(latD,'descend');
    data = data(indLatD,:);
end

%Sort reference coordinates so they're in map orientation:
indLonR = NaN;
indLatR = NaN;
if ~issorted(lonR)
    [lonR, indLonR] = sort(lonR);
end

%Flip reference latitude vector (result is greatest to least):
if ~issorted(flipud(latR)) 
    [latR, indLatR] = sort(latR,'descend');
end

%Check that size of lat and lon matches size of data
if length(latD) ~= length(data(:,1)) || length(lonD) ~= length(data(1,:))
    error('area_int_2D:dataMismatch',['The dimensions of the data do '...
        'not match the dimensions of the corresponding latitude and '...
        'longitude vectors.']);
end

%Initialize output array:
dataR = nan(length(latR),length(lonR));

%If gridding is same, but one is simply outside the other, fill in with
%nans and return
[nLatSame, indUseLatR, indUseLatD] = intersect(round2(latR,4),round2(latD,4));
[nLonSame, indUseLonR, indUseLonD] = intersect(round2(lonR,4),round2(lonD,4));
    nLatSame = numel(nLatSame);
    nLonSame = numel(nLonSame);

%Create counter to tell if grids are actually same but one is contained
%inside the other.  This is much more computationally efficient to solve.
same = 0;
if nLatSame == numel(latR) && nLonSame == numel(lonR) %Reference grid is subset of input grid
    dataR = data(indUseLatD,indUseLonD);
    same = 1;
elseif nLatSame == numel(latD) && nLonSame == numel(lonD) %Input grid is subset of reference grid
    dataR(indUseLatR,indUseLonR) = data;
    same = 1;
end
    
if same == 1
    %Flip output if reference were flipped:
    if ~isnan(indLatR)
        dataR = dataR(indLatR,:);
    end
    if ~isnan(indLonR)
        dataR = dataR(:,indLonR);
    end
    
    return
end


%Define grid at edge of cells:
lonEdgR  = box_edg(lonR);
latEdgR  = box_edg(latR);
lonEdgD  = box_edg(lonD);
latEdgD  = box_edg(latD);

% if latEdgD(1) > 90
%    warning('area_int_2D:latNset',['The current Northern latitude edge is ' ...
%        num2str(latEdgD(1)) ' which is being set to 90.'])
%    latEdgD(1) = 90;
% end
% if latEdgD(end) < -90
%     warning('area_int_2D:latSset',['The current Southern latitude edge is ' ...
%         num2str(latEdgD(end)) ' which is being set to -90.'])
%    latEdgD(end) = -90;
% end

%%Create common data grid:
latEdgC = sort(unique([latEdgD;latEdgR]),'descend');
dLatC = abs(diff(latEdgC));
latC = latEdgC(1:end-1)-0.5*dLatC;

lonEdgC = sort(unique([lonEdgD,lonEdgR]));
dLonC = abs(diff(lonEdgC));
lonC = lonEdgC(1:end-1)-0.5*dLonC;

%Initialize common data array and output array:
dataC = nan(length(latC),length(lonC));

%%Integrate data from common grid to reference grid:
%Calculate areas of grids:
areaD = area_geodata(lonEdgD,latEdgD,'e');
areaD(isnan(data)) = nan;
areaC = area_geodata(lonEdgC,latEdgC,'e');
areaR = area_geodata(lonEdgR,latEdgR,'e');

%%Place all data within common grid:
for ii = 1 : length(latC)
    for jj = 1 : length(lonC)
        indLatC = find(latEdgD(1:end-1) >= latC(ii) & latEdgD(2:end) <= latC(ii) );

        indLonC = find(lonEdgD(2:end) >= lonC(jj) & lonEdgD(1:end-1) <= lonC(jj) );
        
        if ~isempty(indLatC) && ~isempty(indLonC)
            if length(indLatC) == 2 || length(indLonC) == 2
                dataC(ii,jj) = sum2d(areaD(indLatC, indLonC).*data(indLatC, indLonC))./sum2d(areaD(indLatC, indLonC));
            else
                dataC(ii,jj) = data(indLatC, indLonC);
            end
        end
    end
end




areaC(isnan(dataC)) = nan;


%Find all data within each reference grid cell:
for ii = 1 : length(latR)
    for jj = 1 : length(lonR)
        indWR = find( round2(lonEdgC,3) == round2(lonEdgR(jj)  ,3));
        indER = find( round2(lonEdgC,3) == round2(lonEdgR(jj+1),3));
        indNR = find( round2(latEdgC,3) == round2(latEdgR(ii)  ,3));
        indSR = find( round2(latEdgC,3) == round2(latEdgR(ii+1),3));

        if ~isempty(indWR) && ~isempty(indER) && ~isempty(indNR) && ~isempty(indSR)
            if regexpbl(type,'nansum')
                dataR(ii,jj) = nansum(nansum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ nansum(nansum(areaC(indNR:indSR-1,indWR:indER-1)));
            else %If one or more cells of dataC are NaN, output will be NaN
                dataR(ii,jj) = sum(sum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ sum(sum(areaC(indNR:indSR-1,indWR:indER-1)));
            end
%             if regexpbl(type,'nansum')
%                 dataR(ii,jj) = nansum(nansum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ areaR(ii,jj);
%             else %If one or more cells of dataC are NaN, output will be NaN
%                 dataR(ii,jj) = sum(sum(areaC(indNR:indSR-1,indWR:indER-1).*dataC(indNR:indSR-1,indWR:indER-1) )) ./ areaR(ii,jj);
%             end
        else
            error('area_int_2D:noPtsC2R',['For ii = ' num2str(ii) ' and jj = ' num2str(jj) ', no points were selected.'])
        end
    end
end

%If areaR is 0, set dataR to 0 (occurs if latitude extends beyond +90 or
%-90)
dataR(areaR == 0) = 0;

%Flip output if reference were flipped:
if ~isnan(indLatR)
    dataR = dataR(indLatR,:);
end
if ~isnan(indLonR)
    dataR = dataR(:,indLonR);
end

% %Plot output:
% close all
% cSize = 40; %size of circles
% %Pick color bars:
% cMn = min([min(min(data)),min(min(dataR)),min(min(dataC))]);
% cMx = max([max(max(data)),max(max(dataR)),max(max(dataC))]);
% %Plot input data:
% [xD, yD] = meshgrid(lonD,latD);
% figure
% scatter(reshape(xD,1,[]),reshape(yD,1,[]),cSize,reshape(data,1,[]),'fill')
% colorbar
% caxis([cMn, cMx])
% %Plot common array:
% [xC, yC] = meshgrid(lonC,latC);
% figure
% scatter(reshape(xC,1,[]),reshape(yC,1,[]),cSize,reshape(dataC,1,[]),'fill')
% colorbar
% caxis([cMn, cMx])
% % %Plot common array (for reference grid region):
% % [xC, yC] = meshgrid(lonC(indWR:indER-1),latC(indNR:indSR-1));
% % figure
% % scatter(reshape(xC,1,[]),reshape(yC,1,[]),cSize,reshape(dataC(indWR:indER-1,indNR:indSR-1),1,[]),'fill')
% % colorbar
% % caxis([cMn, cMx])
% %Plot reference data:
% [xR, yR] = meshgrid(lonR,latR);
% figure
% scatter(reshape(xR,1,[]),reshape(yR,1,[]),cSize,reshape(dataR,1,[]),'fill')
% %pcolor(xR,yR,dataR)
% colorbar
% caxis([cMn, cMx])
