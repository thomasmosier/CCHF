function [gcmLat, gcmLatUseInd] = NC_lat_use_v3(gcmLat, attLat, sMeta)
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


%Find indices for latitude of interest in NetCDF or GRB file:
latUnitInd = strcmpi(attLat,'units');
[latRow, latCol] = find(latUnitInd == 1);

%Find lat edges of GCM grid:
latEdg = box_edg(gcmLat);

%Check units:
latInfo = attLat{latRow, latCol+1};

%Determine frame to use (around coordinates bounds)
if isfield(sMeta,'frame') 
    frame = sMeta.frame;
else
    frame = 2; 
end

if isfield(sMeta,'crd')
    %Find indices of file corresponding to desired region:
    if ~isempty(regexpi(latInfo,'deg')) && ~isempty(regexpi(latInfo,'north'))
        if gcmLat(2) > gcmLat(1) %Top of array is South
            %Find latitude indices that are within the desired region:
            indS = find(latEdg(1:end-1) <= sMeta.crd(3), 1, 'last');
            if isempty(indS)
                indS = 1; 
            end
            indN = find(latEdg(2:end) >= sMeta.crd(4), 1, 'first');
            if isempty(indN)
                indN = numel(gcmLat); 
            end

            %Indices of GCM cells to use, including double wide border
            gcmLatUseInd = (indS - frame : indN - 1 + frame); %Empirically this provides two centroids outside of requested bounds
        elseif gcmLat(1) > gcmLat(2) %Top of array is North
            %Find latitude indices that are within the desired region:
            indN = find(latEdg(1:end-1) >= sMeta.crd(4), 1, 'last');
            if isempty(indN)
                indN = 1; 
            end
            indS = find(latEdg(2:end) <= sMeta.crd(3), 1, 'first');
            if isempty(indS)
                indS = numel(gcmLat); 
            end
            
            %Indices of GCM cells to use, including border
            gcmLatUseInd = (indN + 1 - frame : indS + frame); %Empirically this provides two centroids outside of requested bounds
        else
            error('NC_lat_use:noChange',['The spacing between latitudes '...
                'in the current NetCDF dataset appears to be 0.  This case has '...
                'not been programmed.']);
        end

    else %No other conditionals written for different units
        error('NC_lat_use:coordSys',['The latitudinal coordinate system '...
            'in the current NetCDF dataset is ' latUnits ', which has not '...
            'been coded for.']);
    end
else
    gcmLatUseInd = (1:numel(gcmLat));
end

%Remove any latitudes above or below North or South Poles, respectively.
indLes = find(gcmLatUseInd < 1);
if ~isempty(indLes)
    gcmLatUseInd(indLes) = 1;
end
indGre = find(gcmLatUseInd > numel(gcmLat));
if ~isempty(indGre)
    gcmLatUseInd(indGre) = numel(gcmLat);
end
%Keep only unique values:
gcmLatUseInd = unique(gcmLatUseInd,'stable');
%List latitudes within target region:
gcmLat = gcmLat(gcmLatUseInd); 
