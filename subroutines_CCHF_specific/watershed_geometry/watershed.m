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

function sHydro = watershed(sHydro,varargin)
%Calculates fdr, fac, and slope (all neighbors and only fdr), aspect, and
%distance between neighboring centroids, and area of each cell:



% rEarth = find_att(varargin{1}.global,'radius_Earth');

rEarth =  6371000;

numDEM = numel(sHydro.dem);

if isfield(sHydro, 'lon')
    varLon = 'lon';
elseif isfield(sHydro, 'longitude')
    varLon = 'longitude';
else
    error('watershed:uknownlongvariable','longitude variable not found.');
end

if isfield(sHydro, 'lat')
    varLat = 'lat';
elseif isfield(sHydro, 'latitude')
    varLat = 'latitude';
else
    error('watershed:uknownlongvariable','longitude variable not found.');
end
    
    
%%FIND ACCUMULATED UPSTREAM AREA, FLOW-DIRECTION, AND SLOPE:
[X, Y] = meshgrid(sHydro.(varLon), sHydro.(varLat));
%Change NaN values in DEM temporarily to very high values:
% nanVal = 9999999;
% sHydro.dem(isnan(sHydro.dem)) = nanVal;
%SLOPE = rise / run (unitless; and not angle)

% sHydro.dem(sHydro.dem == nanVal) = NaN;

szDem = size(sHydro.dem);
%Translate ESRI FDR to sparse matrix:
if isfield(sHydro,'fdrESRI')
    sHydro.fdr = sparse(numel(sHydro.dem),numel(sHydro.dem));
    
    valFDR = [1,2,4,8,16,32,64,128];
    for ii = 1 : numel(sHydro.dem)
        valC = sHydro.fdrESRI(ii);
        
        if ~any(valC == valFDR)
            flagC = 1;  cntr = 0;
            while flagC == 1
                indC = find(valC >= valFDR, 1, 'last');
                if isempty(indC)
                	flagC = 0;
                else
                    cntr = cntr + 1;
                    valC(cntr) = valFDR(indC);
                end
            end
        end
        
        [rC, cC] = ind2sub(szDem, ii);
        for jj = 1 : numel(valC)
            if ~isnan(valC(jj))
                switch valC(jj)
                    case 1 %East
                        r = 0;
                        c = 1;
                    case 2 %Southeast
                        r = 1;
                        c = 1;
                    case 4 %South
                        r = 1;
                        c = 0;
                    case 8 %Southwest
                        r = 1;
                        c = -1;
                    case 16 %West
                        r = 0;
                        c = -1;
                    case 32 %Northwest
                        r = -1;
                        c = -1;
                    case 64 %North
                        r = -1;
                        c = 0;
                    case 128 %Northeast
                        r = -1;
                        c = 1;
                end
            else
                r = NaN;
                c = NaN;
            end
            
            rAbs = rC + r;
            cAbs = cC + c;
            if ~(isnan(rAbs) || rAbs < 1 || rAbs > szDem(1) || cAbs < 1 || cAbs > szDem(2))
                sHydro.fdr(ii,sub2ind(szDem, rAbs, cAbs)) = 1;
            end
        end
    end
    
    [~,~,sHydro.slope] = wflowacc(X,Y,sHydro.dem,'type','single','edges','open','coord','geographic');
else
    % [sHydro.fac, sHydro.fdr, sHydro.slope] = wflowacc(X,Y,sHydro.dem,'type','single','edges','open');
    % [sHydro.fac, sHydro.fdr, sHydro.slope] = wflowacc(X,Y,sHydro.dem,'type','single','edges','open','coord','geographic');
    [sHydro.fac, sHydro.fdr, sHydro.slope] = wflowacc(X,Y,sHydro.dem,'type','single','edges','open','coord','geographic');
    %In FDR: row is cell of origin and col is downstream cell
    %     warning('watershed:ESRINeeded','No fdrESRI field is present');
end

% %SLOPE CORRECTION:
% %SCALE slope by constant facTo(ii)r (to account for DEM resolution: coarser resolution decreases slope)
% sHydro.slope = 3.5*sHydro.slope;
% %BEST
% %Wu, S., Li, J., & Huang, G. H. (2008). A study on DEM-derived primary 
% %topographic attributes for hydrologic applications: sensitivity to 
% %elevation data resolution. Applied Geography, 28(3), 210–223.
% 
% %NOT AS GOOD
% % Thompson, J. A., Bell, J. C., & Butler, C. A. (2001). Digital elevation 
% % model resolution: effects on terrain attribute calculation and 
% % quantitative soil-landscape modeling. Geoderma, 100(1), 67–89.

%Add field called 'fdrSlope', which is the slope in the flow direction
%(used for a few subroutines in CCHF)
sHydro.slopeFdr = -single(full(reshape(sum(sHydro.fdr.*sHydro.slope, 2), szDem))); %slope at index (i,j) is slope from this cell to next cell in flowpath
sHydro.slopeFdrA = atand(sHydro.slopeFdr);


if ~isempty(varargin) && isfield(varargin{1}, 'module') && regexpbl(varargin{1}.module,'Muskingum')
    %%FIND ORDER THAT FLOW MUST BE CALCULATED IN:
    %Define initial order based on elevation ranking:
    [~, flowOrder] = sort(reshape(sHydro.dem,[],1),'descend');
    [flowOrder, ~] = find(sparse((1:numDEM),flowOrder,1,numDEM,numDEM));
    sHydro.flowOrder = uint16(reshape(flowOrder, numel(sHydro.dem(:,1)), numel(sHydro.dem(1,:))));

    %At cells where flow order is backwards, swap:
    [flowFrom, flowTo] = find(sHydro.fdr == 1);

    flagSwap = 1; iterSwap = 0;
    while flagSwap == 1
        iterSwap = iterSwap + 1;

        iSwap = (sHydro.flowOrder(flowFrom) -  sHydro.flowOrder(flowTo)) > 0;

        if sum(iSwap) > 0
            flowTemp = sHydro.flowOrder(flowFrom(iSwap));
            sHydro.flowOrder(flowFrom(iSwap)) = sHydro.flowOrder(flowTo(iSwap));
            sHydro.flowOrder(  flowTo(iSwap)) = flowTemp;
        else
            flagSwap = 0;
        end

        if iterSwap > 100
           break 
        end
    end
end



%%CALCULATE DISTANCE BETWEEN ADJACENT CELLS:
%Populate sHydro.dl as sparse matrix with same properties as sHydro.fdr
%and sHydro.slope:

%A few functions use a grid of latitudes:
[X, sHydro.latGrid] = meshgrid(sHydro.(varLon),sHydro.(varLat));
sHydro.dl = haversine_neighbors(X,sHydro.latGrid, rEarth,'sparse');

%Make field that is distance along flowpath:
sHydro.dlFdr = single(full(reshape(sum(sHydro.fdr.*sHydro.dl, 2), szDem)));
sHydro.dlFdr(sHydro.dlFdr == 0) = nan;


%%CALCULATE ASPECT:
sHydro.aspectFdr = nan(szDem,'single'); %-1 is value assigned to flat surfaces

[indFrom, indTo] = find(sHydro.fdr == 1);
% %Sort:
% [indFrom, orderFrom(ii)] = sort(indFrom,'ascend');
% indTo = indTo(orderFrom(ii));

%Convert indices to row and columns in DEM:
[rFrom, cFrom] = ind2sub(szDem, indFrom);
[  rTo,   cTo] = ind2sub(szDem,   indTo);

for ii = 1 : numel(rFrom)
    if rFrom(ii) - rTo(ii) ==  0 && cTo(ii) - cFrom(ii) ==  0
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = -1;
    elseif rFrom(ii) - rTo(ii) ==  1 && cTo(ii) - cFrom(ii) ==  0
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) =  0; %North
    elseif rFrom(ii) - rTo(ii) == -1 && cTo(ii) - cFrom(ii) ==  0
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = 180; %South
    elseif rFrom(ii) - rTo(ii) ==  0 && cTo(ii) - cFrom(ii) ==  1
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = 270; %West
    elseif rFrom(ii) - rTo(ii) ==  0 && cTo(ii) - cFrom(ii) == -1
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = 90; %East
    elseif rFrom(ii) - rTo(ii) ==  1 && cTo(ii) - cFrom(ii) == -1
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = 45; %Northeast
    elseif rFrom(ii) - rTo(ii) == -1 && cTo(ii) - cFrom(ii) == -1
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = 135; %SouthEast
    elseif rFrom(ii) - rTo(ii) ==  1 && cTo(ii) - cFrom(ii) ==  1
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = 315; %Northwest
    elseif rFrom(ii) - rTo(ii) == -1 && cTo(ii) - cFrom(ii) ==  1
        sHydro.aspectFdr(rFrom(ii), cFrom(ii)) = 225; %Southwest
    end
end

%Set all other values to -1 (flat). This typically occurs for sinks
sHydro.aspectFdr(isnan(sHydro.aspectFdr)) = -1;

%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) ==  0 & cTo(ii) - cFrom(ii) ==  0) =  -1; %Flat
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) ==  1 & cTo(ii) - cFrom(ii) ==  0) =   0; %North
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) == -1 & cTo(ii) - cFrom(ii) ==  0) = 180; %South
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) ==  0 & cTo(ii) - cFrom(ii) ==  1) = 270; %West
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) ==  0 & cTo(ii) - cFrom(ii) == -1) =  90; %East
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) ==  1 & cTo(ii) - cFrom(ii) == -1) =  45; %Northeast
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) == -1 & cTo(ii) - cFrom(ii) == -1) = 135; %SouthEast
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) ==  1 & cTo(ii) - cFrom(ii) ==  1) = 315; %Northwest
%     sHydro.aspectFdr(rFrom(ii) - rTo(ii) == -1 & cTo(ii) - cFrom(ii) ==  1) = 225; %Southwest


% 
% indAspPart = sub2ind([3,3], rFrom(ii) - rTo(ii) + 2, cTo(ii) - cFrom(ii) + 2);
% sHydro.aspectFdr(indFrom) = indAspPart;
% 
% sHydro.aspectFdr(sHydro.aspectFdr == 1) = 315;
% sHydro.aspectFdr(sHydro.aspectFdr == 2) = 270;
% sHydro.aspectFdr(sHydro.aspectFdr == 3) = 225;
% sHydro.aspectFdr(sHydro.aspectFdr == 4) =   0;
% sHydro.aspectFdr(sHydro.aspectFdr == 5) = 180;
% sHydro.aspectFdr(sHydro.aspectFdr == 7) =  45;
% sHydro.aspectFdr(sHydro.aspectFdr == 8) =  90;
% sHydro.aspectFdr(sHydro.aspectFdr == 9) = 135;



%%CALCULATE AREA:
sHydro.area = area_geodata(sHydro.(varLon),sHydro.(varLat));