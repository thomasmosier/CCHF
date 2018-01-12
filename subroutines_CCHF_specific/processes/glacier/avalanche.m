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

function avalanche(sHydro, sMeta)
%Model avalanching of snow when slope is greater than critical angle
%This uses the glacier grid (instead of the main grid)


global sCryo

%These are not currently used:
%They could be used to track snow avalanching at the ice grid scale:
strIcF  = 'iceAvFr';
strIcT  = 'iceAvTo';

%These indicdes track how avalanches (calculated at ice scale) impact
%distribution of snow at main grid scale:
strIndMnF = 'AvIndMainFrom';
strIndMnT = 'AvIndMainTo';
strFrcMnF = 'AvFracMainFrom';
strFrcMnT = 'AvFracmainTo';

szMn = size(sHydro.dem);

% strFromMain = 'avFrMain';
% strToMain   = 'avToMain';

if ~isfield(sCryo,strIcF) ||  ~isfield(sCryo,strIcT) 
    %Angle is in degrees:
    angCrit = find_att(sMeta.global,'avalanche_angle');

    %Ensure critical angle expressed as positive number:
    if angCrit < 0
        angCrit = -angCrit;
    end
    %If expressed as degrees, convert to rise over run:
    if angCrit > 10
        angCrit = tand(angCrit);
    end
    
    %In flow direction grid, row is grid cell that flow is originating from
    %and column is cell the flow is going to
   	[indFrom, indTo] = find(sCryo.iceFdr ~= 0);  
    indFdr = sub2ind(size(sCryo.iceSlope), indFrom, indTo);
    indAval = find(sCryo.iceSlope(indFdr) <= -angCrit | sCryo.iceSlope(indFdr) >= angCrit); 
    indFrom = indFrom(indAval);
    indTo = indTo(indAval);
     
    %These are lists in ice grid frame linking grid cells that avalance and
    %those that are avalanched onto
    sCryo.(strIcF) = indFrom;
    sCryo.(strIcT)   = indTo;
    
    %Find avalanche cases where mass crosses main grid:
    indCross = find(sCryo.ice2main(indFrom) ~= sCryo.ice2main(indTo));
    
    %Only save the indices that cross main grid boundaries:
    %(THESE ARE ACTUALLY THE ONLY CELLS THAT NEED TO BE ADDRESSED EACH TIME STEP)
    indIcXMnF  = indFrom(indCross);
    indIcXMnT    = indTo(indCross);
    indMnFAll = full(sCryo.ice2main(indIcXMnF));
    indMnTAll   = full(sCryo.ice2main(indIcXMnT));
    
    
    %Calculate area-weighting factor, which will be used
    %to determine how change in ice grid cell impacts SWE balance of larger
    %grid cell:
    if ~regexpbl(sMeta.mode, {'calib','valid'})
        tic
        if numel(sCryo.iceSlope(:)) > 10000
            disp('Calculating avalance grids. This may take a few moments.');
        end
    end
    
    %LOOP OVER UNIQUE COMBINATIONS OF MAIN GRID FROM/TO
    [temp, ~, ~] = unique([indMnFAll, indMnTAll],'rows');
    sCryo.(strIndMnF) = temp(:,1);
    sCryo.(strIndMnT) = temp(:,2);
    
    sCryo.(strFrcMnF) = zeros(numel(sCryo.(strIndMnF)), 1);
    sCryo.(strFrcMnT) = zeros(numel(sCryo.(strIndMnF)), 1);
    
    szIce = size(sCryo.iceDem);
    
    for ii = 1 : numel(sCryo.(strIndMnF))
        %Find if there are multiple avalnche paths on ice grid
        %contributing to same main grid cell:
        indCurr = find(sCryo.(strIndMnF)(ii) == indMnFAll & sCryo.(strIndMnT)(ii) == indMnTAll);
        
        
        [rFromMain, cFromMain] = ind2sub(szMn, indMnFAll(indCurr));
        [  rToMain, cToMain]   = ind2sub(szMn, indMnTAll(indCurr));
        
        if numel(indCurr) > 1 && (any(diff(rFromMain) ~= 0) || any(diff(cFromMain) ~= 0) || any(diff(rToMain) ~= 0) || any(diff(cToMain) ~= 0))
           error('avalnche:multMainGrid','The main grid differs but should be unique.'); 
        else
            rFromMain = rFromMain(1);
            cFromMain = cFromMain(1);
            rToMain = rToMain(1);
            cToMain = cToMain(1);
        end
        
        if rFromMain < szMn(1) && cFromMain < szMn(2)
            areaMainFrom = area_geodata(...
                sHydro.lon(cFromMain:cFromMain+1)-0.5*diff(sHydro.lon(cFromMain:cFromMain+1)), ...
                sHydro.lat(rFromMain:rFromMain+1)+0.5*diff(sHydro.lat(rFromMain:rFromMain+1)),'e');
        elseif rFromMain < szMn(1) && cFromMain == szMn(2)
            areaMainFrom = area_geodata(...
                sHydro.lon(cFromMain-1:cFromMain)+0.5*diff(sHydro.lon(cFromMain-1:cFromMain)), ...
                sHydro.lat(rFromMain:rFromMain+1)+0.5*diff(sHydro.lat(rFromMain:rFromMain+1)),'e');
        elseif rFromMain == szMn(1) && cFromMain < szMn(2)
            areaMainFrom = area_geodata(...
                sHydro.lon(cFromMain:cFromMain+1)-0.5*diff(sHydro.lon(cFromMain:cFromMain+1)), ...
                sHydro.lat(rFromMain-1:rFromMain)+1.5*diff(sHydro.lat(rFromMain-1:rFromMain)),'e');
        elseif rFromMain == szMn(1) && cFromMain == szMn(2)
            areaMainFrom = area_geodata(...
                sHydro.lon(cFromMain-1:cFromMain)+0.5*diff(sHydro.lon(cFromMain-1:cFromMain)), ...
                sHydro.lat(rFromMain-1:rFromMain)+1.5*diff(sHydro.lat(rFromMain-1:rFromMain)),'e');
        end
        
        if rToMain < szMn(1) && cToMain < szMn(2)
            areaMainTo = area_geodata(...
                sHydro.lon(cToMain:cToMain+1)-0.5*diff(sHydro.lon(cToMain:cToMain+1)), ...
                sHydro.lat(rToMain:rToMain+1)+0.5*diff(sHydro.lat(rToMain:rToMain+1)),'e');
        elseif rToMain < szMn(1) && cToMain == szMn(2)
            areaMainTo = area_geodata(...
                sHydro.lon(cToMain-1:cToMain)+0.5*diff(sHydro.lon(cToMain-1:cToMain)), ...
                sHydro.lat(rToMain:rToMain+1)+0.5*diff(sHydro.lat(rToMain:rToMain+1)),'e');
        elseif rToMain == szMn(1) && cToMain < szMn(2)
            areaMainTo = area_geodata(...
                sHydro.lon(cToMain:cToMain+1)-0.5*diff(sHydro.lon(cToMain:cToMain+1)), ...
                sHydro.lat(rToMain-1:rToMain)+1.5*diff(sHydro.lat(rToMain-1:rToMain)),'e');
        elseif rToMain == szMn(1) && cToMain == szMn(2)
            areaMainTo = area_geodata(...
                sHydro.lon(cToMain-1:cToMain)+0.5*diff(sHydro.lon(cToMain-1:cToMain)), ...
                sHydro.lat(rToMain-1:rToMain)+1.5*diff(sHydro.lat(rToMain-1:rToMain)),'e');
        end
        
        
        areaIceFrom = zeros(numel(indCurr), 1);   
        areaIceTo = zeros(numel(indCurr), 1); 

        for jj = 1 : numel(indCurr)
            [ rFromIce, cFromIce]  = ind2sub( szIce, indIcXMnF(indCurr(jj)));
            [   rToIce, cToIce]    = ind2sub( szIce, indIcXMnT(indCurr(jj)));

            if rFromIce < szIce(1) && cFromIce < szIce(2)
                areaIceFrom(jj) = area_geodata(...
                    sCryo.iceLon(cFromIce:cFromIce+1)-0.5*diff(sCryo.iceLon(cFromIce:cFromIce+1)), ...
                    sCryo.iceLat(rFromIce:rFromIce+1)+0.5*diff(sCryo.iceLat(rFromIce:rFromIce+1)),'e');
            elseif rFromIce < szIce(1) && cFromIce == szIce(2)
                areaIceFrom(jj) = area_geodata(...
                    sCryo.iceLon(cFromIce-1:cFromIce)+0.5*diff(sCryo.iceLon(cFromIce-1:cFromIce)), ...
                    sCryo.iceLat(rFromIce:rFromIce+1)+0.5*diff(sCryo.iceLat(rFromIce:rFromIce+1)),'e');
            elseif rFromIce == szIce(1) && cFromIce < szIce(2)
                areaIceFrom(jj) = area_geodata(...
                    sCryo.iceLon(cFromIce:cFromIce+1)-0.5*diff(sCryo.iceLon(cFromIce:cFromIce+1)), ...
                    sCryo.iceLat(rFromIce-1:rFromIce)+1.5*diff(sCryo.iceLat(rFromIce-1:rFromIce)),'e');
            elseif rFromIce == szIce(1) && cFromIce == szIce(2)
                areaIceFrom(jj) = area_geodata(...
                    sCryo.iceLon(cFromIce-1:cFromIce)+0.5*diff(sCryo.iceLon(cFromIce-1:cFromIce)), ...
                    sCryo.iceLat(rFromIce-1:rFromIce)+1.5*diff(sCryo.iceLat(rFromIce-1:rFromIce)),'e');
            end

            if rToIce < szIce(1) && cToIce < szIce(2)
                areaIceTo(jj) = area_geodata(...
                    sCryo.iceLon(cToIce:cToIce+1)-0.5*diff(sCryo.iceLon(cToIce:cToIce+1)), ...
                    sCryo.iceLat(rToIce:rToIce+1)+0.5*diff(sCryo.iceLat(rToIce:rToIce+1)),'e');
            elseif rToIce < szIce(1) && cToIce == szIce(2)
                areaIceTo(jj) = area_geodata(...
                    sCryo.iceLon(cToIce-1:cToIce)+0.5*diff(sCryo.iceLon(cToIce-1:cToIce)), ...
                    sCryo.iceLat(rToIce:rToIce+1)+0.5*diff(sCryo.iceLat(rToIce:rToIce+1)),'e');
            elseif rToIce == szIce(1) && cToIce < szIce(2)
                areaIceTo(jj) = area_geodata(...
                    sCryo.iceLon(cToIce:cToIce+1)-0.5*diff(sCryo.iceLon(cToIce:cToIce+1)), ...
                    sCryo.iceLat(rToIce-1:rToIce)+1.5*diff(sCryo.iceLat(rToIce-1:rToIce)),'e');
            elseif rToIce == szIce(1) && cToIce == szIce(2)
                areaIceTo(jj) = area_geodata(...
                    sCryo.iceLon(cToIce-1:cToIce)+0.5*diff(sCryo.iceLon(cToIce-1:cToIce)), ...
                    sCryo.iceLat(rToIce-1:rToIce)+1.5*diff(sCryo.iceLat(rToIce-1:rToIce)),'e');
            end
        end

        sCryo.(strFrcMnF)(ii) = sum(areaIceFrom)/areaMainFrom;
        sCryo.(strFrcMnT)(ii) = sum(areaIceTo)/areaMainTo;
    end
    
    if ~regexpbl(sMeta.mode, {'calib','valid'})
        tAv = round2(toc/60,2);
        if tAv > 0.5
            disp(['It took ' num2str(tAv) ' minutes to calculate the avalanche grids.']);
        end
    end
end


%REDISTRIBUTE AVALANCHED SNOW AT MAIN GRID RESOLUTION:
sCryo.avsnw = sparse(szMn(1), szMn(2));

% strIndMnF = 'AvIndMainFrom';
% strIndMnT = 'AvIndMainTo';
% strFrcMnF = 'AvFracMainFrom';
% strFrcMnT = 'AvFracmainTo';

sCryo.avsnw(sCryo.(strIndMnF)) = -sCryo.(strFrcMnF).*sCryo.snw(sCryo.(strIndMnF));
sCryo.avsnw(sCryo.(strIndMnT)) = -sCryo.avsnw(sCryo.(strIndMnF));

indAval = find(sCryo.avsnw ~= 0);
sCryo.snw(indAval) = sCryo.snw(indAval) + full(sCryo.avsnw(indAval));

