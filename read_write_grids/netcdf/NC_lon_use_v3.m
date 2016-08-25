function [lonTr, indLonUse] = NC_lon_use_v3(lonOrg, attLon, sMeta)
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


%%Find indices for longitude of interest in NetCDF or GRB file:
indLonUnt = strcmpi(attLon,'units');
[rowLon, colLon] = find(indLonUnt == 1);

%Checking units may or may not be important:
lonUnits = attLon{rowLon, colLon+1};

flagLon = zeros(2,1);

%Determine frame to use (around coordinates bounds)
if isfield(sMeta,'frame') 
    frame = sMeta.frame;
else
    frame = 2; 
end

%Find longitude indices that are within the desired region:
%Identifies the case where 'degrees East' are used:
if ~isempty(regexpi(lonUnits,'deg')) && ~isempty(regexpi(lonUnits,'east'))
    if max(lonOrg) > 180 && min(lonOrg) >= 0 %Longitude goes from 0 to 360 degrees east; otherwise assume -180 -> 180
        %Transform to -180 -> 0 -> 180
        indOrg180 = find(lonOrg > 180,1, 'first');
        lonTr = [lonOrg(indOrg180:end)-360, lonOrg(1:indOrg180-1)];

        if isfield(sMeta,'crd')
            lonTrEdg  = box_edg(lonTr);

            %Find indices corresponding to bounds of selection (in transformed 
            %coordinates):
            indW = find(lonTrEdg(1:end-1) < sMeta.crd(1), 1, 'last');
            if isempty(indW) && 2*lonTrEdg(1) - lonTrEdg(2) <= -179 && 2*lonTrEdg(end) - lonTrEdg(end-1) > 179   %Data extends beyond -180
               indW = find(lonTrEdg(1:end-1) -360 < sMeta.crd(1), 1, 'last') - numel(lonOrg) - 1;
            elseif isempty(indW)
                indW = 1;
                flagLon(1) = 1;
            end
            indE = find(lonTrEdg(2:end) > sMeta.crd(2), 1, 'first');
            if isempty(indE) && 2*lonTrEdg(1) - lonTrEdg(2) <= -179 && 2*lonTrEdg(end) - lonTrEdg(end-1) > 179 %Data extends beyond 180
               indE = find(lonTrEdg(2:end) + 360 > sMeta.crd(2), 1, 'first') + numel(lonOrg);
            elseif isempty(indE)
                indE = numel(lonOrg);
                flagLon(2) = 1;
            end

            %Identify indice where transformed coordinates are closest to 0:
            indTr0 = find(lonTr >= 0, 1, 'first');

            %Rearange indices based on GCM (original) orientation:
                %Load double wide boarder.
                %On the "exterior wide" need two more than on the interior side
                %(i.e. 3 vs. 1)
            if any(flagLon)
                indLonUse = (indW : indE);
                warning('NC_Lon_use:small_grid','This case has not been error checked.');
            else
                if indW < indTr0 && indE < indTr0 %Both negative
                    indLonUse = (indOrg180 + indW - 1 - frame : indOrg180 + indE - 1 + frame);
                elseif indW < indTr0 && indE >= indTr0 %western neg and eastern positive
                    indLonUse = [(indOrg180 + indW - 1 - frame : numel(lonOrg)), (1 : indE - indOrg180 + 3 + frame)]; %Empirically this must be 5;
                elseif indW >= indTr0 && indE >= indTr0 %Both positive
                    indLonUse = (indW - indOrg180 - 1 + frame : indE - indOrg180 + 3 + frame);
                else
                    error('NC_lon_use:lonOrientation','Case unexpected and not coded for.');
                end
            end

            %If longitude on edge of grid being used, must wrap:
            indLonUse = bnd_wrap(indLonUse, numel(lonOrg));

            %Select longitude values for output:
            lonTr = lonOrg(indLonUse);
        else
            indLonUse = (1:numel(lonOrg));
        end
        lonTr( lonTr > 180 ) = lonTr( lonTr > 180 ) - 360;
    else %Longitude goes from -180 to 180 degrees east
        if isfield(sMeta,'crd')
            lonOrgEdg  = box_edg(lonOrg);

            indW = find(lonOrgEdg(1:end-1) < sMeta.crd(1), 1, 'last');
            if isempty(indW) && 2*lonOrgEdg(1) - lonOrgEdg(2) <= -179 && 2*lonOrgEdg(end) - lonOrgEdg(end-1) > 179   %Data extends beyond -180
                indW = find(lonOrgEdg(1:end-1) -360 < sMeta.crd(1), 1, 'last') - numel(lonOrg) - 1;
            elseif isempty(indW)
                indW = 1;
                flagLon(1) = 1;
            end
            indE = find(lonOrgEdg(2:end) > sMeta.crd(2), 1, 'first');
            if isempty(indE) && 2*lonOrgEdg(1) - lonOrgEdg(2) <= -179 && 2*lonOrgEdg(end) - lonOrgEdg(end-1) > 179 %Data extends beyond 180
               indE = find(lonOrgEdg(2:end) + 360 > sMeta.crd(2), 1, 'first') + numel(lonOrg); 
            elseif isempty(indE)
                indE = numel(lonOrg);
                flagLon(2) = 1;
            end

            if ~any(flagLon)
                indLonUse = (indW - frame : indE + frame);
            elseif flagLon(1) && ~flagLon(2)
                indLonUse = (indW : indE + frame);
            elseif ~flagLon(1) && flagLon(2)
                 indLonUse = (indW - frame : indE);
            else
                indLonUse = (indW : indE);
            end

            %If longitude on edge of grid being used, must wrap:
            indLonUse = bnd_wrap(indLonUse, numel(lonOrg));
        else
            indLonUse = (1:numel(lonOrg));
        end
        %Select longitude values for output:
        lonTr = lonOrg(indLonUse);
%     else
%         keyboard
%         error('gcm_load:coordOffset',['An unknown longitudinal offset '...
%             'is used in the current GCM time series.' char(10) ...
%             'There is likely nothing wrong except that this case must be coded for.']);
    end
    
            %Bnd Wrap may introduce discontinuity, so remove unnatural jumps
            indRem = find(abs(diff(lonOrg(indLonUse)))/(5*mean(abs(diff(lonOrg(indLonUse))))) > 1);
            if ~isempty(indRem)
                if diff([lonOrg(indLonUse(indRem(1))),lonOrg(indLonUse(indRem(1)+1))]) < 0
                    if indRem(1) > 0.5*numel(indLonUse)
                        indLonUse = indLonUse(1:indRem(1));
                    else
                        indLonUse = indLonUse(indRem(1):end);
                    end
                end
            end

    %If lonTr wraps backwards across -180, adjust values so they extend
    %beyond -180 (necessary for interpolation on flat surface)
    indNeg = find(diff(lonTr) < 0);
    if ~isempty(indNeg)
        if numel(indNeg) == 1
            if indNeg < 0.5*numel(lonTr)
                lonTr(1:indNeg) = lonTr(1:indNeg) - 360;
            else
                lonTr(indNeg+1:end) = lonTr(indNeg+1:end) + 360;
            end      
        elseif numel(indNeg) == 2
            lonTr(1:indNeg(1)) = lonTr(1:indNeg(1)) - 360;
            lonTr(indNeg(2)+1:end) = lonTr(indNeg(2)+1:end) + 360;
        else
            error('NC_lon_use:wrap2Ind','This case not coded for.')
        end
    end
else
    error('gcm_load:coordSys','An unknown longitudinal coordinate system is used in the current GCM time series.');
end