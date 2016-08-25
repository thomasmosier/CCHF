function [gcmTs, attData, attGen, nm] = NC_data_use_v2(path, nm, tInd, latInd, lonInd)
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



%Get dimension names from netCDf file:
dimRaw = ncinfo(path,nm);
dimNC = cell(numel(dimRaw.Dimensions(:)), 1);
for ii = 1 : numel(dimRaw.Dimensions(:))
    dimNC{ii} = dimRaw.Dimensions(ii).Name;
end


%Define desired order:
indD = {tInd, latInd, lonInd};
dimOrdered = {'time';'latitude';'longitude'};
if ~any(strcmpi(dimOrdered{2},dimNC))
    dimOrdered{2} = 'lat';  
end
if ~any(strcmpi(dimOrdered{2},dimNC))
   error('NC_data_use:latNm','Cannot find correct name of latitude dimension.');
end
if ~any(strcmpi(dimOrdered{3},dimNC))
    dimOrdered{3} = 'lon';  
end
if ~any(strcmpi(dimOrdered{3},dimNC))
   error('NC_data_use:lonNm','Cannot find correct name of longitude dimension.');
end


sortDim = (1:3);
%If dimensions not ordered, reorder:
if ~all(strcmpi(dimOrdered,dimNC))
    for ii = 1 : numel(dimOrdered)
        for jj = 1 : numel(dimNC)
            if strcmpi(dimOrdered{ii},dimNC{jj})
                sortDim(ii) = jj;
            end
        end
    end
end


%If data out of order, must load in pieces:
ld = zeros(3,1);
if ~issorted(indD{sortDim(1)}) || any(diff(indD{sortDim(1)}) ~= 1)
    ld(1) = 1;
end
if ~issorted(indD{sortDim(2)}) || any(diff(indD{sortDim(2)}) ~= 1)
    ld(2) = 1;
end
if ~issorted(indD{sortDim(3)}) || any(diff(indD{sortDim(3)}) ~= 1)
    ld(3) = 1;
end


if sum(ld) == 0 %If no pieces need to be re-ordered:
    gcmTs = single(ncread(path, nm, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
elseif sum(ld) == 1
    gcmTs = zeros([numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), numel(indD{sortDim(3)})],'single');
    if ld(1) == 1
        [~, indSort1] = sort(indD{sortDim(1)});
        if all(diff(indSort1) ==1)
            try
                gcmTs = single(ncread(path, nm, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
                catch ME
                    if (strcmp(ME.identifier,'MATLAB:nomem'))
                        indX = round(0.5*indD{sortDim(1)}(1));
                        gcmTs(1:indX,:,:) = single(ncread(path, nm, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [indX, numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
                        gcmTs(indX:end,:,:) = single(ncread(path, nm, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}) - indX, numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
                    end
            end
            gcmTs = gcmTs(indSort1,:,:);
        else
            for ii = 1 : numel(indD{sortDim(1)})
                gcmTs(ii,:,:) = single(ncread(path, nm, [indD{sortDim(1)}(ii), indD{sortDim(2)}(1), indD{sortDim(3)}(1)], [1, numel(indD{sortDim(2)}), numel(indD{sortDim(3)})]));
            end
        end
    elseif ld(2) == 1
        for ii = 1 : numel(indD{sortDim(2)})
            gcmTs(:,ii,:) = single(ncread(path, nm, [indD{sortDim(1)}(1), indD{sortDim(2)}(ii), indD{sortDim(3)}(1)], [numel(indD{sortDim(1)}), 1, numel(indD{sortDim(3)})]));
        end
    elseif ld(3) == 1
        for ii = 1 : numel(indD{sortDim(3)})
            gcmTs(:,:,ii) = single(ncread(path, nm, [indD{sortDim(1)}(1), indD{sortDim(2)}(1), indD{sortDim(3)}(ii)], [numel(indD{sortDim(1)}), numel(indD{sortDim(2)}), 1]));
        end
    end
    
elseif sum(ld) > 1
    gcmTs = zeros([numel(indD{1}), numel(indD{2}), numel(indD{3})],'single');
    warning('NC_data_use:d2sep',[num2str(sum(ld)) ' dimensions are ' ...
        'out of order. This increases NetCDF read time.']);
    for ii = 1 : numel(indD{1})
        for jj = 1 : numel(indD{2})
            for kk = 1 : numel(indD{3})
                vecCurr = [indD{1}(ii), indD{2}(jj), indD{3}(kk)];
                gcmTs(ii,jj,kk) = single(ncread(path, nm, vecCurr(sortDim), [1, 1, 1]));
            end
        end
    end
    sortDim = (1:3);
end
    


%Reorder dimenions if neccesarry:
if ~issorted(sortDim)
    gcmTs = permute(gcmTs,sortDim);
end

%Load general attributes:
attGen = ncinfo(path);
attGen = squeeze(struct2cell(attGen.Attributes))';

%Load data attributes:
attData = ncinfo(path,nm);
attData = squeeze(struct2cell(attData.Attributes))';

end