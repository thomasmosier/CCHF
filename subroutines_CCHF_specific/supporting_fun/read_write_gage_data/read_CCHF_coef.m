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


function [coef, hdr] = read_CCHF_coef(path)

[~, flNm, ~ ] = fileparts(path);
%Check to see if filename appears consistent with expected format:
if ~regexpbl(flNm, {'coefficients','coef'})
    warning('read_CCHF_coef:unexpectedNm',['The file ' char(39) flNm char(39) ' does not appear to be a coefficient file (based on naming).'])
end

fid = fopen(path,'rt');

hdr = cell(0,1);
cntLp = 1;
ii = 1;
while cntLp == 1
    hdr{ii,1}   = fgets(fid);
    if strcmpi(hdr{ii}, char(10)) || isempty(hdr{ii}) || strcmpi(hdr{ii}, [char(13) char(10)])
        rHdr = ii;
        break
    elseif hdr{ii,1} == -1
        break
    end 
    ii = ii + 1;
end
hdr = hdr(1:rHdr-1,:);

labels = fgets(fid);
nCol = numel(regexpi(labels, ',')) + 1;

dataStr = '%s';
for jj = 1 : nCol -1
    dataStr = [dataStr, '\t%s'];
end

cellRaw = textscan(fid, dataStr);

fclose(fid);

coef = cell(numel(cellRaw{1}), numel(cellRaw(1,:)));

for kk = 1 : numel(cellRaw{1})
    for ll = 1 : nCol
        coef{kk,ll} = cellRaw{ll}(kk);
        if iscell(coef{kk,ll})
            if ischar(coef{kk,ll}{1})
                coef{kk,ll} = coef{kk,ll}{1};
                if sum(isletter(coef{kk,ll})) == 0
                    coef{kk,ll} = str2double(coef{kk,ll});
                end
            end

        end
    end
    
    indComma = regexpi(coef{kk,1},',');
    if ~isempty(indComma)
        
        coef{kk,1} = coef{kk,1}(1:indComma(1)-1);
    end
end