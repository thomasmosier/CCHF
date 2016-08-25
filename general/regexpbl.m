function bl = regexpbl(strfull, strsub, varargin)
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

%Function matches strings from 'str2' and 'str1'.
%Default is to return 1 if at least one entry is contained in both.
%'AND' conditional requires all strings from 'str2' be present in 'str1'.

%INPUTS:
%'str1': string or cell array of potential candidate strings
%'str2': string or cell array of strings to match with 'str2'
%'varargin': either 'or' or 'and' (default is 'or').

%OUTPUT:
%'bl': 1 or 0 determined by whether condition was met.

if isempty(strfull)
    bl = 0;
    return
end

chr1 = ischar(strfull);
cll1 = iscell(strfull);
chr2 = ischar(strsub);
cll2 = iscell(strsub);

if chr1 && chr2
    bl = ~isempty(regexpi(strfull, strsub));
    return
elseif ~chr1 && cll1 && chr2 && isempty(varargin)
    bl = 0;
    for ii = 1 : numel(strfull(:))
        if ischar(strfull{ii}) && ~isempty(regexpi(strfull{ii},strsub))
           bl = 1;
           return
        end
    end
    
    return
elseif chr1 && ~chr2 && cll2 && isempty(varargin)
    bl = 0;
    for ii = 1 : numel(strsub(:))
        if ischar(strsub{ii}) && ~isempty(regexpi(strfull,strsub{ii}))
           bl = 1;
           return
        end
    end
    
    return
end

%%HANDLE NORMAL CASES
%Varargin can be used to specify 'or' or 'and' comparison 
if ~isempty(varargin) && ~isempty(varargin{1})
    typ = varargin{1};
else
   typ = 'or'; 
end

if isempty(regexpi(typ,'and')) && isempty(regexpi(typ,'or'))
    warning('regexpbl:type',['Type is set to ' char(39) typ char(39) ...
        ', which is not a valid option.  It is therefore being set to ' ...
        char(39) 'or' char(39) '.']);
    typ = 'or';  
end

if iscell(strfull)
   nIt1 = length(strfull(:));
else
    nIt1 = 1;
    str1Temp = strfull;
    strfull = cell(1,1);
    strfull{1} = str1Temp;
end

if iscell(strsub)
   nIt2 = length(strsub(:));
else
    nIt2 = 1;
    str2Temp = strsub;
    strsub = cell(1,1);
    strsub{1} = str2Temp;
end

%Initialize output:
bl = 0;
blTemp = zeros(length(strsub(:)),1);

%Loop over all combinations of inputs:
for ii = 1 : nIt1
    for jj = 1 : nIt2
        if ~ischar(strfull{ii}) || ~ischar(strsub{jj})
            continue
        elseif ~isempty(regexpi(strfull{ii},strsub{jj})) %|| ~isempty(regexpi(strsub{jj},strfull{ii}))
            blTemp(jj) = 1;
            if strcmpi(typ,'or')
                bl = 1;
                return
            end
        end
    end
end

if strcmpi(typ,'and') && sum(blTemp) >= numel(strsub(:))
    bl = 1; 
end
