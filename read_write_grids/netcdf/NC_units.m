function [units, varargout] = NC_units(att)
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


%INPUTS: 
%att = 2 column cell-array of data attributes, where one of the
%entries in the first column is "units"
%OUTPUTS:
%units = unit string
%varagout = row of units string

indUnit = strcmpi(att,'units');

if ~isempty(indUnit) && sum2d(indUnit) == 1
    [rowUnt, colUnt] = find(indUnit == 1);
    infoUnt = att{rowUnt, colUnt+1};

    %If the units str is strictly non-numeric, take all of it, otherwise take
    %only non-numeric component:
    gcmUnitsInd = regexpi(infoUnt,'\d');
    if ~isempty(gcmUnitsInd)
        units = infoUnt(1:gcmUnitsInd(1)-2);
    else
        units = infoUnt;
    end
elseif ~isempty(indUnit) && sum2d(indUnit) > 1
    error('NC_Units:multEntries',['Multiple entries were found that '...
        'have the term ' char(39) 'units' char(39) '. Unable to select '...
        'correct entry.']);
else
    units = nan;
    rowUnt = nan;
end


if nargout == 2
    varargout{1} = rowUnt;
end