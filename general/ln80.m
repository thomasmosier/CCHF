function str = ln80(str)
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


mxChar = 80;

if length(str) > mxChar
    indLn = 0;
    nIt = 1;

    while length(str) - indLn > mxChar && nIt*mxChar <= length(str)
        indSpc = regexpi(str,char(32));
        indLn = find(indSpc < nIt*mxChar, 1, 'last');
                nIt = nIt + 1;
        if numel(indLn) > 0
            indLn = indSpc(indLn);
            str = [str(1:indLn-1), char(10), str(indLn+1:end)];
        else
            continue
        end
    end
else
    return;
end