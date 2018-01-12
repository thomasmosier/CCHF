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


function arrayDate = str2date(strDate)

%Input string format must be year/month/day/hour (with any symbol
%seperating)

strSymb = '[-\/.,;:|_]';

nDate = numel(regexpi(strDate{1},strSymb)) + 1;

arrayDate = nan(numel(strDate(:)),nDate);


for ii = 1 : numel(strDate(:))
    indCurr = regexpi(strDate{ii},strSymb);
    
    switch nDate
        case 1
           arrayDate(ii) = str2double(strDate{ii});
        case 2
            arrayDate(ii,:) = [...
                str2double(strDate{ii}(1:indCurr(1)-1)), ...
                str2double(strDate{ii}(indCurr(1)+1:end))];
        case 3
            arrayDate(ii,:) = [...
                str2double(strDate{ii}(1:indCurr(1)-1)), ...
                str2double(strDate{ii}(indCurr(1)+1:indCurr(2)-1)), ...
                str2double(strDate{ii}(indCurr(2)+1:end))];
        case 4
            arrayDate(ii,:) = [...
                str2double(strDate{ii}(1:indCurr(1)-1)), ...
                str2double(strDate{ii}(indCurr(1)+1:indCurr(2)-1)), ...
                str2double(strDate{ii}(indCurr(2)+1:indCurr(3)-1)), ...
                str2double(strDate{ii}(indCurr(3)+1:end))];
        otherwise
            error('str2date:nElements',['The date string has ' ...
                num2str(nDate) ' elements.  This case has not been programmed for.'])
    end
            
    
end