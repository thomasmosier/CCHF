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


function fileNm = file_nm(region,var,dateVec)



if isnumeric(dateVec)
    nDates = numel(dateVec(:,1));
    nPrec = numel(dateVec(1,:));
    fileNm = cell(nDates, 1);
    for jj = 1 : nDates
        dateStr = blanks(0);
        for ii = 1 : nPrec
            dateStr = [dateStr, num2str(dateVec(jj,ii))];
            if ii ~= nPrec
                dateStr = [dateStr, '.'];
            end
        end
        fileNm{jj} = [char(region) '_' char(var) '_' dateStr];
    end
elseif iscell(dateVec)
    nDates = numel(dateVec{1}(:,1));
    nPrec = numel(dateVec{1}(1,:));
    fileNm = cell(nDates, 1);
    for jj = 1 : nDates
        dateStr1 = blanks(0);
        for ii = 1 : nPrec
            dateStr1 = [dateStr1, num2str(dateVec{1}(jj,ii))];
            if ii ~= nPrec
                dateStr1 = [dateStr1, '.'];
            end
        end

        dateStr2 = blanks(0);
        for ii = 1 : nPrec
            dateStr2 = [dateStr2, num2str(dateVec{2}(jj,ii))];
            if ii ~= nPrec
                dateStr2 = [dateStr2, '.'];
            end
        end

        fileNm{jj} = [char(region) '_' char(var) '_' dateStr1 'thru' dateStr2];
    end
else
    error('fileNm:unknownType', ['Input of class ' class(dateVec) ' has not been programmed for.']);
end



