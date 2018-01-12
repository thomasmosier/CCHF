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

function sGage = import_Arc_data(path, lonGrid, latGrid)


[~, nmFile, ~ ] = fileparts(path);

szGrid = [numel(latGrid), numel(lonGrid)];

%Read data from Arc attribute table copied into Excel file.
[dataGage, lonData, latData, date] = read_Arc_data(path);

%%Place data into substructure format:
%Initialize structure
sGage = struct;
%Find unique lat, lon points in data:
uniqPairs = unique([lonData,latData],'rows');

%loop over unique grid points to write observation data to sub-structures:
for ii = 1 : numel(uniqPairs(:,1))
    indCurr = find(lonData == uniqPairs(ii,1) & latData == uniqPairs(ii,2));
       
    [~, gageRow] = min(abs(uniqPairs(ii,2)-latGrid));
    [~, gageCol] = min(abs(uniqPairs(ii,1)-lonGrid));
    if uniqPairs(ii,1) > latGrid(1) - 0.5*abs(diff(latGrid(1:2))) || uniqPairs(ii,1) < latGrid(end) + 0.5*abs(diff(latGrid(end-1:end))) 
        warning('SETI_backbone:gagePtRow','The latitutde of the gauge point may be outside the area being modeled.');
    elseif  uniqPairs(ii,2) < lonGrid(1) - 0.5*abs(diff(lonGrid(1:2))) || uniqPairs(ii,2) > lonGrid(end) + 0.5*abs(diff(lonGrid(end-1:end))) 
        warning('SETI_backbone:gagePtCol','The longitude of the gauge point may be outside the area being modeled.');
    end      
        
    nameCurr = ['pt' num2str(round(sub2ind(szGrid, gageRow, gageCol)))];
    
    %Determine the type of data being loaded:
    if regexpbl(nmFile,'stake') || regexpbl(nmFile,{'point','balances'})
        nmData = 'stake';
    elseif regexpbl(nmFile,{'snow','radar'},'and') || regexpbl(nmFile,'swe')
        nmData = 'snowpack';
    else
        error('import_Arc_data:unknownVar',['The variable ' ...
            'corresponding to the data file ' char(39) path char(39) ...
            ' is unknown.  Please program for this variable.'])
    end
    
    %Record data in substructure format:
    sGage.(nameCurr).(nmData) = dataGage(indCurr);
        sGage.(nameCurr).(nmData)(round(sGage.(nameCurr).(nmData)) == -9999) = NaN;
    sGage.(nameCurr).('lon') = lonData(indCurr(1));
    sGage.(nameCurr).('lat') = latData(indCurr(1));
    if iscell(date) && numel(date(:)) == 2
        sGage.(nameCurr).([nmData '_dateStart']) = date{1}(indCurr,:);
        sGage.(nameCurr).([nmData '_dateEnd']) = date{2}(indCurr,:);
    elseif (iscell(date) && numel(date(:)) == 1)
        sGage.(nameCurr).([nmData '_date']) = date{1}(indCurr,:);
    elseif ~iscell(date)
        sGage.(nameCurr).([nmData '_date']) = date(indCurr,:);
    else
        error('import_Arc_data:dateFormat','The date format is unknown');
    end
    
    sGage.(nameCurr).('attributes') = ['Data loaded from ' char(39) path char(39)];
end