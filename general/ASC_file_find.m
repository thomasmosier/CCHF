function [filesUse, varargout] = ASC_file_find(dirTS, varargin)
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


%Designed to work with ASCII data files for which the date is indicated in
%the file as '*year_mnth.asc'



if numel(varargin) == 2
   yrLoad = varargin{1};
   mnthLoad = varargin{2}; 
elseif numel(varargin) == 1
    yrLoad(1) = min(varargin{1}(:,1));
    yrLoad(2) = max(varargin{1}(:,1));
    mnthLoad = varargin{1}(:,2);
else
    error('ASC_file_find:nargin','An incorrect number of input arguments was provided.');
end
mnthsUniq = unique(mnthLoad);
nMnths = numel(mnthsUniq);


if all(~isnan(mnthLoad)) && all(~isnan(yrLoad)) 
    if numel(mnthLoad) ~= numel(yrLoad) && nMnths == numel(mnthLoad)  
        if numel(yrLoad) == 2
            nYrsLd  = yrLoad(2) - yrLoad(1) + 1;
        else
            nYrsLd = numel(unique(yrLoad));
        end

        ldTs = nan(nYrsLd*nMnths,2);

        for ii = 1 : nMnths
            ldTs((ii-1)*nYrsLd+1 : ii*nYrsLd,1) = (yrLoad(1):yrLoad(end))';
            ldTs((ii-1)*nYrsLd+1 : ii*nYrsLd,2) = mnthsUniq(ii);
        end
    elseif numel(yrLoad) == numel(mnthLoad)
        ldTs(:,1) = yrLoad;
        ldTs(:,2) = mnthLoad;
    else
        error('ASC_file_find:dateFormat','The input data format is unknown.')
    end
else
   ldTs = nan(1,2);  
end

%Find root filename and date range:
fileTemp = dir([dirTS filesep '*.asc']);
filesTS = cell(length(fileTemp),1); 
for ii = 1 : length(fileTemp)
    filesTS(ii) = cellstr(fileTemp(ii).name);
end

%Find dates of all files:
indUnd = regexpi(filesTS{1},'_|\.');
nDiv = 0;
for ii = numel(indUnd)-1 : -1 : 1
    if ~isnan(str2double(filesTS{1}(indUnd(ii)+1:indUnd(ii+1)-1)))
        nDiv = nDiv + 1; 
    end
end
fileDate = nan(length(filesTS), nDiv);
for jj = 1:length(filesTS)
    indSep = regexpi(filesTS{jj},'_|\.');
    
    for ii = 1 : nDiv
        indC = numel(indSep) - nDiv + ii - 1;
        fileDate(jj,ii) = str2double(filesTS{jj}(indSep(indC)+1:indSep(indC+1)-1)); 
    end
end


%Find dates requested:
[~, indUse, ~] = intersect(fileDate,ldTs,'rows');

%Find file names of all files to use:
filesUse = cell(numel(indUse),1);
for ii = 1 : numel(indUse)
    filesUse(ii) = fullfile(dirTS, cellstr(filesTS(indUse(ii))));
end
fileDate = fileDate(indUse,:);
   
if numel(indUse) ~= sum(~isnan(ldTs(:,1)))
    warning('asc_file_find:nWrong',[num2str(numel(indUse)) ...
        ' have been loaded from the folder ' char(39) strrep(dirTS,filesep,':') char(39) ...
        ' but the reqested number to load is ' num2str(numel(ldTs(:,1))) '.']);
end

%Put data in order:
[fileDate, indSort] = sortrows(fileDate);
filesUse(indSort);


if numel(varargin) == 1
    [fileDate, ~, indUse] = intersect(varargin{1},fileDate,'rows');
    filesUse = filesUse(indUse);
end

if nargout > 1
   varargout{1} = fileDate; 
   if nargout > 2
       if numel(indUse) > 0
            [~, varargout{2}, metaESRI] = read_ESRI(filesUse{1});
       else
           varargout{2} = nan(6,1);
       end
       
        if nargout > 3
            if numel(indUse) > 0
                varargout{3} = metaESRI;
            else
                varargout{3} = cell(6,1);
            end
        end
   end
end
