function [filesUse, varargout] = GPCC_file_find(dirTS, yrLoad, mnthLoad)
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



mnthsUniq = unique(mnthLoad);
nMnths = numel(mnthsUniq);

if numel(mnthLoad) ~= numel(yrLoad) && nMnths == numel(mnthLoad)  
    if numel(yrLoad) == 2
        nYrsLd  = yrLoad(2) - yrLoad(1) + 1;
    else
        nYrsLd = numel(unique(yrLoad));
    end
    
    ldTs = nan(nYrsLd*nMnths,2);
        
    for ii = 1 : nMnths
        ldTs((ii-1)*nYrsLd+ii : ii*nYrsLd,1) = (yrLoad(1):yrLoad(end))';
        ldTs((ii-1)*nYrsLd+ii : ii*nYrsLd,2) = mnthsUniq(ii);
    end
elseif numel(yrLoad) == numel(mnthLoad)
    ldTs(:,1) = yrLoad;
    ldTs(:,2) = mnthLoad;
else
    error('ASC_file_find:dateFormat','The input data format is unknown.')
end

%Find all GPCC time-series files in directory:
filesTS = dir([dirTS filesep 'gpcc_*_degree_*']);

indSep = regexpi(filesTS(1).name,'_');

%Find time indices of all files:
fileDate = nan(numel(filesTS(:,1)),2);
for jj = 1:length(filesTS)
    fileDate(jj,1) = str2double(filesTS(jj).name(indSep(end)+3:indSep(end)+6)); 
    fileDate(jj,2) = str2double(filesTS(jj).name(indSep(end)+1:indSep(end)+2)); 
end

%Find dates requested:
[~, indUse, ~] = intersect(fileDate,ldTs,'rows');

%Find file names of all files to use:
filesUse = cell(numel(indUse),1);
for ii = 1 : numel(indUse)
    filesUse(ii) = cellstr(filesTS(indUse(ii)).name);
end
      
if numel(indUse) ~= sum(~isnan(ldTs(:,1)))
    error('GPCC_file_find:nWrong',[num2str(numel(indUse)) ...
        ' have been loaded from the folder ' char(39) dirTS char(39) ...
        ' but the reqested number to load is ' num2str(numel(ldTs(:,1))) '.']);
end

%Extra output arguments:
if nargout > 1
   varargout{1} = ldTs; 
   if nargout > 2
       if numel(indUse) > 0
            [~, varargout{2}, metaESRI, ~, ~] = read_GPCC_ASCII(filesUse{1});
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