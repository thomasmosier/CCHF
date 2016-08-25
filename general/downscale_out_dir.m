function dirOutput = downscale_out_dir(sPath,sMeta,sDownscale)
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


%Find current index:
indUse = curr_ind(sMeta);

%Find Name of time series directory:
if ~isempty(sPath.projMod{indUse}) %GCM (downscaling projected)
    [dirTS,~,~] = fileparts(char(sPath.projMod{indUse}));

    %Find name of projection GCM used:
    [~, nmFutGcm,~] = fileparts(sPath.projMod{indUse});
    indFutGcmUnd = regexpi(nmFutGcm,'_');
    if ~isempty(indFutGcmUnd) && length(indFutGcmUnd) > 1
        nmFutGcm = nmFutGcm(indFutGcmUnd(2)+1:indFutGcmUnd(end)-1);
    end

%     if ~isempty(sPath.hisTS{indUse})
%         %Find name of historical GCM used:
%         [~, nmHisGcm,~] = fileparts(sPath.futGCM{indUse});
%         indHisGcmUnd = regexpi(nmHisGcm,'_');
%         if ~isempty(indFutGcmUnd) && length(indFutGcmUnd) > 1
%             nmHisGcm = nmHisGcm(indHisGcmUnd(2)+1:indHisGcmUnd(end)-1);
%         end
%     end
elseif isempty(sPath.projMod{indUse}) && ~isempty(sPath.hisMod{indUse}) %gridded historical (downscaling hindcast)
   dirTS = char(sPath.hisMod{indUse});    
else
    error('downscale:noTSselected','It appears no gridded time-series was selected.');
end


%Name of output subdirectory:
dirDownscale = [sDownscale.method '_' sDownscale.intrp '_' sMeta.currVar '_' sMeta.region];
if ~isempty(sPath.projMod{indUse})
    dirDownscale = [dirDownscale '_' nmFutGcm];
end

%Ensure output directory unique (doesn't write over existing data):
cntrDir = 0;
indSep = regexpi(dirTS,filesep);
while exist(fullfile(dirTS(1:indSep(end)-1), dirDownscale),'dir')
    cntrDir = cntrDir + 1;
    if cntrDir == 1
        dirDownscale = [dirDownscale,'_' num2str(cntrDir)];
    else
        indC = regexpi(dirDownscale,'_');
        dirDownscale = [dirDownscale(1:indC(end)), num2str(cntrDir)];
    end
end
%Display name of output path:
disp(['The output path is ' char(39) ...
    fullfile(dirTS(1:indSep(end)-1), dirDownscale) char(39) '.' char(10)]);

%Create output directory
dirOutput = mk_output_dir(dirTS, dirDownscale);