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

function [dirOutput, strModules] = CCHF_out_dir(sPath, sMeta)



%Find Name of directory to use:
if isfield(sPath, 'resume') && ~isempty(sPath.resume)
    dirTS = sPath.resume;
elseif ~isempty(sPath.dem) %If DEM path included, use this
    [dirTS,~,~] = fileparts(char(sPath.dem));
    [dirTS,~, ~] = fileparts(dirTS);
elseif ~isempty(sPath.pre) %Elseif, use precipitation path
   [dirTS,~,~] = fileparts(char(sPath.pre)); 
   [dirTS,~, ~] = fileparts(dirTS);
else
    error('downscale:noData','It appears no input data were selected.');
end

%Get module string from sMeta, then edit and write
strModules = blanks(0);
for ii = 1 : numel(sMeta.module(:,1))
    if ii ~= 1
        strModules = [strModules '-'];
    end
    
    strModules = [strModules sMeta.module{ii,1} '_' sMeta.module{ii,2}];
end


%Remove modules from string if not used in current version (e.g. in simple degree index):
if regexpbl(strModules, 'heat-simple') || regexpbl(strModules, 'energy-simple')
	indTrans = regexpi(strModules,'trans');
	if ~isempty(indTrans)
        indUnd = regexpi(strModules,'_');
		indRem = find(indUnd > indTrans, 1, 'first');
		strModules = [strModules(1:indTrans-1) strModules(indUnd(indRem)+1:end)];
	end
	indAlbedo = regexpi(strModules,'albe');
	if ~isempty(indAlbedo)
        indUnd = regexpi(strModules,'_');
        indRem = find(indUnd > indAlbedo, 1, 'first');
        strModules = [strModules(1:indAlbedo-1) strModules(indUnd(indRem)+1:end)];
	end
	indToa = regexpi(strModules,'toa');
	if ~isempty(indToa)
        indUnd = regexpi(strModules,'_');
		indRem = find(indUnd > indToa, 1, 'first');
		strModules = [strModules(1:indToa-1) strModules(indUnd(indRem)+1:end)];
	end
end

%Only keep names of first four modules:
indUnd = regexpi(strModules,'_');
strModules = strModules(1:indUnd(min(4,numel(indUnd)))-1);

%Name of output subdirectory:
dirCCHF = [sMeta.runType '_' ...
    num2str(sMeta.dateStart(1)) 'thru' num2str(sMeta.dateEnd(1)) ...
    '-' strModules];


%Ensure output directory unique (doesn't write over existing data):
%Exception: if 'resume' is being used...
if isfield(sPath, 'resume') && strcmpi(dirTS, sPath.resume)
    dirOutput = dirTS;
    
    %Display name of output path:
	disp(['The output path is ' char(39) dirOutput char(39) '.' char(10) ...
	'This existing path is being used because CCHF is resuming a previous run.']);
else
    cntrDir = 0;
    indSep = regexpi(dirTS,filesep);
    while exist(fullfile(dirTS, dirCCHF),'dir')
        cntrDir = cntrDir + 1;
        if cntrDir == 1
            dirCCHF = [dirCCHF,'_' num2str(cntrDir)];
        else
            indC = regexpi(dirCCHF,'_');
            dirCCHF = [dirCCHF(1:indC(end)), num2str(cntrDir)];
        end
    end
    
    dirOutput = fullfile(dirTS, dirCCHF);
    %Create output directory
    if ~exist(dirOutput, 'dir')
        mkdir(dirOutput);
    end
    
    %Display name of output path:
	disp(['The output path is ' char(39) dirOutput char(39) '.' char(10)]);
end


