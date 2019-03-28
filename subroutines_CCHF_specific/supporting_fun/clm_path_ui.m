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

function sPath = clm_path_ui(sPath, sMeta, varargin)

if ~isempty(varargin(:))
    region = varargin{1};
else
    if ischar(sMeta.region)
        region = sMeta.region;
    elseif iscell(sMeta.region) && numel(sMeta.region(:)) == 1
        region = sMeta.region{1};
    else
        error('clm_path_ui:unknownRegion','The function was not able to detect the region.');
    end
end

if isfield(sPath, 'dem') %Use dem path (if available)
    [startPath, ~, ~] = fileparts(sPath.dem);
else
 	nmsTemp = fieldnames(sPath);
    if ~isempty(nmsTemp) %Pick first available field
        [startPath, ~, ~] = fileparts(sPath.(nmsTemp{1}));
    else
        startPath = pwd;
    end
end
    



for ii = 1 : numel(sMeta.varLd)
    uiwait(msgbox(sprintf(['Select the folder containing ' sMeta.varLdDisp{ii} ...
        ' time-series data for ' region '.\n']), ...
        '(Click OK to Proceed)','modal'));
    sPath.(sMeta.varLd{ii}) = uigetdir(startPath,['Select the folder containing ' sMeta.varLdDisp{ii} ...
        ' time-series data for ' region]);
    disp([char(39) sPath.(sMeta.varLd{ii}) char(39) ' has been chosen as the ' sMeta.varLdDisp{ii} ' time-series data.']);
    
    %Update search path:
    [startPath, ~, ~] = fileparts(sPath.(sMeta.varLd{ii}));
%     indStart = regexpi(startPath,filesep);
%     startPath = startPath(1:indStart(end)-1);
end
