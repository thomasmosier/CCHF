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

function varargout = snlq_density(varargin)

%holdCap = decimal between 0 and 1


global sCryo

%WITHOUT HOLDING CAPACITY
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'sn_hold_density', 0,   900,    550, 'snlq_density','cryo'}); %Units of depth melt
    return
else
    holdRho = find_att(varargin{1}.coef,'sn_hold_density','no_warning'); 
end

if ~isfield(sCryo,'rhosn')
   error('snlq_density:noRhoField',['This snow liquid holding '...
       'cpacity function requires the field ' char(39) 'sCryo.rhosn' ...
       char(39) ', which is created in the densification module.']); 
end

%Initialize melt release array:
if ~isfield(sCryo,'snlr')
    sCryo.snlr = zeros(size(sCryo.snw), 'single');
else
    sCryo.snlr(:) = 0; 
end

%%Release liquid in excess of snow holding capacity:
%Find holding capacity (percentage of solid snow):
indRelease = find(sCryo.rhosn > holdRho);

if ~isempty(indRelease)
    %Fraction of water to drain (based on density):
    frac = 
    %Amount of release equals exceedance of liquid water holding capacity:
    sCryo.snlr(indRelease) = frac.*sCryo.lwsnl(indRelease);
    sCryo.sndwe(indRelease) = sCryo.sndwe(indRelease) - sCryo.snlr(indRelease);
    %Remove drained water from snowpack liquid water content:
    sCryo.lwsnl(indRelease) = (1 - frac).*sCryo.lwsnl(indRelease);
end
