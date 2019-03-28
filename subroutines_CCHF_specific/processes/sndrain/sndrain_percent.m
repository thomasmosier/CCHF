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

function varargout = sndrain_percent(varargin)

%holdCap = decimal between 0 and 1


global sCryo

%WITHOUT HOLDING CAPACITY
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'sn_hold', 0,   .1,    0.077, 'sndrain_percent', 'cryo'}); %Units of depth melt
    %Calibrate with 'land' routines because it is more sensitive to runoff
    %generation and flow than cryosphere
    return
else
    holdCap = find_att(varargin{1}.coef,'sn_hold');  
end


%Initialize snowmelt release array:
sCryo.snlr = zeros(size(sCryo.snw),'single');


%Find holding capacity (percentage of solid snow):
sCryo.snlh = holdCap*sCryo.snw;


%Release snowpack liquid in excess of snowpack holding capacity:
indRelease = find(sCryo.snlw > sCryo.snlh);
if ~isempty(indRelease)
    %Amount of release equals exceedance of liquid water holding capacity:
    sCryo.snlr(indRelease) = sCryo.snlw(indRelease) - sCryo.snlh(indRelease);

    %Remove drained water from snowpack liquid water content:
    sCryo.snlw(indRelease) = sCryo.snlh(indRelease);
end
