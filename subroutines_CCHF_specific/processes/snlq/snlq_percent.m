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

function varargout = snlq_percent(varargin)

%holdCap = decimal between 0 and 1


global sCryo

%WITHOUT HOLDING CAPACITY
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'sn_hold', 0,   .1,    0.077, 'snlq_percent','cryo'}); %Units of depth melt
    return
else
    holdCap = find_att(varargin{1}.coef,'sn_hold','no_warning'); 
%     if isempty(holdCap)
%        holdCap = 0; 
%     end
%     a = 0.03;  
end


%Initialize melt release array:
sCryo.snlr = zeros(size(sCryo.snw),'single');

%%Release liquid in excess of snow holding capacity:
%Find holding capacity (percentage of solid snow):
sCryo.snlh = holdCap*sCryo.snw;
%Find indices where it's exceeded:
indRelease = find(sCryo.lwsnl > sCryo.snlh);
if ~isempty(indRelease)
    %Amount of release equals exceedance of liquid water holding capacity:
    sCryo.snlr(indRelease) = sCryo.lwsnl(indRelease) - sCryo.snlh(indRelease);
    sCryo.sndwe(indRelease) = sCryo.sndwe(indRelease) - sCryo.snlr(indRelease);
    %Remove drained water from snowpack liquid water content:
    sCryo.lwsnl(indRelease) = sCryo.snlh(indRelease);
end
