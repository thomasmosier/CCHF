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

function varargout = scx_max(varargin)
%Assumes snow is uniformly distributed over grid cell


global sCryo

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    
    return
else
    sMeta = varargin{1};
end

optThresh = find_att(sMeta.global,'SCA_SWE_opt_depth'); 

%reset each time step
sCryo.scx = zeros(size(sCryo.snw),'single');

sCryo.scx(sCryo.snw >= optThresh) = 1;