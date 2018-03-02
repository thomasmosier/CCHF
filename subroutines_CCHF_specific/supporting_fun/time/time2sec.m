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

function dt = time2sec(n,unit,varargin)
%Calculate number of seconds in given unit of time.

if regexpbl(unit,{'day','daily'})
    dt = 86400;
elseif regexpbl(unit,'month')
    if ~isempty(varargin(:))
        dt = 86400*eomday(varargin{1}(1),varargin{1}(2));
    else
        error('time2sec:curr_Date','The current year and month must be supplied since monthly time step used.');
    end
elseif regexpbl(unit,'hour')
    dt = 3600;
else
    error('time2sec:dt',[char(39) unit char(39) ' is an unknown time step.']);
end

dt = n*dt;