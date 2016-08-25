function ind = bnd_wrap(ind,length)
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

%Assumes coordinates on a cylinder and wraps indices around.  Similar to 
%the mathemtical modulus operator but applied to coordinate indexing.

%If longitude on edge of grid being used, must wrap:
indNeg = find(ind < 1);
if ~isempty(indNeg)
    ind(indNeg) = length + ind(indNeg);
end
indOvr = find(ind > length);
if ~isempty(indOvr)
    ind(indOvr) = ind(indOvr) - length;
end