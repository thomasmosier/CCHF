function out = ceil2(num,dec,varargin)
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

%num = number to be rounded
%dec = precision of rounding
%varargin = factor*10^dec to round to (e.g. if num = -2700, dec = -2, and 
    %varargin = 2, then out = -2800) 

if ~isempty(varargin) && isnumeric(varargin{1})
   fact = varargin{1};
else
    fact = 1;
end
    
out = ceil(10.^(dec)*num./fact).*fact./10.^(dec);