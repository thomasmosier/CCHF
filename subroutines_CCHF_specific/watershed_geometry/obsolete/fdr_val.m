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

function val = fdr_val(ind)
%indMin = can have values 1 through 9.  Presumes that view is center cell 
%of a 3x3 grid and that labelling 1-9 follows matlab order (i.e. 
%defined by "reshape((1:9),[3,3])")

%key:
%to self = 0
%E  =   1
%SE =   2
%S  =   4
%SW =   8
%W  =  16
%NW =  32
%N  =  64
%NE = 128

%Write output:
switch ind
    case 1
        val = 32; %NW
    case 2
        val = 16; %W
    case 3
        val = 8; %SW
    case 4
        val = 64; %N
    case 5
        val = 0; %Flows to self
    case 6
        val = 4; %S
    case 7
        val = 128; %NE
    case 8
        val = 1; %E
    case 9
        val = 2; %SE
end