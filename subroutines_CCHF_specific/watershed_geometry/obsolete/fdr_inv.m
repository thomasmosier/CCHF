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

function out = fdr_inv(val,varargin)
%Val = values output by fdr_val
%Indice = vector indice of 3x3 grid (i.e. indice if 3x3 grid converted to
%reshape((1:9),[3 3])


%key:
    %Self = 0;
    %E  =   1
    %SE =   2
    %S  =   4
    %SW =   8
    %W  =  16
    %NW =  32
    %N  =  64
    %NE = 128
    
if isempty(varargin(:))
    switch val
        case 1
            out = 8;
        case 2
            out = 9;
        case 4
            out = 6;
        case 8
            out = 3;
        case 16
            out = 2;
        case 32
            out = 1;
        case 64
            out = 4;
        case 128
            out = 7;
        case 0
            out = 5;
    end
else
    
    switch val
        case 1
            out = [ 0,  1];
        case 2
            out = [ 1,  1];
        case 4
            out = [ 1,  0];
        case 8
            out = [ 1, -1];
        case 16
            out = [ 0, -1];
        case 32
            out = [-1, -1];
        case 64
            out = [-1,  0];
        case 128
            out = [-1,  1];
        case 0
            out = [ 0,  0];
    end
end