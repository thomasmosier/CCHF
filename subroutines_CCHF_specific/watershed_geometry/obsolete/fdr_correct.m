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

function fdr = fdr_correct(fdr)

fdr(        fdr <  1.5)       =   1;
fdr( 1.5 <= fdr & fdr < 3.0 ) =   2;
fdr( 3.0 <= fdr & fdr < 6.0 ) =   4;
fdr( 6.0 <= fdr & fdr < 12.0) =   8;
fdr(12.0 <= fdr & fdr < 24.0) =  16;
fdr(24.0 <= fdr & fdr < 48.0) =  32;
fdr(48.0 <= fdr & fdr < 96.0) =  64;
fdr(96.0 <= fdr       )       = 128;
