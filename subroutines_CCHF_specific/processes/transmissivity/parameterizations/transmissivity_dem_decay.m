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

function trans = transmissivity_dem_decay(tasRng, elev, a, b, c, d)

%This is a simple exponential decay function with a similar, but simplified
%form relative to Coops et al. (2000).

%Coops, N. C., Waring, R. H., & Moncrieff, J. B. (2000). Estimating mean 
%monthly incident solar radiation on horizontal and inclined slopes from 
%mean monthly temperatures extremes. International Journal of 
%Biometeorology, 44(4), 204–211.

 
tClear = a+b*elev/1000;

%If elevation is 2d and tasrng is 3d, make 3d version of tClear
if ismatrix(tClear) && ndims(tasRng) == 3 && isequal(size(tClear), size(squeeze(tasRng(1,:,:))))
    temp = nan([1,size(tClear)], 'single');
    temp(1,:,:) = tClear;
    tClear = repmat(temp, numel(tasRng(:,1,1)), 1, 1);
elseif ndims(tasRng) == 3 && ndims(tClear) == 3 && ~isequal(size(tClear), size(tasRng))
    error('transmissivityCoops:arraySzDiff', 'Both input arrays are 3d, but their sizes are different.');
end

trans = tClear.*(1 - exp(-c*0.01*tasRng.^d));

trans(trans > 1) = 1;
trans(trans < 0) = 0;

