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

function slope = slope_Liston(dem,fdr,lon,lat)
%Slope is the angle of inclination along the flowpath


slope = nan(size(dem));

if sum(size(lon) == 1) > 0 && sum(size(lon) > 1) > 0
    [lon, lat] = meshgrid(lon, lat);
end

tic
for jj = 2 : numel(dem(:,1)) - 1
   
%    disp([num2str(jj) ', ' num2str(toc)]);
   
   for ii = 2 : numel(dem(1,:)) - 1
       crd = fdr_inv_mult(fdr(jj,ii),2);
       
        if numel(crd) == 2
            slope(jj,ii) = atand(slope_2pts(dem(jj,ii) - dem(jj+crd(1),ii+crd(2)), [lon(jj,ii), lat(jj,ii)], [lon(jj+crd(1),ii+crd(2)), lat(jj+crd(1),ii+crd(2))]));
        elseif numel(crd) == 4
            if sum(abs(crd(1,:) - crd(2,:)) > 1)
                slope(jj,ii) = 0;
            else
                slope(jj,ii) = mean([ ...
                    atand(slope_2pts(dem(jj,ii) - dem(jj+crd(1,1),ii+crd(1,2)), [lon(jj,ii), lat(jj,ii)], [lon(jj+crd(1,1),ii+crd(1,2)), lat(jj+crd(1,1),ii+crd(1,2))])), ...
                    atand(slope_2pts(dem(jj,ii) - dem(jj+crd(2,1),ii+crd(2,2)), [lon(jj,ii), lat(jj,ii)], [lon(jj+crd(2,1),ii+crd(2,2)), lat(jj+crd(2,1),ii+crd(2,2))])) ...
                    ]);
            end
        else
            slope(jj,ii) = 0;
        end
           
   end
end