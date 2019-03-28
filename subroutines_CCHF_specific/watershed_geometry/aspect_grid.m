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

function aspect = aspect_grid(fdr)
%'aspect' is the slope's azimuth (degrees clockwise from North)



aspect = nan(size(fdr));

for jj = 2 : numel(fdr(:,1)) - 1
	for ii = 2 : numel(fdr(1,:)) - 1     
        ind = fdr_inv_mult(fdr(jj,ii));
        
        out = nan(numel(ind),1);

        for kk = 1 : numel(ind)
            switch ind(kk)
            %E  =   1;      8
            %SE =   2;      9
            %S  =   4;      6
            %SW =   8;      3
            %W  =  16;      2
            %NW =  32;      1
            %N  =  64;      4
            %NE = 128;      7
                case 8
                    out = 90;
                case 9
                    out = 135;
                case 6
                    out = 180;
                case 3
                    out = 225;
                case 2
                    out = 270;
                case 1
                    out = 315;
                case 4
                    out = 0;
                case 7
                    out = 45;
                case 0
                    out = 5;
            end
        end

        aspect(jj,ii) = nanmean(out);
	end
end