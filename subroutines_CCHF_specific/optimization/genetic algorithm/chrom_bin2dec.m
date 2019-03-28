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


function parameters = chrom_bin2dec(chromosomes, prmBnds, paramLeng)

if iscell(prmBnds)
    prmBnds = cell2mat(prmBnds);
end

%1st indice is parameter set (i.e. individual), 2nd indice is parameter
parameters = zeros( numel(chromosomes(:,1)), length(prmBnds(:,1)));

nPts = 2.^paramLeng - 1;

ind = 1;
for ii = 1 : length(prmBnds(:,1))
    if isequal(prmBnds(ii,1), prmBnds(ii,2))
       parameters(:,ii) = prmBnds(ii,1);
    else
        for jj = 1 : numel(chromosomes(:))
            parameters(jj,ii) = prmBnds(ii,1) ...
                + (prmBnds(ii,2) - prmBnds(ii,1))/nPts(ii)*bin2dec( chromosomes{jj}(ind : ind + paramLeng(ii)-1) ); 
        end
        ind = ind + paramLeng(ii); 
    end
    
end

