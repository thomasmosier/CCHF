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

function chrom = chrom_dec2bin(parameters, prmBnds, paramLeng)

chrom = cell(length(parameters(:,1)),1);
[ chrom{:} ]= deal( blanks( sum(paramLeng)));

if iscell(prmBnds)
    prmBnds = cell2mat(prmBnds);
end

nPts = 2.^paramLeng - 1;

nBit = 1;
for ii = 1 : numel(paramLeng)  %Loop over each parameter
    
    if prmBnds(ii,2) - prmBnds(ii,1) == 0
        continue;
    else
        
        for jj = 1 : length( parameters(:,1) )  %Loop over each member of parameter population
            chrom{jj}(nBit:nBit + paramLeng(ii) - 1) = dec2bin((parameters(jj,ii) - prmBnds(ii,1))*nPts(ii)/(prmBnds(ii,2) - prmBnds(ii,1)),paramLeng(ii));
        end
        nBit = nBit + paramLeng(ii);
    end
end