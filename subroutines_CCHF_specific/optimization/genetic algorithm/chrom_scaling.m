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

function paramLeng = chrom_scaling(prmSpecs, varargin)

if ~isempty(varargin)
   if numel(varargin{1}) == 1
       nPts = varargin{1}*ones(numel(prmSpecs(:,1)),1);
   elseif isequal(numel(varargin{1}), numel(prmSpecs(:,1)))
       nPts = varargin{1};
   else
       error('chrom_scaling:nDelta','The delta scaling parameter has an unexpected number of entries.')
   end
else
    nPts = 100*ones(numel(prmSpecs(:,1)),1);
end

% if iscell(prmSpecs)
%     prmSpecs = cell2mat(prmSpecs);
% end

paramLeng = ceil(log(nPts)/log(2));
nPts = 2.^paramLeng(:,1) - 1;
    paramLeng( nPts == 1 ) = 0;