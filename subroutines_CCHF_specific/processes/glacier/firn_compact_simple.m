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

function firn_compact_simple()


% if isempty(varargin(:))
%     varargout{1} = cell(0,6);
%     
%     return
% end


%DOES NOT USE:
%Li, J., & Zwally, H. J. (2011). Modeling of firn compaction for estimating 
%ice-sheet mass change from observed ice-sheet elevation change. Annals of 
%Glaciology, 52(59), 1–7.


global sCryo


%Simply sets snow threshold. Any swe greater than that, becomes ice.
threshold = 30; %(meters of water equivalent)


sCryo.frn2ic = zeros(size(sCryo.snw));

indCap = find(sCryo.snw > threshold);

if ~isempty(indCap)
    sCryo.frn2ic(indCap) = sCryo.snw(indCap) - threshold;
    
    sCryo.icdwe(indCap) = sCryo.icdwe(indCap) + sCryo.frn2ic(indCap);
    sCryo.sndwe(indCap) = sCryo.sndwe(indCap) - sCryo.frn2ic(indCap);
    sCryo.snw(indCap) = threshold;
end
