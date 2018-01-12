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

function varargout = glacier_slide(varargin)

global sCryo

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    return
else
    sMeta = varargin{1};
end

dt = time2sec(1,sMeta.dt,sMeta.currTime); 


if ~isfield(sCryo,'iceDl')
    rEarth = find_att(sMeta.global,'radius_Earth');
    
    dl = haversine_neighbors(sCryo.iceLon, sCryo.iceLat, rEarth,'sparse');

    %Make field that is distance along flowpath:
    sCryo.iceDl = single(full(reshape(sum(sCryo.iceFdr.*dl, 2), size(sCryo.iceDem))));
    sCryo.iceDl(sCryo.iceDl == 0) = nan; 
end

if ~isfield(sCryo, 'iceFdrM')
   [~, sCryo.iceFdrM, ~] = ...
                wflowacc(sCryo.iceLon, sCryo.iceLat, sCryo.iceDem, ...
                'type','multi', 'edges','open', 'coord','geographic', 'routeflats', 'yes'); 
end

%Set velocity threshold:
velThresh = 0.01/3.154*10^7; %(1 cm per year; m/s)

if any2d(sCryo.iceVel > velThresh)
    frac = dt*full(sCryo.iceVel)./full(sCryo.iceDl); %fractional displacement (unitless)
        frac(isnan(frac)) = 0;
    
    %Add ice from uphill cells:
    dispIce = frac(:).*full(sCryo.ice(:)); 
    sCryo.iceWE(:) = sCryo.iceWE(:) + sCryo.iceFdrM'*dispIce; %Units of meters
    %Remove ice from uphill cell:
    sCryo.iceWE = 1-dispIce; %Units of meters

    %Add displaced ice to ice change field:
    sCryo.iceDWE(:) = sCryo.iceDWE(:) + sCryo.iceFdrM'*dispIce(:) - dispIce(:);
end