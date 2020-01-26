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

function varargout = icalbedo_constant(varargin)

global sCryo


if isempty(varargin(:))
    varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'albedo_fresh', 0.5, 1, 0.71, 'albedo_Pellicciotti','cryo'});

    return
else
	sMeta = varargin{1};
end

if ~isfield(sCryo, 'icalb')
    %Get albedo of ice
    aIce = find_att(varargin{1}.global,'albedo_ice');

    %Initialize albedo array:
    sCryo.icalb = aIce*ones(size(sCryo.snalb));
        sCryo.icalb(isnan(sCryo.icx)) = nan;

    %Set different albedo for grid cells with debris:
    if isfield(sCryo, 'icdbr')
        aDebris = find_att(sMeta.global,'albedo_debris');
        dDebris = find_att(sMeta.global, 'debris_threshold');

        sCryo.icalb(~isnan(sCryo.icdbr) & sCryo.icdbr > dDebris) = aDebris;
    end
    
    
    %If lake inventory is available, modify glacier albedo values at lake:
    fldLake = '';
    if isfield(sCryo, 'icpndx')
        fldLake = 'icpndx';
    elseif isfield(sCryo, 'iclk')
        fldLake = 'iclk';
    end

    if ~isempty(fldLake)
        aLake = find_att(sMeta.global,'albedo_water');
        
        if any2d((sCryo.(fldLake) < 0 | sCryo.(fldLake) > 1) & ~isnan(sCryo.(fldLake)))
            error('icalbedoConstant:lakeNotFraction', ['The lake inventory ' ...
                'contains values outside the range 0 and 1. It is expected '...
                'that the values should represent fractional area of lake ' ...
                'coverage within each grid cell.']);
        else
            sCryo.icalb = aLake*sCryo.(fldLake) + sCryo.icalb.*(1 - sCryo.(fldLake));
        end
    end
end