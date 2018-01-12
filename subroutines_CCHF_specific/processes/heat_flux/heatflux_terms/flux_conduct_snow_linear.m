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

function argout = flux_conduct_snow_linear(tsis, varargin)

%Positive values indicate heat flowing from internal snowpack to snow-air interface
%Negative values indicate heat flowing from snow-air interface to internal snowpack


global sCryo

if isempty(varargin(:))
	argout = cell(0,6);
    argout = cat(1, argout, {'cond_pwr', -2, 2, 1, 'heat_conduction_snow','cryo'});  
    argout = cat(1,argout, ['tsn',cell(1,5)]);
    return
else
    c = find_att(varargin{1}.coef,'cond_pwr');
    
    if numel(varargin(:)) > 1
       mode = varargin{2};
    else
        mode = '';
    end
end

%Set maximum value for conduction:
condtMax = 500; %Units = Watts per meter squared


if ~isfield(sCryo,'rhosn')
    sCryo.rhosn = 400*ones(size(sCryo.snw),'single');
end

if ~regexpbl(mode,'deriv') %Normal mode:
    if isfield(sCryo, 'tsn')
        argout = -10^(c).*(tsis - sCryo.tsn) ./ sCryo.snw;
    else
        error('heat_conduct_snow:missingTmp','Surface or internal snow temperature not defined.');
    end
else %Derivative of above:
    argout = -10^(c) ./ sCryo.snw;
end

%Ensure no odd values (due to there not being any snow:
argout(isnan(argout)) = 0;
argout(sCryo.snw == 0) = 0;
argout(argout > condtMax) = condtMax;
argout(argout < -condtMax) = -condtMax;

%If temperatures equal, then no conduction
argout(sCryo.tsn == tsis) = 0;