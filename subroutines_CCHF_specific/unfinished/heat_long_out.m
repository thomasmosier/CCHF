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

function varargout = heat_long_out(varargin)
%R_long from atm is parameterized using Deacon (1970)

global sCryo sAtm

%VERSION WITH ONLY ONE FITTING PARAMETER
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {    'lw_out_pwr', -2,  2, -0.8, 'heat_long_out','cryo'});
%     varargout{1}(2,:) = {'emm_sky', 0, 1, 0.8, 'heat_long_simple', 'cryo'};
    return
else
    lwPwr = find_att(varargin{1}.coef,'lw_out_pwr');
%     emm = find_att(varargin{1}.coef,'emm_sky'); 
end

% stefan = 5.67*10^(-8);
if ~isfield(sCryo,'sntas')
    sCryo.sntas = single(squeeze(sAtm.tas));
    sCryo.sntas(sCryo.sntas > 0 == 0);
end


sCryo.heatLR = 10^(lwPwr)*single(-5.67*10^(-8)*(sCryo.sntas + 273.15.^4));
sCryo.heatLR(sCryo.solid == 0 | sCryo.ice == 0 | isnan(sCryo.ice)) = 0;

sCryo.sntas = squeeze(sAtm.tas);
sCryo.sntas(sCryo.sntas > 0 == 0);
% 
% sCryo.heatLR = (10^scaleOut)*5.67*10^(-8)*(eSky*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) + 273.15).^4 ...
%     - (sCryo.tasp + 273.15).^4);
