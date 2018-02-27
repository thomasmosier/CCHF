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

function argout = flux_long_c_emm(tsis, varargin)
%Assumes emissivity from clear-sky = 0.7 and with rain = 1
%References:
%Hock, R. (2005). Glacier melt: a review of processes and their modelling. 
%Progress in Physical Geography, 29(3), 362–391.

%Plüss, C., & Ohmura, A. (1997). Longwave radiation on snow-covered 
%mountainous surfaces. Journal of Applied Meteorology, 36(6), 818–824.

global sAtm

%VERSION WITH ONLY ONE FITTING PARAMETER
if isempty(varargin(:))
	argout = cell(0,6);
%     argout{1} = cat(1,argout{1}, {'lw_pwr', -2, 2, 0, 'heat_long_simple', 'cryo'});

    argout = cat(1,argout, ['tsn', cell(1,5)]); %This line signals that snow temperature is needed in process representation.
    return
% else
%     scaleOut = find_att(varargin{1}.coef,'lw_pwr_tsn'); 
else
    sMeta = varargin{1};
    
    if numel(varargin(:)) > 1
       mode = varargin{2};
    else
        mode = '';
    end
end


%Find emmissivity values
eSkyClear = find_att(sMeta.global, 'emm_clear');
 

if ~isfield(sAtm,'emm')
   sAtm.emm = eSkyClear*ones(size(sAtm.rain)); 
end    

% stefan = 5.67*10^(-8);

if ~regexpbl(mode,'deriv') %Normal evaluation:
    argout = 5.67*10^(-8)*(sAtm.emm.*(squeeze(sAtm.tas(sAtm.indtas,:,:)) + 273.15).^4 ...
        - (tsis + 273.15).^4);
else %Derivative of above:
    argout = -4*5.67*10^(-8)*(tsis + 273.15).^3;
end

argout(isnan(argout)) = 0;