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

function argout = flux_long_Hock(varargin)
%R_long from atm is parameterized using Deacon (1970)

global sCryo sAtm

%VERSION WITH ONLY ONE FITTING PARAMETER
if isempty(varargin(:))
	argout = cell(0,6);
    argout = cat(1,argout, {'lw_pwr_tsn', -2, 2, 0, 'flux_long_Hock', 'cryo'});
    argout = cat(1,argout, ['tsn',cell(1,5)]);
    return
else
    scaleOut = find_att(varargin{1}.coef,'lw_pwr'); 
end

eSky = 0.7; %Emmissivity of sky
% stefan = 5.67*10^(-8);

%Primary reference is:
%Hock, R. (2005). Glacier melt: a review of processes and their modelling. 
%Progress in Physical Geography, 29(3), 362–391.

argout = (10^scaleOut)*5.67*10^(-8)*(eSky*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) + 273.15).^4 ...
    - (sCryo.tsn + 273.15).^4);