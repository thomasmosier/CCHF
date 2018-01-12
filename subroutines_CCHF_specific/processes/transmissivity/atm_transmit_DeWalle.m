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

function varargout = atm_transmit_DeWalle(varargin)

%Eq. 6.7 on pg. 151 of DeWalle and Rango 2008
 
global sAtm

%VERSION WITHOUT FITTING PARAMETERS:
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1},{'atm_lin_scalar', 0, 30, 16.2, 'atm_transmit_DeWalle','cryo'});
    varargout{1} = cat(1,varargout{1},{'clear_trans', 0.5, 1, 0.65, 'atm_transmit_DeWalle','cryo'});
 
    return
else
    a = find_att(varargin{1}.coef,'atm_lin_scalar'); 
    b = find_att(varargin{1}.coef,'clear_trans'); 
end



%Calculate current temperature range and
if isfield(sAtm,'tasmin') && isfield(sAtm,'tasmax')
    sAtm.tasrng = squeeze(sAtm.tasmax(sAtm.indtasmax,:,:)) - squeeze(sAtm.tasmin(sAtm.indtasmin,:,:)); 
else
    error('atm_transmit_DeWalle:tasMaxMin','tasmax and tasmin are not present, but are needed for present function.');
end
    
    
sAtm.rstran =  b*(1-exp(-0.01*a*sAtm.tasrng.^2.4));

sAtm.rstran(sAtm.rstran > 1) = 1;
sAtm.rstran(sAtm.rstran < 0) = 0; %This line shouldn't be necessary



