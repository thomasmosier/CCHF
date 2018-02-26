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

function varargout = atmtrans_dem_decay(sHydro,varargin)

%This is a simple exponential decay function with a similar, but simplified
%form relative to Coops et al. (2000).

%Coops, N. C., Waring, R. H., & Moncrieff, J. B. (2000). Estimating mean 
%monthly incident solar radiation on horizontal and inclined slopes from 
%mean monthly temperatures extremes. International Journal of 
%Biometeorology, 44(4), 204–211.

 
global sAtm

%VERSION WITHOUT FITTING PARAMETERS:
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'trans_clear_intercept',   0, 1, 0.71, 'atm_tansmit_dem_decay','cryo'});
    varargout{1} = cat(1,varargout{1}, {'trans_clear_dem',   0, 0.08, 0.023, 'atm_tansmit_dem_decay','cryo'});
    varargout{1} = cat(1,varargout{1}, {'trans_decay_rate',   1, 7, 3.4, 'atm_tansmit_dem_decay','cryo'});
    varargout{1} = cat(1,varargout{1}, {'trans_decay_power',   1, 4, 1.6, 'atm_tansmit_dem_decay','cryo'});
    return
else
    a = find_att(varargin{1}.coef,'trans_clear_intercept'); 
    b = find_att(varargin{1}.coef,'trans_clear_dem'); 
    c = find_att(varargin{1}.coef,'trans_decay_rate'); 
    d = find_att(varargin{1}.coef,'trans_decay_power'); 
end


%Calculate current temperature range and
if isfield(sAtm,'tasmin') && isfield(sAtm,'tasmax')
    sAtm.tasrng = squeeze(sAtm.tasmax(sAtm.indtasmax,:,:)) - squeeze(sAtm.tasmin(sAtm.indtasmin,:,:)); 
else
    error('atm_transmit_Coops:tasMaxMin','tasmax and tasmin are not present, but are needed for present function.');
end

sAtm.rstran = transmissivity_dem_decay(sAtm.tasrng, sHydro.dem, a, b, c, d);

