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

function varargout = toa_rad_DeWalle(date, sHydro, varargin)

%Formulation comes from
%DeWalle, D. R., & Rango, A. (2008). Principles of snow hydrology. 
%Cambridge University Press. (Appendix B)

global sAtm
%OUTPUT SCALING FITTING PARAMETER:
if isempty(varargin(:))
	varargout{1} = cell(0,6);
 
    return
else
    sMeta = varargin{1};
end


%NO FITTING PARAMETERS:
% if isempty(varargin(:))
% 	argout = cell(0,5);
%     return
% else
%     sMeta = varargin{1};  
% end

%Solar Constant (W/m^2):
I0 =  find_att(varargin{1}.global, 'solar_constant'); 

[~, indUse, ~] = unique(date(:,2:end), 'rows');
dateUse = date(indUse,:);

sAtm.rsdt = solar_daily_Dewalle(sHydro.latGrid, sHydro.aspectFdr, sHydro.slopeFdrA, 'slope', dateUse, sMeta.dt, 'solar_constant', I0);
sAtm.datersdt = dateUse;