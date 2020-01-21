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

function varargout = snalbedo_Brock(varargin)

%Primary reference:
%Brock, B. W., Willis, I. C., & Sharp, M. J. (2000). Measurement and 
%parameterization of albedo variations at Haut Glacier d’Arolla, 
%Switzerland. Journal of Glaciology, 46(155), 675–688.

%See albedo discussion in:
%Pellicciotti, F., Brock, B., Strasser, U., Burlando, P., Funk, M., & 
%Corripio, J. (2005). An enhanced temperature-index glacier melt model 
%including the shortwave radiation balance: development and testing for 
%Haut Glacier d'Arolla, Switzerland. Journal of Glaciology, 51(175), 
%573-587.

global sCryo sAtm



if isempty(varargin(:))
    varargout{1} = snow_albedo_Brock(0.1, 0.8, 0.5); %Doesn't matter because this is only for parameter mode
    return
end

%NO FITTING PARAMETERS
aUnder = find_att(varargin{1}.global,'albedo_ground'); %Underlying or debris albedo (Unitless)

aFresh = find_att(varargin{1}.global,'albedo_snow_fresh'); %Albedo of fresh snow
aOld   = find_att(varargin{1}.global,'albedo_snow_old');  %Albedo of old/melting snow

%FITTING PARAMETERS:
% if isempty(varargin(:))
%     argout = snow_albedo_Brock(sSnowpack, aUnder); 
%     argout(end+1,:) = {'albedo_ice', 0.2, 0.7, 0.45, 'albedo_Brock'}; 
%     return
% else  
%     aIce = find_att(varargin{1}.coef,'albedo_ice');
% end

if ~isfield(sAtm,'tasmaxc')
    sAtm.tasmaxc = zeros(size(sCryo.snw),'single');
end

%Update accumulated daily maximum air temperature since last snow event (used in albedo calculation):
sAtm.tasmaxc = sAtm.tasmaxc + squeeze(sAtm.tasmax(sAtm.indtasmax,:,:));
sAtm.tasmaxc( sAtm.prsn > 0 ) = 0;
sAtm.tasmaxc(sAtm.tasmaxc < 0) = 0;

snow_albedo_Brock(aUnder, aFresh, aOld, varargin{1});

sCryo.snalb(isnan(sCryo.icx)) = nan;