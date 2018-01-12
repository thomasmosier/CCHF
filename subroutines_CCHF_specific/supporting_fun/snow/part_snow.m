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

function varargout = part_snow(varargin)
%GETSNOWFALL partitions precip into snow based on the temperature using a
%linear threshold ramp function.
%   INPUTS:
%           ppt     =   timeseries of precip
%           temp    =   timeseries of temperature 
%           ts      =   temperature below which all precip is snow
%           tr      =   temperature above which all precip is rain
%             
%   OUTPUTS:
%           snow    =   the simulated snowfall 
%             
%   Author: Matt Cooper

global sAtm

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    return
else
    sMeta = varargin{1};
end

% if isempty(varargin(:))
% 	argout = cell(2,5);
%     argout(1,:) = {'tSnow', -5, 0, -2, 'part_snow'};
%     argout(2,:) = {'tRain', 0, 5, 3, 'part_snow'};
%     return
% else
%     ts = find_att(varargin{1}.coef,'tSnow'); 
%     tr = find_att(varargin{1}.coef,'tRain');
% end

if regexpbl(sMeta.dt,'month')
   ts = find_att(sMeta.global,'snow_ramp_sn_mon'); %-7.9;
   tr = find_att(sMeta.global,'snow_ramp_rn_mon'); %5.2;
elseif regexpbl(sMeta.dt,{'day','daily','hour'})
   ts = find_att(sMeta.global,'snow_ramp_sn_day'); %-3
   tr = find_att(sMeta.global,'snow_ramp_rn_day'); %1;
end

%Initialize/reset snow field:
sAtm.prsn = zeros(size(squeeze(sAtm.pr(sAtm.indpr,:,:))));

%Cells with entirely snow
indSnow = find(squeeze(sAtm.tas(sAtm.indtas,:,:)) <= ts);
sAtm.prsn(indSnow) = squeeze(sAtm.pr(sAtm.indpr, indSnow ));
%Cells with some snow
indTrans = find(squeeze(sAtm.tas(sAtm.indtas,:,:)) > ts & squeeze(sAtm.tas(sAtm.indtas,:,:)) <= tr);
sAtm.prsn(indTrans) = (tr - squeeze(sAtm.tas(sAtm.indtas, indTrans)))/(tr-ts) ...
    .* squeeze(sAtm.pr(sAtm.indpr, indTrans));
%Remainder:
sAtm.rain = squeeze(sAtm.pr(sAtm.indpr,:,:)) - sAtm.prsn;

