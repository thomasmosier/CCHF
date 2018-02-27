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

function varargout = partition_ramp(varargin)
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
%     varargout{1} = cat(1,varargout{1}, {'albedo_fresh', 0.5, 1, 0.71, 'albedo_Pellicciotti','cryo'});
%     varargout{1} = cat(1,varargout{1}, {'albedo_decay', 0.0, 0.5, 0.11, 'albedo_Pellicciotti','cryo'});
%     varargout{1} = cat(1,varargout{1}, {  'albedo_ice', 0.2, 0.7, 0.45, 'albedo_Pellicciotti','cryo'});
    
    return
else
%     aFresh = find_att(varargin{1}.coef,'albedo_fresh'); %Albedo of fresh, deep snow
%     aDecay = find_att(varargin{1}.coef,'albedo_decay');  
%     aIce = find_att(varargin{1}.coef,'albedo_ice');
    sMeta = varargin{1};
end

if regexpbl(sMeta.dt,'month')
    %Values from comparison to global climate model snow fields
   ts = -7.9;
   tr = 5.2;
elseif regexpbl(sMeta.dt,{'day','daily','hour'})
    %Values from visual inspection of figures in
    %Dai, A. (2008). Temperature and pressure dependence of the rain-snow 
    %phase transition over land and ocean. Geophysical Research Letters, 
    %35(12), n/a-n/a. https://doi.org/10.1029/2008GL033295
   ts = -1; %-1
   tr = 3; %3;
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

