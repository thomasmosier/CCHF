% Copyright 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the SETI Hydrologic Modelling Package (referred to 
% as the 'SETI Package').
% 
% The SETI Package is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% The SETI Package is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the SETI Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function varargout = PET_Hargreaves(date, sHydro, varargin)

%See document:
%Hargreaves, G. H., & Samani, Z. A. (1985). Reference Crop 
%Evapotranspiration From Ambient Air Temperature. American Society of 
%Agricultural Engineers.

global sAtm sLand

if isempty(varargin(:))
	varargout{1} = cell(1,6);
    varargout{1}(1,:) = {'PET_pwr', -3, 2, 0, 'PET_Hargreaves','land'};
    return
else
    sMeta = varargin{1};
    a = find_att(varargin{1}.coef, 'PET_pwr');
end


[~, indDate] = unique(date(:,2:end),'rows');
date = date(indDate,:);

if regexpbl(sMeta.dt,{'day','daily'}) %Calculate for day:
    if numel(date(1,:)) == 3
        sLand.datepet = date;
    else
        error('PET_Hargreaves:dateNumel',['The date vector has ' num2str(numel(date(1,:))) ' elements but 3 are needed.']);
    end
elseif regexpbl(sMeta.dt,'month') %if month time-step, calculate PET for each day and average to monthly value
    error('PET_Hargreaves:monthTs',[ char(39) 'PET_Hargreaves' char(39) ...
        ' can only be used if the time step is daily instead of monthly.']);
else
    error('PET_Hargreaves:tUnit',['Current time units are ' sMeta.dt ' but days or months are required.  Program for this.'])
end
    

if ~isfield(sHydro,'latGrid')
    [~, sHydro.latGrid] = meshgrid(sHydro.lon,sHydro.lat);
end

if ~isfield(sLand,'rsdt')
    error('PET_Hargreaves:toaRaad',[ char(39) 'PET_Hargreaves' char(39) ...
        ' can only be used if top-of-atmosphere radiation is used.']);
end
if ~isfield(sAtm, 'tasRange')
    error('PET_Hargreaves:tasRange',[ char(39) 'PET_Hargreaves' char(39) ...
        ' can only be used if minimum and maximum temperature are loaded.']);
end

cConvert = 0.000010861;

%Find top of the atmosphere radiation:
indRSDT = find(ismember(sLand.datersdt(:,2:end),varargin{1}.dateCurr(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('PET_Hargreaves:noPETcurr','No PET grid was calculated for the current time-step');
end

sLand.pet = (10^a)*cConvert*squeeze(sLand.rsdt(indRSDT,:,:)) ...
    .*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) + 17.8).*sAtm.tasRange.^(0.5); %Units of PET in above are m/day


