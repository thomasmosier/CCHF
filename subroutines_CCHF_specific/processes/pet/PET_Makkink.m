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


function varargout = PET_Makkink(date, sHydro, varargin)

%See document:
%Makkink, G. F. (1957). Testing the Penman formula by means of lysimeters. 
%J. Inst. Water Eng, 11(3), 277-288.

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
        error('PET_Makkink:dateNumel',['The date vector has ' num2str(numel(date(1,:))) ' elements but 3 are needed.']);
    end
elseif regexpbl(sMeta.dt,'month') %if month time-step, calculate PET for each day and average to monthly value
    error('PET_Hargreaves:monthTs',[ char(39) 'PET_Hargreaves' char(39) ...
        ' can only be used if the time step is daily instead of monthly.']);
else
    error('PET_Makkink:tUnit',['Current time units are ' sMeta.dt ' but days or months are required.  Program for this.'])
end
    


if ~isfield(sLand,'rsdt')
    error('PET_Makkink:toaRaad',[ char(39) 'PET_Hargreaves' char(39) ...
        ' can only be used if top-of-atmosphere radiation is calculated.']);
end



%Find top of the atmosphere radiation:
indRSDT = find(ismember(sLand.datersdt(:,2:end),varargin{1}.currTime(2:end), 'rows') == 1);
if isempty(indRSDT)
    error('PET_Makkink:noPETcurr','No PET grid was calculated for the current time-step');
end

%Slope of the vapor pressure curve (Staghellini Equation):
delta = 0.04145*exp(0.06088*squeeze(sAtm.tas(sAtm.indtas,:,:)));

psych = 0.000665*101.3*((293-0.0065*sHydro.dem)/293).^5.26; %kPa / degC


%Find dt (seconds):
dt = time2sec(1,sMeta.dt,sMeta.currTime); %Used for distributing previous time-steps streamflow (units = seconds)

raadEn = dt*sAtm.rstran.*squeeze(sLand.rsdt(indRSDT,:,:)); %Watts/m^2


sLand.pet = (10^a)*((0.61/2.45) * (delta ./ (delta + psych)).*raadEn - 0.12)/1000; %Units of PET in above are m/day


