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

function varargout = PET_Hammon(date, sHydro, varargin)

%Wi, S., Yang, Y. C. E., Steinschneider, S., Khalil, A., & Brown, C. M. 
%(2015). Calibration approaches for distributed hydrologic models  in 
%poorly gaged basins: implication for streamflow projections under 
%climate change. Hydrol. Earth Syst. Sci., 19(2), 857–876. 
%http://doi.org/10.5194/hess-19-857-2015

%Old version: document 'Hammon_PET_equations' in 'ET' subfolder of SETI model.
%Output PET is in units of m / dt (i.e. m per hour, m per day, or m per month)
global sAtm sLand

if isempty(varargin(:))
	varargout{1} = cell(1,6);
    varargout{1}(1,:) = {'PET_pwr', -3, 2, -0.956, 'PET_Hammon','land'};
    return
else
    sMeta = varargin{1};
    a = find_att(varargin{1}.coef, 'PET_pwr');
end


[~, indDate] = unique(date(:,2:end),'rows');
date = date(indDate,:);

if regexpbl(sMeta.dt,{'day','daily'}) %Calculate for day:
    if numel(date(1,:)) == 3
        dayJ = Julian_day(date,'ofyear');
        PETdate = date;
    else
        error('PET_Hammon:dateNumel',['The date vector has ' num2str(numel(date(1,:))) ' elements but 3 are needed.']);
    end
elseif regexpbl(sMeta.dt,'month') %if month time-step, calculate PET for each day and average to monthly value
    dayJ = nan(31*numel(date(:,1)));
    PETdate = nan([31*numel(date(:,1)),3]);
    cntr = 0;
    for ii = 1 : numel(date(:,1))
        for jj = 1 : eomday(date(ii,1),date(ii,2))
            cntr = cntr + 1;
            PETdate(cntr,:) = [date(ii,1:2),jj];
            dayJ(cntr) = Julian_day(PETdate(cntr,:),'ofyear');
        end
    end
    dayJ = dayJ(1:cntr);
    PETdate = PETdate(1:cntr,:);
else
    error('PET_Hammon:tUnit',['Current time units are ' sMeta.dt ' but days or months are required.  Program for this.'])
end
    
% PET = 0.165*216.7*(daylight hrs / 12 hrs) ...
%     *(6.108*exp((17.27*squeeze(sTemp.data(indTemp,:,:)))/(squeeze(sTemp.data(indTemp,:,:))+237.3))) ...
%     /(squeeze(sTemp.data(indTemp,:,:))+273.3);

if ~isfield(sHydro,'latGrid')
    [~, sHydro.latGrid] = meshgrid(sHydro.lon,sHydro.lat);
end

d = 1 + 0.033*cos(2*pi*dayJ/365);
PET = nan([numel(d),size(sHydro.latGrid)]);
for ii = 1 : numel(d)
    N = daylight_hrs(sHydro, dayJ(ii))/12;

    PET(ii,:,:) = (10^a)*0.001*18.2*N.*exp(17.27*squeeze(sAtm.tas(sAtm.indtas,:,:))./(squeeze(sAtm.tas(sAtm.indtas,:,:))+273.3)) ...
        ./(squeeze(sAtm.tas(sAtm.indtas,:,:))+273.3); %Units of PET in above are m/day
%    PET(ii,:,:) = (10^a)*218.3946*N.*exp((17.27*squeeze(sAtm.tas(sAtm.indtas,:,:)))./(squeeze(sAtm.tas(sAtm.indtas,:,:))+273.3)) ...
%        ./((squeeze(sAtm.tas(sAtm.indtas,:,:))+273.3)*1000); %Units of PET in above are m/day
end

%Convert from m/day to m/(delta time-step)
if regexpbl(sMeta.dt,{'day','daily'})
    if numel(PETdate(1,:)) == 3
        sLand.datepet = PETdate;
        sLand.pet = PET;
    elseif numel(PETdate(1,:)) ~= 3
        error('PET_Hammon:tUnit',['Current time units are ' ...
            sMeta.dt ' but date resolution is ' ...
            num2str(numel(sLand.datepet(1,:))) ', which is not expected.']);
    end
elseif regexpbl(sMeta.dt,'month')
    if numel(PETdate(1,:)) == 3
        sLand.datepet = date;
        sLand.pet = nan([numel(date(:,1)), size(squeeze(sAtm.tas(1,:,:)))],'single');
        for ii = 1 : numel(date(:,1))
            sLand.pet(ii,:,:) = squeeze(sum(PET(ismember(PETdate(:,1:2), date(ii,1:2),'rows'),:,:),1)); 
        end
    elseif numel(PETdate(1,:)) ~= 2
        error('PET_Hammon:tUnit',['Current time units are ' ...
            sMeta.dt ' but date resolution is ' ...
            num2str(numel(sLand.datepet(1,:))) ', which is not expected.']);
    end
end

