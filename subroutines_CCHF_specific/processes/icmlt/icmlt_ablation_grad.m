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

function varargout = icmlt_ablation_grad(sHydro, varargin)
%%MELT ICE:
    
global sCryo sAtm


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'ice_tas_offset' , -5,   5,  0, 'icmlt_ablation_grad','cryo'});  %Unitless
    varargout{1} = cat(1,varargout{1}, {'ice_ab_grad'   , 0, 0.3, 0.01, 'icmlt_ablation_grad','cryo'});  %Units of m (melt) / m (elevation)
    
    if isfield(sCryo,'icdbr')
        varargout{1} = cat(1,varargout{1}, {'debris_scalar' , 0, 2, 1, 'icmlt_ablation_grad','cryo'});  %Unitless
    end
    
    return
else
    tasScl = find_att(varargin{1}.coef,'ice_tas_offset');  %Unitless
    abGrad = find_att(varargin{1}.coef,'ice_ab_grad');  %
    
%     if isfield(sCryo,'tsn')
%         tsnPwr = find_att(varargin{1}.coef,'tsn_pwr');
%     end
    
    if isfield(sCryo,'icdbr')
        debScl = find_att(varargin{1}.coef,'debris_scalar');  %Unitless
    end
    
    sMeta = varargin{1};
end

% if ~isfield(sCryo, tic)
%     sCryo.tic = nan(size(sCryo.snw),'single');
% end


%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'icx')
    indIce = find( sCryo.icx ~= 0);
else
    indIce = [];
end

%Only calculate ice melt once per year (necessary assumption in ablation
%gradient model)

indDtCurr = find(ismember(sMeta.dateRun(:,1),sMeta.dateCurr(1), 'rows') == 1);
if isempty(indDtCurr)
    error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
end


%Reset values on first day of year
if isequal(sMeta.dateRun(indDtCurr(1),2:end), sMeta.dateCurr(2:end))
    sCryo.tasaccum = zeros(size(sCryo.snw),'single');
        sCryo.tasaccum(isnan(sCryo.icx)) = nan;
    sCryo.icdwe = zeros(size(sCryo.snw),'single');
        sCryo.icdwe(isnan(sCryo.icx)) = nan;
    sCryo.iclr = zeros(size(sCryo.snw),'single'); %Ice liquid release 
        sCryo.iclr(isnan(sCryo.icx)) = nan;
end
    
%Calculate accumulated air temperature every day of year:
sCryo.tasaccum = sCryo.tasaccum + squeeze(sAtm.tas(sAtm.indtas,:,:));

%Calculate year's melt on last day of year:
if isequal(sMeta.dateRun(indDtCurr(end),2:end), sMeta.dateCurr(2:end))
    %Average acculumated temperature:
    sCryo.tasaccum = sCryo.tasaccum / numel(indDtCurr);
    
    indIceMlt = intersect(find(sCryo.tasaccum >= tasScl), indIce);
    
    if ~isempty(indIceMlt)
        %Find max elevation of melting ice (Treat this as approximate ELA
        %elevation)
        %Use search window
        szWind = 20;
        elevMax = sHydro.dem;
        elevMax(setdiff((1:numel(sHydro.dem)), indIceMlt)) = 0;
        elevMax = min_2d_window(elevMax, szWind);

        %Find elevation difference (used in ablation gradient calculation):
        elevDiff = elevMax - sHydro.dem;
       
        %Calculate melt based on ablation gradient:
        sCryo.icdwe(indIceMlt) = abGrad*(elevDiff(indIceMlt)/100);
        
        %Physical constraints:
        sCryo.lhicme(isnan(sCryo.lhicme)) = 0;
        sCryo.lhicme(indIceMlt(sCryo.lhicme(indIceMlt) < 0)) = 0;
    
        %Ice release:
        sCryo.iclr(indIceMlt) = sCryo.icx(indIceMlt).*sCryo.icdwe(indIceMlt);
        
        %Update ice water equivalent based on melt:
        sCryo.icwe(indIceMlt) = sCryo.icwe(indIceMlt) - sCryo.lhicme(indIceMlt);
            sCryo.icwe(indIceMlt(sCryo.icwe(indIceMlt) < 0)) = 0;
    end
end

%Set change to nan outside region
sCryo.icdwe(isnan(sCryo.icx)) = nan;