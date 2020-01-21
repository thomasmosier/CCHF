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

function varargout = runoff_direct(varargin)

global sCryo sLand

if isempty(varargin(:))
	varargout{1} = cell(1,6);
    varargout{1}(1,:) = {'runoff_factor', 0, 1, 0.928, 'direct_runoff', 'routing'}; %percent runoff (unitless)
        
    return
else
    sMeta = varargin{1};
    rf = find_att(sMeta.coef,'runoff_factor');   
end


%Reset runoff field:
sLand.mrro(:) = 0;
    sLand.mrro(isnan(sCryo.icx)) = 0;

indX = [];
if isfield(sCryo, 'ice')
   indX = find(sCryo.icx > 0);
end

%Account for PET
if isfield(sLand,'pet')
    %Set ET so that it at most equals available liquid:
    indPET = find(ismember(sLand.datepet(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
    if isempty(indPET)
        error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
    end
    
    sLand.et = squeeze(sLand.pet(indPET,:,:));
        sLand.et(isnan(sCryo.icx)) = 0;
    
    %Set ET to 0 when snowpack or glaciers exist:
    sLand.et(sCryo.snw > 0) = 0;
    if ~isempty(indX)
        sLand.et(indX) = sCryo.icx(indX).*sLand.et(indX);
    end
else
    error('direct_runoff:noPET','There should be a potential evapotranspiration field.');
end

%Calculate total input water as snowmelt plus rain (will account for ET and direct runoff later):
%**Rain added to sCryo.release in snow melt functions
sLand.mrro = rf*(sLand.rnrf + sCryo.snlr + sCryo.iclr - sLand.et);
    sLand.mrro(sLand.mrro < 0 ) = 0;
