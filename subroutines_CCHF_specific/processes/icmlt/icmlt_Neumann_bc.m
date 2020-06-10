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

function varargout = icmlt_Neumann_bc(varargin)
%%MELT ICE:
    
global sCryo


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'tmp_depth', -2, 2, 0, 'icmlt_Neumann_bc','cryo'});
    
    return
else
	tmpDepth = find_att(varargin{1}.coef,'tmp_depth');
    sMeta = varargin{1};
end


kGlac = find_att(sMeta.global,'thermal_conduct_ice'); 
densW = find_att(sMeta.global,'density_water'); 
cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3


%Calculate conduction from bottom boundary condition of glacier upwards:
%conduct = -k * temp grad / distance
sCryo.hficbc = -kGlac*13/tmpDepth;


sCryo.lhicme  = zeros(size(sCryo.snw),'single');
    sCryo.lhicme(isnan(sCryo.icx)) = nan;
sCryo.icdwe   = zeros(size(sCryo.snw),'single');
    sCryo.icdwe(isnan(sCryo.icx)) = nan;
sCryo.iclr    = zeros(size(sCryo.snw),'single'); %Ice liquid release
    sCryo.iclr(isnan(sCryo.icx)) = nan;


%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'icx')
    eps = 0.001;
    indIceMlt = find( sCryo.icx ~= 0 & sCryo.hfneti > 0 & sCryo.snw < eps & sCryo.lhsnme < eps);
else
    indIceMlt = [];
end

%Calculate changes using residual energy at locations with ice:
if ~isempty(indIceMlt)
    %Calculate ice melt, which involves different degree-index factor

    sCryo.lhicme(indIceMlt) = time2sec(1,sMeta.dt,sMeta.dateCurr)*sCryo.hfneti(indIceMlt)/cLate;
        sCryo.lhpme(indIceMlt) = 0;
        
    %Physical constraints:
    sCryo.lhicme(isnan(sCryo.lhicme)) = 0;
    sCryo.lhicme(indIceMlt(sCryo.lhicme(indIceMlt) < 0)) = 0;
    
    sCryo.iclr(indIceMlt) = sCryo.icx(indIceMlt).*sCryo.lhicme(indIceMlt);
    sCryo.icdwe(indIceMlt) = -sCryo.lhicme(indIceMlt);
    
    %Update ice water equivalent based on melt:
    sCryo.icwe(indIceMlt) = sCryo.icwe(indIceMlt) - sCryo.lhicme(indIceMlt);
        sCryo.icwe(indIceMlt(sCryo.icwe(indIceMlt) < 0)) = 0;
end


%Set change to nan outside region
sCryo.icdwe(isnan(sCryo.icx)) = nan;

% %Set ice Temp (locations where no snow):
% sCryo.tic(sCryo.icx ~= 0 & sCryo.snw == 0) = 0;