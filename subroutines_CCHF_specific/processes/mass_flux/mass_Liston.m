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

function varargout = mass_Liston(varargin)

    
global sCryo


if isempty(varargin(:))
	varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'gen_melt_pwr' , -2, 2, 0, 'mass_Liston','cryo'});
%     varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr' , -2, 2, -0.2, 'mass_Liston','cryo'});  %Units of mm/hr/Deg C
%     varargout{1} = cat(1,varargout{1}, {'tsn_pwr', -2,   2,    0, 'mass_Liston','cryo'});
    
%     if isfield(sCryo,'icdbr')
%         varargout{1} = cat(1,varargout{1}, {'debris_scalar' , 0, 2, 1, 'mass_Liston','cryo'});  %Unitless
%     end
    
    return
else
%     dgi = find_att(varargin{1}.coef,'gen_melt_pwr');  %Units of mm/hr/Deg C
%     dii = find_att(varargin{1}.coef,'ice_melt_pwr');  %Units of mm/hr/Deg C
    
%     if isfield(sCryo,'tsn')
%         tsnPwr = find_att(varargin{1}.coef,'tsn_pwr');
%     end
    
%     if isfield(sCryo,'icdbr')
%         debScl = find_att(varargin{1}.coef,'debris_scalar');  %Unitless
%     end
    
    sMeta = varargin{1};
end


%This conversion complements a similar one in 'ETI_Pellicciotti' and
%'simple_degree'. Makes units consistent with SETI
densW = find_att(sMeta.global,'density_water'); 
cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3
cSens = densW*find_att(sMeta.global,'heat_cap_snow');
minTemp = find_att(sMeta.global,'snow_temp_min'); 

%Initialize release field:
sCryo.lhsnme = zeros(size(sCryo.snw),'single');
if ~isfield(sCryo, 'lwsnl')
    sCryo.lwsnl = zeros(size(sCryo.snw),'single');
end



% if isequal(sMeta.dateCurr(1:2),[2000, 1])
%     keyboard
% end

%Calculate melt potential using simple degree indec formulation (units of m): 
sCryo.lhpme = time2sec(1,sMeta.dt,sMeta.dateCurr)*sCryo.hfnet/cLate;
sCryo.lhpme(sCryo.lhpme < 0) = 0;

sCryo.lhsnme = sCryo.lhpme;
sCryo.lhsnme = min(sCryo.lhsnme, sCryo.snw);

%Add melted snow to 'release' field and remove from 'solid':
sCryo.lwsnl = sCryo.lwsnl + sCryo.lhsnme;
sCryo.snw = sCryo.snw - sCryo.lhsnme;
sCryo.lhpme = sCryo.lhpme - sCryo.lhsnme;


%Modify snowpack temperature:
if isfield(sCryo,'hfsnc')
    %Negative means that heat is conducted from inner area to surface
    sCryo.tsn = sCryo.tsn - sCryo.hfsnc./(cSens*sCryo.snw);
else
   error('mass_enbal:noConduct','The enbal conservation of energy and mass formulation requires a heat conduction term.'); 
end


%ADJUST SNOW INTERNAL TEMP AND OTHER PARAMETERS:
%Set netagive solid snow values to 0:
sCryo.snw(sCryo.snw < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.lwsnl(sCryo.lwsnl < 0 ) = 0;

sCryo.tsn(sCryo.tsn > 0) = 0;
sCryo.tsn(sCryo.tsn < minTemp) = minTemp; 
sCryo.tsn(sCryo.lwsnl > 0.005*sCryo.snw) = 0;
sCryo.tsn(sCryo.snw == 0) = nan; 

% %Refreeze liquid in snow if internal temperature < 0
% %Don't freeze rain, because that water is already in snow liquid content
% %field
% indLiqFrz = find(sCryo.tsn < 0 & sCryo.lwsnl > 0 & sCryo.snw > 0); %Indices where there is positive cold content and there is liquid water in snow
% if ~isempty(indLiqFrz)
%     %Indices where there's more liquid water than energy to freeze
%     indFrzMax = find(sCryo.lwsnl(indLiqFrz) > -(cSens/cLate)*sCryo.snw(indLiqFrz).*sCryo.tsn(indLiqFrz) );
%     if ~isempty(indFrzMax)
%         frzMax = -(cSens/cLate)*sCryo.snw(indLiqFrz(indFrzMax)).*sCryo.tsn(indLiqFrz(indFrzMax));
%         sCryo.snw(indLiqFrz(indFrzMax)) = sCryo.snw(indLiqFrz(indFrzMax)) + frzMax;
%         sCryo.lwsnl(indLiqFrz(indFrzMax)) = sCryo.lwsnl(indLiqFrz(indFrzMax)) - frzMax;
%         sCryo.tsn(indLiqFrz(indFrzMax)) = 0;
%     end
%     
%     %Indices where there's sufficient energy to freeze all liquid:
%     if isempty(indFrzMax)
%         indFrzNMax = indLiqFrz;
%     else
%         indFrzNMax = setdiff(indLiqFrz, indLiqFrz(indFrzMax));
%     end
%     if ~isempty(indFrzNMax)
%         sCryo.snw(indFrzNMax) = sCryo.snw(indFrzNMax) + sCryo.lwsnl(indFrzNMax);
%         sCryo.tsn(indFrzNMax) = sCryo.tsn(indFrzNMax) + (cLate/cSens)*sCryo.lwsnl(indFrzNMax)./sCryo.snw(indFrzNMax);
%         sCryo.lwsnl(indFrzNMax) = 0;
%     end
% end





% %%MELT ICE:
% %if snowpack has ice field, allow additional melt at those locations:
% %Only finds indices where ice exists and where melt potential greater than
% %snowpack:
% if isfield(sCryo,'icx')
%     indIceMlt = find( sCryo.icx ~= 0 & sCryo.lhpme > 0);
% else
%     indIceMlt = [];
% end
% 
% %Calculate changes using residual energy at locations with ice:
% if ~isempty(indIceMlt)
%     %Calculate ice melt, which involves different degree-index factor
%     ratioMelt = sCryo.hfneti(indIceMlt)./sCryo.hfnet(indIceMlt);
%     maxRat = 4;
%     ratioMelt(ratioMelt > maxRat) = maxRat;
%     ratioMelt(ratioMelt < 0) = 0; %Can't have negative ice energy
%     sCryo.lhicme(indIceMlt) = 10^(dii)*sCryo.lhpme(indIceMlt).*ratioMelt;
%         sCryo.lhpme(indIceMlt) = 0;
%         
%     %Scale melt for debris covered glacier cells:
%     if isfield(sCryo,'icdbr')
%         indDebrisMlt = intersect(indIceMlt, find(~isnan(sCryo.icdbr))); 
%     	sCryo.lhicme(indDebrisMlt) = debScl*sCryo.lhicme(indDebrisMlt);
%     end
%     sCryo.lhicme( isnan(sCryo.lhicme)) = 0;
%     sCryo.iclr(indIceMlt) = sCryo.icx(indIceMlt).*sCryo.lhicme(indIceMlt);
%     sCryo.icdwe(indIceMlt) = -sCryo.lhicme(indIceMlt);
% end
% 
% 
% %ADJUST SNOW INTERNAL TEMP AND OTHER PARAMETERS:
% %Set netagive solid snow values to 0:
% sCryo.snw(sCryo.snw < 0 ) = 0;
% %Set negative snow liquid values to 0:
% sCryo.lwsnl(sCryo.lwsnl < 0 ) = 0;
% 
% sCryo.tsn(sCryo.tsn > 0) = 0;
% sCryo.tsn(sCryo.tsn < minTemp) = minTemp; 
% sCryo.tsn(sCryo.lwsnl > 0.005*sCryo.snw) = 0;
% sCryo.tsn(sCryo.snw == 0) = nan; 
% 
% 
% %Set ice Temp (locations where no snow):
% indBIce = find(sCryo.icx ~= 0 & sCryo.snw == 0);
%     indBIce3 = ind2_3d(size(sAtm.tas),indBIce, sAtm.indtas);
% %     sCryo.tsn(indBIce) = 0;
% %     sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tas(indBIce3)));
% %Temprature of bare ice is min(0, air temp):
% if isfield(sAtm,'tasmin')
%     sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tasmin(indBIce3)));
% else
%     sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tas(indBIce3)));
% end