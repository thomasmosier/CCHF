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

function varargout = snmass_step(varargin)

    
global sCryo

if isempty(varargin(:))
	varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'tmp_melt_offset', 0,   4,    1, 'heat_melt_threshold','cryo'}); %Units of depth melt
        
    return
else
%     tmpOffset = find_att(varargin{1}.coef,'tmp_melt_offset'); 
    sMeta = varargin{1};
end

%This conversion complements a similar one in 'ETI_Pellicciotti' and
%'simple_degree'. Makes units consistent with SETI
densW = find_att(sMeta.global,'density_water'); 
cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3

%Initialize release field:
sCryo.lhsnme = zeros(size(sCryo.snw),'single');
% sCryo.lhicme  = zeros(size(sCryo.snw),'single');
if ~isfield(sCryo, 'snlw')
    sCryo.snlw = zeros(size(sCryo.snw),'single');
end



%Calculate melt potential using simple degree indec formulation (units of m): 
sCryo.lhpme = time2sec(1,sMeta.dt,sMeta.dateCurr)*sCryo.hfnet/cLate;
% % %Set melt potential to zero at cells where temperature less than threshold:
% sCryo.lhpme(sCryo.lhpme < 0 ) = 0;


%Refreeze liquid water in snow if energy input is negative and there is
%liquid water
indLiqFrz = find(sCryo.lhpme < 0 & sCryo.snlw > 0 & sCryo.snw > 0); 
if ~isempty(indLiqFrz)
    %Melt potential is negative at these grid cells
    
    %Indices where there's not sufficient energy to freeze all liquid water
    indFrzMax = find(sCryo.snlw(indLiqFrz) > -sCryo.lhpme(indLiqFrz) );
    if ~isempty(indFrzMax)
        sCryo.snw(  indLiqFrz(indFrzMax)) =   sCryo.snw(indLiqFrz(indFrzMax)) - sCryo.lhpme(indLiqFrz(indFrzMax));
        sCryo.snlw( indLiqFrz(indFrzMax)) =  sCryo.snlw(indLiqFrz(indFrzMax)) + sCryo.lhpme(indLiqFrz(indFrzMax));
        sCryo.sndwe(indLiqFrz(indFrzMax)) = sCryo.sndwe(indLiqFrz(indFrzMax)) - sCryo.lhpme(indLiqFrz(indFrzMax));
        sCryo.lhpme(indLiqFrz(indFrzMax)) = 0;
    end
    
    %Indices where there is sufficient energy to freeze all liquid:
    if isempty(indFrzMax)
        indFrzNMax = indLiqFrz;
    else
        indFrzNMax = setdiff(indLiqFrz, indLiqFrz(indFrzMax));
    end
    if ~isempty(indFrzNMax)
        sCryo.snw(  indFrzNMax) =   sCryo.snw(indFrzNMax) + sCryo.snlw(indFrzNMax);
        sCryo.sndwe(indFrzNMax) = sCryo.sndwe(indFrzNMax) + sCryo.snlw(indFrzNMax);
        sCryo.lhpme(indFrzNMax) = sCryo.lhpme(indFrzNMax) + sCryo.snlw(indFrzNMax);
        sCryo.snlw(indFrzNMax) = 0;
    end
end


%Melt potential all goes into melting snow but is limited by amount of snow
%available
%The remaining input heat energy goes into melting snow
indMelt = find( sCryo.lhpme > 0);
if ~isempty(indMelt)
    %Grid cells where melt is limited by amount of snow
    indMeltAll = find(sCryo.lhpme(indMelt) >= sCryo.snw(indMelt));
    sCryo.lhsnme(indMelt(indMeltAll)) = sCryo.snw(indMelt(indMeltAll));
    
    %Grid cells where melt is limited by energy
    if isempty(indMeltAll)
        indMeltSome = indMelt;
    else
        indMeltSome = setdiff(indMelt, indMelt(indMeltAll));
    end
    if ~isempty(indMeltSome)
        sCryo.lhsnme(indMeltSome) = sCryo.lhpme(indMeltSome);
    end
    
    %Balance values (including adding melted snow to liquid water content):
    sCryo.lhpme(indMelt) = sCryo.lhpme(indMelt) - sCryo.lhsnme(indMelt);
    sCryo.snlw( indMelt) = sCryo.snlw( indMelt) + sCryo.lhsnme(indMelt);
    sCryo.snw(  indMelt) = sCryo.snw(  indMelt) - sCryo.lhsnme(indMelt);
    sCryo.sndwe(indMelt) = sCryo.sndwe(indMelt) - sCryo.lhsnme(indMelt);
end


%Set negative solid snow values to 0:
sCryo.snw(sCryo.snw < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.snlw(sCryo.snlw < 0 ) = 0;