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

function varargout = mass_step_tas(varargin)

    
global sCryo sAtm

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'tmp_melt_offset', 0,   4,    1, 'heat_melt_threshold','cryo'}); %Units of depth melt
        
    return
else
    tmpOffset = find_att(varargin{1}.coef,'tmp_melt_offset'); 
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
% %Set melt potential to zero at cells where temperature less than threshold:
% sCryo.lhpme(squeeze(sAtm.tas(sAtm.indtas,:,:)) <= tMelt ) = 0;
sCryo.lhpme(squeeze(sAtm.tas(sAtm.indtas,:,:)) <= tmpOffset ) = 0;
sCryo.lhpme(sCryo.lhpme < 0 ) = 0;



%Melt potential all goes into melting snow but is limited by amount of snow
%available
indMelt = find( sCryo.lhpme > 0);
if ~isempty(indMelt)
    sCryo.lhsnme(indMelt) = sCryo.lhpme(indMelt);
    indMaxSnow = find(sCryo.lhsnme(indMelt) > sCryo.snw(indMelt));
    sCryo.lhsnme(indMelt(indMaxSnow)) = sCryo.snw(indMelt(indMaxSnow));
    
    sCryo.lhpme(indMelt) = sCryo.lhpme(indMelt) - sCryo.lhsnme(indMelt);
    
    %Add melted snow to 'liquid' field and remove from 'solid':
    sCryo.snlw(indMelt) = sCryo.snlw(indMelt) + sCryo.lhsnme(indMelt);
    sCryo.snw(indMelt) = sCryo.snw(indMelt) - sCryo.lhsnme(indMelt);
end


%Set netagive solid snow values to 0:
sCryo.snw(sCryo.snw < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.snlw(sCryo.snlw < 0 ) = 0;