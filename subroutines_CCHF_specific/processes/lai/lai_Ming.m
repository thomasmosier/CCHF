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

function varargout = lai_Ming(varargin)
%Light absorbing impurities, method by Ming et al. (2009)

%See:
%Ming, J., Xiao, C., Cachier, H., Qin, D., Qin, X., Li, Z., & Pu, J. 
%(2009). Black Carbon (BC) in the snow of glaciers in west China and its 
%potential effects on albedos. Atmospheric Research, 92(1), 114?123. 
%https://doi.org/10.1016/j.atmosres.2008.09.007


global sCryo sAtm

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
end

%Calculate total deposition
deposition = sAtm.BC(sAtm.('indBC'),:,:) + sAtm.DUST(sAtm.('indBC'),:,:);


dens = 0.3;

%Determine classification of grid cells that impact distribution of lai
indSnwfall = find(sCryo.dswe >= 0.1*dens);
indMlt = find(sCryo.dswe <= -0.1*dens);
indNeut = setdiff((1:numel(sCryo.dswe)), union(indSnwfall, indMlt));

%At Netural cells, deposition adds to concentration
sCryo.laiUpper(indNeut) = sCryo.laiUpper(indNeut) + deposition(indNeut);
%At snowfall cells, deposition adds to lower layer and any deposited BC in
%upper layer transfer to lower layer
sCryo.laiLower(indSnwfall) = sCryo.laiLower(indSnwfall) + deposition(indSnwfall);
sCryo.laiLower(indSnwfall) = sCryo.laiLower(indSnwfall) + sCryo.laiUpper(indSnwfall);
sCryo.laiUpper(indSnwfall) = 0;
%At snowmelt cells, deposition adds to upper layer and there is trans
sCryo.laiUpper(indMlt) = sCryo.laiUpper(indMlt) + deposition(indMlt);
laiTrans = sCryo.laiLower(indMlt).*(-sCryo.dswe(indMlt)/sCryo.swe(indMlt));
laiTrans(laiTrans < 0) = 0;
sCryo.laiUpper(indMlt) = sCryo.laiUpper(indMlt) + laiTrans;
sCryo.laiLower(indMlt) = sCryo.laiLower(indMlt) - laiTrans;


%Calculate albedo reduction factor based on concentration of lai in upper 
%layer and Ming's equition
sCryo.laiReduce = 0.0757*log10(sCryo.laiUpper)+0.0575;