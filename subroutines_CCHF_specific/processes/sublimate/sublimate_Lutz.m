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

function varargout = sublimate_Lutz(sHydro, varargin)

    
%Follows: 
%Lutz, A. F., Immerzeel, W. W., Kraaijenbrink, P. D. A., Shrestha, A. B., 
%& Bierkens, M. F. P. (2016). Climate Change Impacts on the Upper Indus 
%Hydrology: Sources, Shifts and Extremes. PLOS ONE, 11(11), e0165630. 
%https://doi.org/10.1371/journal.pone.0165630

global sCryo

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'sublimate_scl', 0,   2,    1, 'sublimate_Lutz','cryo'}); %Units of depth melt
        
    return
else
    scl = find_att(varargin{1}.coef,'sublimate_scl'); 
    sMeta = varargin{1};
end


%Initialize release field:
% sCryo.lhsnsb = zeros(size(sCryo.snw),'single');
sCryo.snsb = zeros(size(sCryo.snw),'single');

elevSubMag = sHydro.dem - 3000;
elevSubMag(elevSubMag < 0) = 0;

constant = 0.0015*10^(-3); %units = m/day


%Table 1 (Lutz et al., 2016) 
sCryo.snsb = scl*constant*(time2sec(1,sMeta.dt,sMeta.dateCurr)/86400)*elevSubMag;

sCryo.snsb(sCryo.snsb < 0) = 0;

sCryo.snw = sCryo.snw - sCryo.snsb;
sCryo.snw(sCryo.snw < 0) = 0;

