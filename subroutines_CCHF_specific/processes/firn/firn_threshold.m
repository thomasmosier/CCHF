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

function varargout = firn_threshold(varargin)


if isempty(varargin(:))
    varargout{1} = cell(0,6);
    varargout{1} = cat(1, varargout{1}, {'firn_depth', 5, 50, 30, 'firn_threshold','cryo'}); %Unitless scalar
%     argout = cell(1,6);
%     argout(1,:) = {'timelag_pwr', 0.5, 1.5, 0.9561, 'tLag_Clark', 'routing'}; %Unitless scalar
    return
else
    sMeta = varargin{1};
    threshold = find_att(sMeta.coef,'firn_depth'); %Unitless scalar
end


%DOES NOT USE:
%Li, J., & Zwally, H. J. (2011). Modeling of firn compaction for estimating 
%ice-sheet mass change from observed ice-sheet elevation change. Annals of 
%Glaciology, 52(59), 1–7.


global sCryo

%Initialize firn2ic grid
sCryo.frn2ic = zeros(size(sCryo.snw));

%Find indices greater than threshold:
indCap = find(sCryo.snw > threshold);

if ~isempty(indCap)
    %Calculate
    sCryo.frn2ic(indCap) = sCryo.snw(indCap) - threshold;
    
    sCryo.icwe(indCap)  = sCryo.icwe(indCap) - sCryo.frn2ic(indCap);
    
    sCryo.icdwe(indCap) = sCryo.icdwe(indCap) + sCryo.frn2ic(indCap);
    sCryo.sndwe(indCap) = sCryo.sndwe(indCap) - sCryo.frn2ic(indCap);
    sCryo.snw(indCap) = threshold;
end
