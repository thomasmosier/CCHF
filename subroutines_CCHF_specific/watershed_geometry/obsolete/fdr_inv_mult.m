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

function ind = fdr_inv_mult(val,varargin)
%Val = values output by fdr_val
%Indice = vector indice of 3x3 grid (i.e. indice if 3x3 grid converted to
%reshape((1:9),[3 3])



arcVal = [1,2,4,8,16,32,64,128];
%key: Direction, FDR Val, indice on 3x3
            %E  =   1;      8
            %SE =   2;      9
            %S  =   4;      6
            %SW =   8;      3
            %W  =  16;      2
            %NW =  32;      1
            %N  =  64;      4
            %NE = 128;      7

ind = [];

if isnan(val)
    return
elseif sum(val == arcVal) > 0
    ind = fdr_inv(val);
else
    while sum(val == arcVal) == 0
        indCurr = find(val >= arcVal, 1,'last');
        ind(end+1) = fdr_inv(arcVal(indCurr));
        val = val - arcVal(indCurr);
    end
    ind(end+1) = fdr_inv(arcVal(arcVal == val));
end


if ~isempty(varargin(:))
        crd = nan(numel(ind),2);
    
    for kk = 1 : numel(ind) 
        switch ind(kk)
            case 8
                crd(kk,:) = [ 0,  1];
            case 9
                crd(kk,:) = [ 1,  1];
            case 6
                crd(kk,:) = [ 1,  0];
            case 3
                crd(kk,:) = [ 1, -1];
            case 2
                crd(kk,:) = [ 0, -1];
            case 1
                crd(kk,:) = [-1, -1];
            case 4
                crd(kk,:) = [-1,  0];
            case 7
                crd(kk,:) = [-1,  1];
            case 5
                crd(kk,:) = [ 0,  0];
        end
    end
    
    ind = crd;
end

%%Slightly Slower version:
% ind = nan(8,1);
% cntr = 1;
% while sum(val == arcVal) == 0
% 
%     indCurr = find(val >= arcVal, 1,'first');
%     ind(cntr) = fdr_inv(arcVal(indCurr));
%     val = val - arcVal(indCurr);
%     
%     cntr = cntr + 1;
% end
% ind(cntr:8) = [];