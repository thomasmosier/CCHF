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

function facO = fac(dem, varargin)

%If flow-direction grid not input, calculate.
if length(varargin(:)) == 0
    fdrO = fdr(dem);
else
    fdrO = varargin{1};
end

%Initalize flow accumulation grid:
facO = zeros(size(fdrO));
%Calculate dimensions of DEM (used in converting vector indices to matrix)
dimDem = size(dem);

%Sort DEM, so that processing begins with first element and proceeds to
%lowest:
[~, ind] = sort(dem(:),'descend');

for ii = 1 : numel(dem)
    [r,c] = ind2sub(dimDem,ind(ii));
    
    if r == 1 || r == dimDem(1) || c == 1 || c == dimDem(2)
       continue
    end
    %key:
    %E  =   1
    %SE =   2
    %S  =   4
    %SW =   8
    %W  =  16
    %NW =  32
    %N  =  64
    %NE = 128
    switch fdrO(r,c)
        case 1
            facO(r  ,c+1) = facO(r,c) + 1;
        case 2
            facO(r+1,c+1) = facO(r,c) + 1;
        case 4
            facO(r+1,c  ) = facO(r,c) + 1;
        case 8
            facO(r+1,c-1) = facO(r,c) + 1;
        case 16
            facO(r  ,c-1) = facO(r,c) + 1;
        case 32
            facO(r-1,c-1) = facO(r,c) + 1;
        case 64
            facO(r-1,c  ) = facO(r,c) + 1;
        case 128
            facO(r-1,c+1) = facO(r,c) + 1;
    end
end

facO(isnan(dem)) = NaN;

% 
% cntr = 0;
% for ii = 2 : length(grid(1,:))-1 %Loop over columns
%     for jj = 2 : length(grid(:,1))-1
%     %key:
%     %E  =   1
%     %SE =   2
%     %S  =   4
%     %SW =   8
%     %W  =  16
%     %NW =  32
%     %N  =  64
%     %NE = 128
%     switch grid(jj,ii)
%         case 1
%             facO(jj  ,ii+1) = facO(jj,ii) + 1;
%         case 2
%             facO(jj+1,ii+1) = facO(jj,ii) + 1;
%         case 4
%             facO(jj+1,ii  ) = facO(jj,ii) + 1;
%         case 8
%             facO(jj+1,ii-1) = facO(jj,ii) + 1;
%         case 16
%             facO(jj  ,ii-1) = facO(jj,ii) + 1;
%         case 32
%             facO(jj-1,ii-1) = facO(jj,ii) + 1;
%         case 64
%             facO(jj-1,ii  ) = facO(jj,ii) + 1;
%         case 128
%             facO(jj-1,ii+1) = facO(jj,ii) + 1;
%     end
%     
%     cntr = cntr + 1;
%     if mod(cntr,100) == 0
%         disp(['FAC loop has run ' num2str(cntr) 'iterations.'])
%     end
%         
%         
%         
%     end
% end
