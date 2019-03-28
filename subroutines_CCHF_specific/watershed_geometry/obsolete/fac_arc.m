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

function facO = fac_arc(dem, varargin)


%The result of Flow Accumulation is a raster of accumulated flow to each 
%cell, as determined by accumulating the weight for all cells that flow 
%into each downslope cell.

%Cells of undefined flow direction will only receive flow; they will not 
%contribute to any downstream flow. A cell is considered to have an 
%undefined flow direction if its value in the flow direction raster is 
%anything other than 1, 2, 4, 8, 16, 32, 64, or 128.

%The accumulated flow is based on the number of cells flowing into each 
%cell in the output raster. The current processing cell is not considered 
%in this accumulation.

%Output cells with a flow accumulation of zero are local topographic highs 
%and can be used to identify ridges.

%If the input flow direction raster is not created with the Flow Direction 
%tool, there is a chance that the defined flow could loop. If the flow 
%direction does loop, Flow Accumulation will go into an infinite loop and 
%never finish.

%Find sinks and use sink-filled DEM in processing
[~, dem] = is_sink(dem);

%If flow-direction grid not input, calculate.
if numel(varargin(:)) == 0
    fdrO = fdr_arc(dem);
else
    fdrO = varargin{1};
end

%Initalize flow accumulation grid:
facO = zeros(size(fdrO));
%Calculate dimensions of DEM (used in converting vector indices to matrix)
dimDem = size(dem);

%Sort DEM, so that processing begins with first element and proceeds to
%lowest:
[~, indDem] = sort(dem(:),'descend');

%Remove all nan elements:
indDem(isnan(dem(indDem))) = [];

%Loop over all non-nan elements
for ii = 1 : numel(indDem)
    [r,c] = ind2sub(dimDem,indDem(ii));
    
    %If flow is to off-grid cell, skip
    if r == 1 && (indFl(1) == 1 || indFl(1) == 4 || indFl(1) == 7 )
        continue
    end
    if r == dimDem(1) && (indFl(1) == 3 || indFl(1) == 6 || indFl(1) == 9 )
        continue
    end
    if c == 1 && (indFl(1) == 1 || indFl(1) == 2 || indFl(1) == 3 )
        continue
    end
    if c == dimDem(2) && (indFl(1) == 7 || indFl(1) == 8 || indFl(1) == 9 )
       continue
    end
    
    indSame = find(reshape(dem(r-1:r+1,c-1:c+1),[],1) == dem(r,c));
    if length(indSame) > 1
        indSame(indSame == 5) = [];
        
        for kk = 1 : numel(indSame)
            if ( indSame(kk) == 1 && fdr(indSame(kk)) == 2 ) || ( indSame(kk) == 2 && fdr(indSame(kk)) == 1 ) || ( indSame(kk) == 3 && fdr(indSame(kk)) == 128 ) || ( indSame(kk) == 4 && fdr(indSame(kk)) == 4 ) || (indSame(kk) == 6 && fdr(indSame(kk)) == 64) || (indSame(kk) == 7 && fdr(indSame(kk)) == 8) || ( indSame(kk) == 8 && fdr(indSame(kk)) == 16 ) || ( indSame(kk) == 9 && fdr(indSame(kk)) == 32 )
                indDem(ii:ii+1) = indDem(ii+1:ii);
                [r,c] = ind2sub(dimDem,indDem(ii));    
            end
        end
    end
        
    
%     if r == 11 && c == 15
%        keyboard 
%     end
    
    %Return indice of 3x3 grid (corresponding to reshape((1:9),3,3) )
    indFl = fdr_inv_arc(fdrO(r,c));
    

    
    if numel(indFl) == 1
        switch indFl
        %key: Direction, FDR Val, indice on 3x3
            %E  =   1;      8
            %SE =   2;      9
            %S  =   4;      6
            %SW =   8;      3
            %W  =  16;      2
            %NW =  32;      1
            %N  =  64;      4
            %NE = 128;      7
            case 8 %E
                facO(r  ,c+1) = facO(r  ,c+1) + facO(r,c) + 1;
            case 9 %SE
                facO(r+1,c+1) = facO(r+1,c+1) + facO(r,c) + 1;
            case 6 %S
                facO(r+1,c  ) = facO(r+1,c  ) + facO(r,c) + 1;
            case 3 %SW
                facO(r+1,c-1) = facO(r+1,c-1) + facO(r,c) + 1;
            case 2 %W
                facO(r  ,c-1) = facO(r  ,c-1) + facO(r,c) + 1;
            case 1 %NW
                facO(r-1,c-1) = facO(r-1,c-1) + facO(r,c) + 1;
            case 4 %N
                facO(r-1,c  ) = facO(r-1,c  ) + facO(r,c) + 1;
            case 7 %NE
                facO(r-1,c+1) = facO(r-1,c+1) + facO(r,c) + 1;
        end
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
