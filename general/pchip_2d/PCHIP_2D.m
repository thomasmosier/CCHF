function Zi = PCHIP_2D(Xa,Ya,Za,Xi,Yi)
% Copyright 2013, 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the GlobalClimateData Downscaling Package.
% 
% The GlobalClimateData Downscaling Package is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The Downscaling Package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%Uses Matlab's built in 1D pchip function in 2D.  Input and output formats
%same as for Matlab's 'interp2' function.

warning('OFF','MATLAB:chckxy:IgnoreNaN');

%Reform Xa, Ya, Xi, and Yi (Xa to colum vectors and Ya to row vectors):
if length(Xa(:,1)) > 1 && length(Xa(1,:)) == 1 
    Xa = Xa';
end
if length(Xi(:,1)) > 1 && length(Xi(1,:)) == 1 
    Xi = Xi';
end
if length(Ya(1,:)) > 1 && length(Ya(:,1)) == 1 
    Ya = Ya';
end
if length(Yi(1,:)) > 1 && length(Yi(:,1)) == 1 
    Yi = Yi';
end

%Ensure X vectors go from W to E and Y vectors from N to S:
if Xa(1,1) > Xa(1,2)
    Xa = fliplr(Xa);
    Za = fliplr(Za);
end
if Ya(1,1) < Ya(2,1)
    Ya = flipud(Ya);
    Za = flipud(Za);
end
flipOut = [0,0];
if Xi(1,1) > Xi(1,2)
    Xi = fliplr(Xi);
    flipOut(1) = 1;
end
if Yi(1,1) < Yi(2,1)
    Yi = flipud(Yi);
    flipOut(2) = 1;
end

%Initialize temporary matrix (row ind of original, columns of interp ind):
ZiTemp = nan( length(Xa(:,1)), length(Xi(1,:)) );

for jj = 1 : length( Za(:,1) )  %Interpolate over each row of input data

    %interpolate, using Matlab's built in 'pchip', which is 1D. 
    %at each of these points and save as row of output data.
    if length(Za(jj,:)) - sum(isnan(Za(jj,:))) < 2
        %Cannot perform interpolation with less than two points
        ZiTemp(jj, : ) = NaN; 
        warning('PCHIP_2D:NaN','Interpolation output for current row is NaN because less than two input points in the row are non-NaN.');
    else
        ZiTemp(jj, : ) = pchip( Xa(1,:), Za(jj,:), Xi(1,:));
    end
end

%Initialize output matrix:
Zi = nan(length(Yi(:,1)), length(Xi(1,:)));

%Interpolate each column of ZiTemp:
for kk = 1 : length( ZiTemp(1,:) )
    if length(ZiTemp(:,kk)) - sum(isnan(ZiTemp(:,kk))) < 2
        Zi( : , kk) = NaN;
        warning('PCHIP_2D:NaN','Interpolation output for current column is NaN because less than two input points in the column are non-NaN.');
    else
        Zi( : , kk) = pchip( Ya(:,1), ZiTemp(:,kk), Yi(:,1) );
    end
end

%If Xi or Yi were changed, flip output:
if flipOut(1) == 1
    Zi = fliplr(Zi);
end
if flipOut(2) == 1
    Zi = flipud(Zi);
end

warning('ON','MATLAB:chckxy:IgnoreNaN');
end