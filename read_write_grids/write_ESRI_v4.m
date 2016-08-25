function write_ESRI_v4(data, hdr, fullFilepath, nDec)
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


%Writes ASCII text file using ESRI format
%'fullFilepath' defines the path that the file will be written to, including extension,
%'hdr' is a six element vector with header info following ESRI ASCII convention
%'data' is a rectangular data array
%'nDec' defines decimal point precision

%Check that data only has two dimensions:
if length(size(data)) > 2
   error('write_ESRI:dim3',['The data to be written to ' char(39) ...
       fullFilepath char(39) ' has more than two dimensions, ' ...
       'which is not allowed.']);
end

%Convert from NaN to numberic tag:
data( isnan(data) ) = hdr(6);

%Make output path and open file:
[filePath, ~, ext] = fileparts(fullFilepath);
if ~exist(filePath,'dir')
   mkdir(filePath); 
end
if isempty(ext) || ~regexpbl(ext,{'asc','txt'},'or')
    fullFilepath = [fullFilepath '.asc'];
end
if isunix
    fullFilepath = [filesep fullFilepath];
end
fiddata = fopen(fullFilepath, 'w');
    
%write header:
fprintf(fiddata,'%s %i\r\n',        'ncols', length(data(1,:) ) );
fprintf(fiddata,'%s %i\r\n',        'nrows', length(data(:,1) ));
fprintf(fiddata,'%s %.13f\r\n', 'xllcorner', hdr(3) );
fprintf(fiddata,'%s %.13f\r\n', 'yllcorner', hdr(4) );
fprintf(fiddata,'%s %.13f\r\n',  'cellsize', hdr(5) );
fprintf(fiddata,'%s %i\r\n', 'NODATA_value', hdr(6) );

%Create precision string for data:
if isnumeric(nDec) 
    if nDec > 0
        precStr = ['%.' num2str(nDec) 'f '];
    elseif nDec == 0
        precStr = '%i ';
    else
        error('write_ESRI:precNeg',['The number of decimals requested '...
            'is ' num2str(nDec) ', which is not possible.']);
    end
else
    error('wrtESRI:precision',['The value entered for the variable '...
        'specifying precision is not a number.']);
end
%round data to requested precision:
data = round2(data,nDec);

%Write data:
for ii = 1 : length(data(:,1))
    if issparse(data)
        fprintf(fiddata,precStr, full(data(ii,:)));
    else
        fprintf(fiddata,precStr, data(ii,:));
    end
    fprintf(fiddata,'\n');
end

fclose(fiddata);

       