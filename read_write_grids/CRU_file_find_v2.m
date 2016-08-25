function [filesUseTS, nFiles] = CRU_file_find_v2(dirTS, mnth, yrData, metVar)
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


%This function locates CRU time-series files in a directory and selects
%which to use based upon naming convection and input dates.
%Version 2 has been amended to work with CRU NetCDF files in addition to
%ASCII files.

%If NC CRU files and program works as expected, 'nFiles' = 1; however, if ASCII files 'nFiles' is number of files that match the  

yrMin = yrData(1);
yrMax = yrData(2);
nYrs  = yrData(3);

if ~isempty(dir([dirTS filesep 'cru_*.nc']))
    
    filesTS = dir([dirTS filesep 'cru' '*' '.nc']);

    for ii = length(filesTS):-1:1
        if isempty(regexpi(filesTS(ii).name, metVar))
            filesTS(ii) = [];
        end
    end
    
    filesUseTS = cell(1,1);
    
    nFiles = length(filesTS);
    
    if nFiles > 1
        filesUseTS = cellstr(filesTS(1).name);
        nFiles = 1;
       
        warning('CRU_file_find:NCmult',['Multiple CRU NetCDF files ' ...
            'were detected for ' metVar '.  ' cellstr(filesUseTS{1}) ...
            ' has been chosen for use.  Check that this is correct.']);
    elseif nFiles == 1
        filesUseTS = cellstr(filesTS.name);
        nFiles = 1;
    else
        error('CRU_file_find:NCempty',['CRU NetCDF files were detected;'...
            'however, none seem to be for the correct variable.']);
    end
    
    
elseif ~isempty(dir([dirTS filesep 'cru' '*' '.asc']))
    %Find files named in the CRU V2-V3 format for the specified month(ii)
    filesTS = dir([dirTS filesep '*_' num2str(mnth) '.asc']);

    %This finds all of the files in 'filesCRU' that fall within our desired 
    %range of years.
        kk = 1;     %kk is counter used inside of loop
        filesUseTS = cell(nYrs,1);    %filesUseTs is a list of all 
                                        %time series files to be 
                                        %downscaled for the current
                                        %month.
    for jj = 1:length(filesTS)
        if str2double(filesTS(jj).name(27:30)) >= yrMin && str2double(filesTS(jj).name(27:30)) <= yrMax
            if strcmpi(filesTS(jj).name(23:25), climVar) || ( strcmpi(climVar,'tmn') && strcmpi(filesTS(jj).name(1:3), 'CRU') &&  strcmpi(filesTS(jj).name(23:25), 'tmp') )
                filesUseTS(kk) = cellstr(filesTS(jj).name);
                kk = kk + 1;
            end
        end
    end

    nFiles = kk - 1; 
else
    error('CRU_file_find:empty','No CRU time-series files were identified in either NetCDF or ASCII format.');
end


