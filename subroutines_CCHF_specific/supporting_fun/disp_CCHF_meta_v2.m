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

function disp_CCHF_meta_v2(sPath, sMeta)
%Displays select information from 'sPath{ii}', 'sMeta', and 'sHydro' to help
%ascertain and record parameters used in current model run.


%DISPLAY RUN TYPE AND DURATION
disp(['The CCHF hydrologic model is being run in ' sMeta.runType ' mode.'...
    char(10) 'The model run will begin ' num2str(sMeta.dateStart(2)) ...
    '/' num2str(sMeta.dateStart(1)) ' and will end ' num2str(sMeta.dateEnd(2)) ...
    '/' num2str(sMeta.dateEnd(1)) ', with a spinup period of ' ...
    num2str(sMeta.spinup) ' months prior to the start date.']);

%DISPLAY INPUT DATA
for ii = 1 : numel(sPath(:))
    disp(['FOR SITE ' num2str(ii) ' of ' num2str(numel(sPath(:))) ':']);
    disp(['Digital elevation model (DEM) data will be loaded from ' char(39) sPath{ii}.dem char(39) '.']);
    disp(['Precipitation time-series data will be loaded from ' char(39) sPath{ii}.pr char(39) '.']);
    disp(['Mean temperature time-series data will be loaded from ' char(39) sPath{ii}.tas char(39) '.']);
    if isfield(sPath{ii},'tasmin') && ~isempty(sPath{ii}.tasmin)
        disp(['Minimum temperature time-series data will be loaded from ' char(39) sPath{ii}.tasmin char(39) '.']);
    end
    if isfield(sPath{ii},'tasmax') && ~isempty(sPath{ii}.tasmax)
        disp(['Maximum temperature time-series data will be loaded from ' char(39) sPath{ii}.tasmax char(39) '.']);
    end
end
    
    
% %DISPLAY MODEL ROUTINE NAMES AND COEFFICIENTS USED
% disp(['The version of the CCHF model being implemented is ' char(39) sHydro.method char(39) '.']);
