% Copyright 2016 Thomas M. Mosier

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
% along with the Downscabolling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function sGageCurr = import_geodetic_data(path, lon, lat, mask, varargin)

%Output is in units of water equivalent
unitsOut = 'meters WE';

varMB = 'geodetic';
flagWrt = 0;


if ~isempty(varargin(:)) && ~isempty(varargin{1}) && regexpbl(varargin{1}, {'interp','intrp','dist', 'avg','mean','average'})
    strType = varargin{1}; 
else
	%Require user input of average or spatially distributed MODIS SCA
    %measurements:
    uiwait(msgbox(sprintf(['Enter a character string to indicate whether '...
        'to average all geodetic observations over the domain (' char(39) ...
        'avg' char(39) ') or interpolate to model domain (' char(39) ...
        'interp' char(39) ').\n']), '(Click OK to Proceed)','modal'));
    strType = input(['Load ' char(39) 'avg' char(39) ' or ' char(39) 'interp' char(39) ' geodetic data?' char(10)],'s');
end



if ~isempty(mask)
    indNanMask = find(isnan(mask));
else
    warning('read_gagedata:maskNeeded',['Loading geodetic glacier mass balance data '...
        'requires a mask as an optional input argument. The ' ...
        'mask must be the spatial grid of the model and nan at '...
        'locations where the model will not be evaluated.']);
    
    indNanMask = [];
end


%Loading method depends on file type:
[fold, file, ext] = fileparts(path);
if regexpbl(ext, {'.asc','.txt','tif'})
    flagWrt = 1;
    %Elevation change or mass change?
    uiwait(msgbox(sprintf(['Enter a character string to indicate whether '...
            'the observations are in units of elevation change (' char(39) ...
            'elev' char(39) ') or mass change (' char(39) ...
            'mass' char(39) ').\n Data by Tobias Bolch are typically elevation.']), '(Click OK to Proceed)','modal'));
    strObsType = input(['Observation in units of ' char(39) 'elev' char(39) ' or ' char(39) 'mass' char(39) '?' char(10)],'s');

    if regexpbl(strObsType, 'elev')
        uiwait(msgbox(sprintf(['Enter the density to use for conversion (units of kg/m).\n' ...
            'Tobias Bolch recommends a value of 850 +/- 60 kg/m']), '(Click OK to Proceed)','modal'));
        iceDens = input(['Enter assumed density of ice/snow' char(10)]);
    else
        iceDens = 1000;
    end
    
    %Start dates and end dates:
    uiwait(msgbox(sprintf(['Enter the start date of the observation in the format ' ...
        char(39) 'dd/mm/yyyy' char(39) '.']), '(Click OK to Proceed)','modal'));
    strDateStart = input(['Enter start date: dd/mm/yyyy' char(10)], 's');

    uiwait(msgbox(sprintf(['Enter the ending date of the observation in the format ' ...
        char(39) 'dd/mm/yyyy' char(39) '.']), '(Click OK to Proceed)','modal'));
    strDateEnd = input(['Enter end date: dd/mm/yyyy' char(10)], 's');

    datesMB{1} = [str2double(strDateStart(7:10)), str2double(strDateStart(4:5)), str2double(strDateStart(1:2))];
    datesMB{2} = [str2double(  strDateEnd(7:10)), str2double(  strDateEnd(4:5)), str2double(  strDateEnd(1:2))];
    
    if regexpbl(ext, {'.asc','.txt'})
        [dataMB,hdrESRI,metaESRI] = read_ESRI(path);
        [latMB, lonMB] = ESRI_hdr2geo(hdrESRI,metaESRI);
    elseif regexpbl(ext, 'tif')
        varTemp = 'tempmb';
        sData = read_geotif(path, varTemp, nan(1,2), nan(1,2), 0, 'out');
        
        dataMB = sData.(varTemp);
        latMB = sData.('latitude');
        lonMB = sData.('longitude');
    end
    %Diagnostics:
%     figure; imagesc(dataMB); colorbar;
%     figure; histogram(dataMB)
%     min(dataMB(:))
%     max(dataMB(:))

    if ~isnan(iceDens)
        dataMB = (iceDens/1000)*dataMB;
    end
        
    [lonMeshMB, latMeshMB] = meshgrid(lonMB, latMB);
elseif regexpbl(ext, '.nc')
    sData = read_grid_NC(path);
    if any(size(sData.lat) == 1) && any(size(sData.lon) == 1)
        [lonMeshMB, latMeshMB] = meshgrid(sData.lon, sData.lat);
    else
        lonMeshMB = sData.lon;
        latMeshMB = sData.lat;
    end
    
    if isfield(sData, varMB)
        varMbLoad = varMB;
    elseif isfield(sData, 'geodetic')
        varMbLoad = 'geodetic';
    else
        error('importGeodeticData:unknownNetCDFFeild','The NetCDf field with geodetic data is not known.');
    end
    
    dataMB = sData.(varMbLoad);
    
    date = nc_days_2_date(sData, 'time_bnds');
    
    datesMB{1} = date(1,:);
    datesMB{2} = date(end,:);
    
    %Check units:
    if isfield(sData, ['att' varMbLoad])
        unitsIn = find_att(sData.(['att' varMbLoad]), 'units');
    else
        error('importGeodeticData:noUnits', ['No units are present in the geodetic mass balance netCDF file at ' path]);
    end
    
    if regexpbl(unitsIn, {'z', 'elev','ddem'})
        uiwait(msgbox(sprintf(['Enter the density to use for conversion (units of kg/m).\n' ...
            'Tobias Bolch recommends a value of 850 +/- 60 kg/m']), '(Click OK to Proceed)','modal'));
        iceDens = input(['Enter assumed density of ice/snow' char(10)]);
        
        if isnumeric() && ~isnan(iceDens)
            dataMB = (iceDens/1000)*dataMB;
        else
            error('importGeodeticData:iceDensityType','Ice density that was entered is incorrect.');
        end
    elseif ~regexpbl(unitsIn, {'w.e.', 'we'})
        error('importGeodeticData:unknownUnits', ['Geodetic mass balance units of ' unitsIn ' have not been programmed for']);
    else
    end
else
    error('import_geodetic_data:unknownType', ['The file type ' ext ...
        ' has not been programmed for.']);
end

%Test: Write dataMB to file:
% nDec = 2;
% hdrOut = ESRI_hdr(lonMB, latMB, 'corner');
% [foldWrt, fileWrt, ext] = fileparts(path);
% pathOutCurr = fullfile(foldWrt, [fileWrt, '_we_read_test.asc']);
% if exist(pathOutCurr, 'file')
%    delete(pathOutCurr); 
% end
% write_ESRI_v4(dataMB, hdrOut, pathOutCurr, nDec);  

%Create meshgrid for input lat/lon
[lonMesh, latMesh] = meshgrid(lon, lat);

%Interpolate Geodetic data to model grid:
if ~isequal(round2(lonMesh,3), round2(lonMeshMB,3)) || ~isequal(round2(latMesh,3), round2(latMeshMB,3))
    flagWrt = 1;
    
    %INTEGRATE OVER NEW COORDINATES:
    %NOTE: Set values to 0, then spatially integrate. This is equivalent to
    %the 'sCryo.iclr' field, which multiplies ice melt potential times
    %fractional glacier coverage
%     warning('import_geodetic_data:areaInt', ['The program is preparing ' ...
%         'to use an area weighted interpolation scheme to condition the '...
%         'geodetic mass balance data.' char(10) 'This process takes quite some time. ' ...
%         'The conditioned output will be written to netCDF file to avoid this processing step in future model runs.']);
%     
    dataMBGrid = geodata_area_wgt(lonMeshMB(1,:), latMeshMB(:,1), dataMB, lon, lat, 'nansum');

%FOR TESTING:
%     nDec = 2;
%     hdrOut = ESRI_hdr(lon, lat, 'corner');
%     [foldWrt, fileWrt, ext] = fileparts(path);
%     pathOutCurr = fullfile(foldWrt, [fileWrt, '_we_interp_test.asc']);
%     if exist(pathOutCurr, 'file')
%        delete(pathOutCurr); 
%     end
%     write_ESRI_v4(dataMBGrid, hdrOut, pathOutCurr, nDec); 

else
    dataMBGrid = single(dataMB);
end

%Set interpolated data to nan if mask present:
if ~isempty(indNanMask)
   dataMB(indNanMask) = nan;
end


if flagWrt == 1
    pathOut = fullfile(fold, [file, 'interp2grid', '.nc']);
    
    if exist(pathOut, 'file')
       delete(pathOut); 
    end
    
    disp(['The spatially integrated geodetic mass balance grid is being written to ' ...
        pathOut char(10) 'Use this path to save time in future model runs.'])
    
    timeBndsVec = [datesMB{1}; datesMB{2}];
    timeVec = round(nanmean([timeBndsVec(1,:); timeBndsVec(2,:)]));
    print_grid_NC_v2(pathOut, dataMBGrid, varMB, lon, lat, timeVec, timeBndsVec, 3, unitsOut);
else
    pathOut = path;
end


%Create output from MODIS time-series (either averaged or distributed):
if regexpbl(strType, {'avg','mean','average'})
     sGageCurr = struct(...
        'lon',  'avg',...
        'lat',  'avg',...
        varMB, mean2d(dataMB),...
        [varMB '_dateStart'], datesMB{1}, ...
        [varMB '_dateEnd'], datesMB{2} ...
        );
elseif regexpbl(strType, {'interp','intrp','dist'})
    %Ensure output has three dimensions:
    if ismatrix(dataMBGrid)
        if any(size(dataMBGrid) == 1)
            warning('importGeodeticData:unexpectedGeodeticDim','The size of at least one dimension of the geodetic grid is 1, which is not expected.');
        else
            dataTemp = dataMBGrid;
            dataMBGrid = nan([1, size(dataTemp)], 'single');
            dataMBGrid(1,:,:) = dataTemp;
        end
    end
    
    sGageCurr = struct(...
    	'lon',  'all',...
        'lat',  'all',...
        varMB, dataMBGrid,... 
        [varMB '_path'], pathOut,...
        [varMB '_dateStart'], datesMB{1}, ...
        [varMB '_dateEnd'], datesMB{2} ...
        ); 
else
   error('import_geodetic_data:unknownInterpretationType', [strType ...
       ' is an unknown geodetic mass balance data interpretation type.']);
end