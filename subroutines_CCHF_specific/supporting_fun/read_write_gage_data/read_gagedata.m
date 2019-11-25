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

function [sObs, varargout] = read_gagedata(path, lon, lat, varargin)

%Can handle cell array of paths.  Will load each as substructure of output
%and then combine based on location.

if ischar(path)
   path = {path}; 
end

time = nan(2,4);
mask = [];
if ~isempty(varargin(:))
    for ii = 1 : numel(varargin)
        if ischar(varargin{ii})
            switch varargin{ii}
%                 case 'box'
%                     lon = varargin{ii+1};
%                     lat = varargin{ii+2};
%                     coord = varargin{ii+1}; %[west, east, south, north]
                case 'time'
                    time = varargin{ii+1};
                case 'mask'
                    mask = varargin{ii+1};
            end
        end
    end
end

szGrid = [numel(lat),numel(lon)]; 

%Initialize output structure:
sObs = struct;
obsTypes = cell(0,1); %Record each datatype loaded (only used for writing CCHF formatted gagedata)


flagWrt = 0;
for ii = 1 : numel(path(:))
    [~, fileNm, ext] = fileparts(path{ii});
    
    if regexpbl(fileNm, 'cchf')
        [sGageCurr, varCurr] = read_CCHF_gagedata(path{ii}, lon, lat, mask);
        obsTypes(end+1 : end+numel(varCurr(:))) = varCurr;
        disp('finished loading cchf gagedata')
    elseif regexpbl(ext, {'.csv','.txt','.xls'}) && regexpbl(fileNm, {'flow','stream','discharge'})
        if regexpbl(fileNm, 'usgs')
            unitsIn = 'ft3/s';
        else
            unitsIn = 'm3/s';
        end
        
        warning('read_gagedata:flowUnits', ['It is being assumed that the units of ' fileNm ' are ' unitsIn '.']);
        
        flagWrt = 1; %If non-CCHF data being input, combine all data into one CCHF file after loading all.
        
        [dataGage, dateGage] = read_ts_table_data(path{ii}, unitsIn);

        %Check order of date array and reanarrange (year, month, day):
        indDate = [find(max(dateGage) > 50), find(max(dateGage) <= 12), find(max(dateGage) > 12 & max(dateGage) < 50)];
        dateGage = dateGage(:,indDate); 
        
        %User inputs coordinates of station since this not in USGS metadata:
        lonCurr = input(['Enter the longitude to the data at ' char(39) path{ii} char(39) ':' char(10)]);
        latCurr = input(['Enter the latitude to the data at ' char(39) path{ii} char(39) ':' char(10)]);

        varFlow = 'flow';
        sGageCurr = struct(...
        	'lon',  lonCurr,...
            'lat',  latCurr,...
            varFlow, dataGage,...
            [varFlow '_date'], dateGage);
        obsTypes{end+1} = 'flow';
    elseif  regexpbl(fileNm, 'arc') %POINT SHAPEFILE data copied from ArcMap (CSV file?)
        flagWrt = 1; %If non-CCHF data being input, combine all data into one CCHF file after loading all.
        
        if regexpbl(ext,'xls')
            sGageCurr = import_Arc_data(path{ii}, lon, lat);
        else
            error('read_Arc_geodata:unknownExt',['Currently the only '...
                'Arc data must be exported into an Excel spreadsheet.'...
                char(10) 'Program more file format types.']);
        end
        
        if regexpbl(path{ii},'stake')
            obsTypes{end+1} = 'stake';
        elseif regexpbl(path{ii},{'snow','radar'},'and') || regexpbl(path{ii},'swe')
            obsTypes{end+1} = 'radar';
        end
    elseif ~isempty(regexp(path{ii}, 'MOD')) && regexpbl(path{ii}, {'hdf','nc'}) %MODIS DATA
        flagWrt = 1;
        
        if all2d(isnan(lon)) && all2d(isnan(lat)) 
            warning('read_geodata:noCoord',['MODIS data not being '...
                'loaded because the coordinate input variable was not provided.']);
            continue
        end
        
        sGageCurr = import_MODIS_data(path{ii}, lon, lat, mask);
    	obsTypes{end+1} = 'casi'; %Stands for Covered Area Snow and Ice
    elseif regexpbl(ext, {'.asc','.txt', '.nc'}) && regexpbl(fileNm, {'geodetic','dDEM', 'Bolch'})
        flagWrt = 1;
        
        if all2d(isnan(lon)) && all2d(isnan(lat)) 
            warning('read_geodata:noCoord',['Geodetic mass balance data not being '...
                'loaded because the coordinate input variable was not provided.']);
            continue
        end
        
        sGageCurr = import_geodetic_data(path{ii}, lon, lat, mask);
        obsTypes{end+1} = 'geodetic'; %Mass balance at main grid
    else
        error('read_CCHF_geodata:unknownType',['Filename ' char(39) fileNm char(39) ...
            ' is not a recognized input type. Program more observation data type. '...
            'See if/else-if statements in this function for list of programmed data identifiers.']);
    end


    %If sGageCurr is not named based on 'pt', 'all', 'avg', add layer
    %struct or if doesn't have lat/lon fields
    namesPtsWrt = fieldnames(sGageCurr);
    if any(strcmpi(namesPtsWrt,'lon'))
        if isnumeric(sGageCurr.lon)
            pt = gage_pt_on_grid(sGageCurr.lon, sGageCurr.lat, lon, lat);
            strWrite = ['pt' num2str(pt)];
        elseif ischar(sGageCurr.lon) && strcmpi(sGageCurr.lon, 'all')
            strWrite = 'all';
        elseif ischar(sGageCurr.lon) && strcmpi(sGageCurr.lon, 'avg')
            strWrite = 'avg';
        else
            error('read_gagedata:lonType','The current lon type has not been programmed for.');
        end
        
        temp = sGageCurr;
        sGageCurr = struct;
        sGageCurr.(strWrite) = temp;
    end
    
    namesPtsWrt = fieldnames(sGageCurr);
    %Copy sCurrGage to observation structure array to be output from function
    for jj = 1 : numel(namesPtsWrt)
        namesObs = fieldnames(sObs);
        
        %Current name has a lon field and it is numeric (I.e. observation
        %exists at a point or set of points)
        if ~any(strcmpi(namesPtsWrt{jj}, namesObs)) %No such field exists, add
            sObs.(namesPtsWrt{jj}) = sGageCurr.(namesPtsWrt{jj});
        else %Append subfields
            fieldsAppend = fieldnames(sGageCurr.(namesPtsWrt{jj}));
            for kk = 1 : numel(fieldsAppend)
                fieldsObs = fieldnames(sObs.(namesPtsWrt{jj}));
               if ~any(strcmpi(fieldsAppend{kk}, fieldsObs))
                   sObs.(namesPtsWrt{jj}).(fieldsAppend{kk}) = sGageCurr.(namesPtsWrt{jj}).(fieldsAppend{kk});
               else
                   if ~any(strcmpi(fieldsAppend{kk}, {'lat', 'lon', 'latitude', 'longitude'}))
                       error('read_geodtata:unknownField', ...
                           ['The field ' fieldsAppend{kk} ' has not been programmed for.'])
                   end
               end
            end
        end
    end
end


%Write to file:
%This must come before cropping data to specific time bounds of current
%model run:
if flagWrt == 1
    if iscell(path) && numel(path(:)) > 1
        %Find common output directory:
        pathTemp = char(path(:));
        indSame = find(diff(pathTemp,2) == 0);
        if numel(indSame) < 2
            [dirOut, ~, ~] = fileparts(pathTemp);
            [dirOut, ~, ~] = fileparts(dirOut); %This moves the output directory one level up
        else
            indRtSame = find(diff(indSame) ~= 1, 1, 'first');
            if isempty(indRtSame)
                dirOut = path{1}(1:indSame(end));
            else
		if isequal(path{1}(indSame(indRtSame)), filesep)
		    dirOut = path{1}(1:indSame(indRtSame));
		else
		    dirOut = path{1}(1:indSame(indRtSame)+1);
		end
            end
        end
    else
        pathTemp = char(path);
        [dirOut, ~, ~ ] = fileparts(pathTemp);
    end
    
    strTypes = blanks(0);
    for ii = 1 : numel(obsTypes)
        strTypes = [strTypes, '-', obsTypes{ii}];
    end

    ccfhFileRt = ['CCHF_formatted_observations' strTypes];
    pathOutTxt = fullfile(dirOut, [ccfhFileRt '.txt']);
    pathOutMat = fullfile(dirOut, [ccfhFileRt '.mat']);
    
    disp(['All the input gage data are being written to ' char(39) ...
        strrep(pathOutTxt, '\', '/') char(39) '.' char(10) 'This file can '...
        'be used in place of all of the current gage files.' char(10) ...
        'This file can be moved to a preferred directory for future use.']);
    save(pathOutMat, 'sObs', 'obsTypes');
    write_CCHF_gagedata(pathOutTxt, sObs);
end

if nargout > 1
    if flagWrt == 1
        varargout{1} = pathOutTxt;
    else
        varargout{1} = '';
    end
end


%Crop data to specific time bounds of interest
if all2d(~isnan(time))
    date = cell(2,1);
    dateRef = [1,1,1,1];
    cal = 'gregorian';
    
    %Make sure time vectors are in order:
    timeChk = days_since(dateRef,time,cal);
    if timeChk(1) > timeChk(2)
        time = flipud(time);
    end
    
    if numel(time(1,:)) == 2
        time(1,3) = 1;
        time(2,3) = eomday(time(2,1),time(2,2));
    elseif numel(time(1,:)) == 1
        time(1,2:3) = [1, 1];
        time(2,2:3) = [12, eomday(time(2,1),12)];
    end
        
    ptsCurr = fieldnames(sObs);
    for ii = 1 : numel(ptsCurr)
        nmsAll = fieldnames(sObs.(ptsCurr{ii}));

        %Remove non-data or date fields
        for kk = numel(nmsAll) : -1 : 1
            if regexpbl(nmsAll{kk},{'lat', 'lon', 'att', 'path'})
                nmsAll(kk) = [];
            end
        end
        
        %Identify data fields:
        varCurr = cell(0,1);
        for kk = 1 : numel(obsTypes(:))
            if any(strcmpi(obsTypes{kk}, nmsAll))
                varCurr{end+1} = obsTypes{kk};
            end
        end
        
        %Loop over variables at current observation point:
        for kk = 1 : numel(varCurr)
            nmsCurr = nmsAll;
            for zz = numel(nmsCurr) : -1 : 1
                if ~regexpbl(nmsCurr{zz}, varCurr{kk})
                    nmsCurr(zz) = [];
                end
            end

            %Find date fields:
            dateFld = [varCurr{kk} '_date'];
            dateStFld = [varCurr{kk} '_dateStart'];
            dateEnFld = [varCurr{kk} '_dateEnd'];
            if any(strcmpi(nmsCurr, dateFld))
                date{1} = sObs.(ptsCurr{ii}).(dateFld);
                if numel(date(:)) > 1
                    date(2:end) = [];
                end
            elseif any(strcmpi(nmsCurr, dateStFld)) && any(strcmpi(nmsCurr, dateEnFld))
                date{1} = sObs.(ptsCurr{ii}).(dateStFld);
                date{2} = sObs.(ptsCurr{ii}).(dateEnFld);
            else
                if ~regexpbl(varCurr{kk}, 'flag')
                    warning('read_gagedate:dateFind',['No date field found in sObs at pt ' ptsCurr{ii} ' for field ' varCurr{kk} '.']);
                end
            end

            if numel(date{1}(1,:)) < numel(time(1,:))
                timeCurr = time( : , 1:numel(date{1}(1,:)));
            else
                timeCurr = time;
            end

            if numel(date{1}(1,:)) < numel(dateRef)
                dateRefCurr = dateRef(1 : numel(date{1}(1,:)));
            else
                dateRefCurr = dateRef;
            end


            %Find indices of dates to keep:
            if numel(date) == 1
                daysObs = days_since(dateRefCurr, date{1},cal);
                daysTime = days_since(dateRefCurr, timeCurr,cal);

                indKp = find(daysObs >= daysTime(1) & daysObs <= daysTime(2));
            else %Includes start and end date:
                daysObsStart = days_since(dateRefCurr, date{1},cal);
                daysObsEnd = days_since(dateRefCurr, date{2},cal);
                daysTime = days_since(dateRefCurr, timeCurr,cal);

                indKp = find(daysObsStart >= daysTime(1) & daysObsEnd <= daysTime(2));
            end

            indTBnds = [min(indKp) max(indKp)];
            if ~isempty(indTBnds)
                for jj = 1 : numel(nmsCurr)
                    if isfield(sObs.(ptsCurr{ii}), nmsCurr{jj}) && all(numel(sObs.(ptsCurr{ii}).(nmsCurr{jj})(:,1)) > indTBnds)
                        szCurr = size(sObs.(ptsCurr{ii}).(nmsCurr{jj}));
                        if numel(szCurr) == 2 && numel(lat) ~= szCurr(1)
                            sObs.(ptsCurr{ii}).(nmsCurr{jj}) ...
                                = sObs.(ptsCurr{ii}).(nmsCurr{jj})(indTBnds(1):indTBnds(2),:);
                        elseif numel(szCurr) == 3
                            sObs.(ptsCurr{ii}).(nmsCurr{jj}) ...
                                = sObs.(ptsCurr{ii}).(nmsCurr{jj})(indTBnds(1):indTBnds(2),:,:);
                        elseif numel(szCurr) == 2 && szCurr(1) == numel(lat) && szCurr(2) == numel(lon) 
                            sObs.(ptsCurr{ii}).(nmsCurr{jj}) = sObs.(ptsCurr{ii}).(nmsCurr{jj});
                        else
                            warning('read_gagedate:dim', ['The field ' char(39) nmsCurr{jj} ...
                                char(39) ' has ' num2str(numel(size(sObs.(ptsCurr{ii}).(nmsCurr{jj})))) ...
                                ' dimensions, which is unexpected.']);
                        end
                    end
                end
            end
        end
    end
end
