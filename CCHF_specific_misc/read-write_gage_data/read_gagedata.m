function sObs = read_gagedata(path, lon, lat, varargin)

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
dataType = cell(0,1); %Record each datatype loaded (only used for writing CCHF formatted gagedata)


flagWrt = 0;
for ii = 1 : numel(path(:))
    [~, fileNm, ext] = fileparts(path{ii});
    
    if regexpbl(fileNm, 'cchf')
        sGageCurr = read_CCHF_gagedata(path{ii}, lon, lat, mask);
    elseif regexpbl(fileNm, 'usgs')
        flagWrt = 1; %If non-CCHF data being input, combine all data into one CCHF file after loading all.
        
        [dataGage, dateGage] = read_USGS_data(path{ii});

        %Check order of date array and reanarrange (year, month, day):
        indDate = [find(max(dateGage) > 50), find(max(dateGage) <= 12), find(max(dateGage) > 12 & max(dateGage) < 50)];
        dateGage = dateGage(:,indDate); 
        
        %User inputs coordinates of station since this not in USGS metadata:
        lonCurr = input(['Enter the longitude to the data at ' char(39) path{ii} char(39) ':' char(10)]);
        latCurr = input(['Enter the latitude to the data at ' char(39) path{ii} char(39) ':' char(10)]);

         sGageCurr = struct(...
            'lon',  lonCurr,...
            'lat',  latCurr,...
            'flow', dataGage,...
            'date', dateGage);
        dataType{end+1} = 'flow';
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
            dataType{end+1} = 'stake';
        elseif regexpbl(path{ii},{'snow','radar'},'and') || regexpbl(path{ii},'swe')
            dataType{end+1} = 'radar';
        end
    elseif regexpbl(path{ii}, 'MOD') && regexpbl(path{ii}, {'hdf','nc'}) %MODIS DATA
        flagWrt = 1;
        
        if all2d(isnan(lon)) && all2d(isnan(lat)) 
            warning('read_geodata:noCoord',['MODIS data not being '...
                'loaded because the coordinate input variable was not provided.']);
            continue
        end
        
        sGageCurr = import_MODIS_data(path{ii}, lon, lat, mask);
    	dataType{end+1} = 'casi'; %Stands for Covered Area Snow and Ice
    else
        error('read_CCHF_geodata:unknownType',['Currently the only '...
            'data that can be read must have ' char(39) 'USGS' char(39) ...
            ' or ' char(39) 'Arc' char(39) ...
            'in the title.  Program more gage data types.']);
    end


    namesPtsWrt = fieldnames(sGageCurr);
    
    %If current observation structure, 'sGageCurr', has fieldnames lat and
    %lon: add a substructure layer to the output structure;
    %otherwise: directly copy current structure as substructure 
    %('sGageCurr' composed of multiple substructures with 'lat', 'lon', etc.
    if ~any(strcmpi(namesPtsWrt,'lon')) && any(strcmpi(fieldnames(sGageCurr.(namesPtsWrt{1})), 'lon'))     
        %Transfer current gage substructure to output 
        for jj = 1 : numel(namesPtsWrt)
            namesObs = fieldnames(sObs);
            
            if ~regexpbl(namesObs,namesPtsWrt{jj})
                sObs.(namesPtsWrt{jj}) = sGageCurr.(namesPtsWrt{jj});
            else %substructure fieldname already exists, Append data:
                fieldsAppend = fieldnames(sObs.(namesPtsWrt{jj}));
                for kk = 1 : numel(fieldsAppend)
                    if ~isfield(sObs.(namesPtsWrt{jj}),fieldsAppend{kk})
                        sObs.(namesPtsWrt{jj}).(fieldsAppend{kk}) = sGageCurr.(namesPtsWrt{jj}).(fieldsAppend{kk});
                    else
                        sObs.(namesPtsWrt{jj}).(fieldsAppend{kk})(end+1:end+numel(sGageCurr.(namesPtsWrt{jj}).(fieldsAppend{kk}))) ...
                            = sGageCurr.(namesPtsWrt{jj}).(fieldsAppend{kk});
                    end
                end
            end
        end
        
%         for jj = 1 : numel(namesPtsWrt)
%             namesObs = fieldnames(sObs);
%             if ~regexpbl(namesObs,namesPtsWrt{jj})
%                 sObs.(namesPtsWrt{jj}) = sGageCurr.(namesPtsWrt{jj});
%             else %substructure fieldname already exists, create new unique fieldname
%                 cntrNewNm = 1;
%                 while any(strcmpi(namesObs,['pt' num2str(cntrNewNm)]))
%                     cntrNewNm = cntrNewNm + 1;
%                 end
%                 
%                 sObs.(['pt' num2str(cntrNewNm)]) = sGageCurr.(namesPtsWrt{jj});
%             end
%         end
    else %current gagedata lat and lon fields
        if ischar(sGageCurr.lat) && ischar(sGageCurr.lon) %lat and lon are strings and values are AVG
            if strcmpi(sGageCurr.lat,'avg') && strcmpi(sGageCurr.lon,'avg')
                if ~any(strcmpi(fieldnames(sObs),'avg'))
                   sObs.avg = sGageCurr;
                else
                    for jj = 1 : numel(namesPtsWrt)
                        if any2d(strcmpi(sObs.avg, namesPtsWrt{jj}))
                           warning('read_gagedata:multFieldsAvg', ...
                               ['There are multiple observations '...
                               'averaged over the spatial domain with '...
                               'the name ' char(39) namesPtsWrt{jj} ...
                               char(39) '. The second will overwrite the first.']);
                        end
                        sObs.avg.(namesPtsWrt{jj}) = sGageCurr.(namesPtsWrt{jj});
                    end
                    
%                     warning('read_gagedata:multFieldsAvg',...
%                         ['Multiple gagedata are area-averages. '...
%                         'Currently there is no check that the dates '...
%                         'between the two fields match.']);
                end
            elseif strcmpi(sGageCurr.lat,'all') && strcmpi(sGageCurr.lon,'all')
                if ~any(strcmpi(fieldnames(sObs),'all'))
                   sObs.all = sGageCurr;
                else
                    for jj = 1 : numel(namesPtsWrt)
                        if any2d(strcmpi(sObs.avg, namesPtsWrt{jj}))
                           warning('read_gagedata:multFieldsAll', ...
                               ['There are multiple observations '...
                               'distributed over the entire spatial domain with '...
                               'the name ' char(39) namesPtsWrt{jj} ...
                               char(39) '. The second will overwrite the first.']);
                        end
                        
                        sObs.all.(namesPtsWrt{jj}) = sGageCurr.(namesPtsWrt{jj});
                    end
                    
%                     warning('read_gagedata:multFieldsAll',...
%                         ['Multiple gagedata are for the entire domain. '...
%                         'Currently there is no check that the dates '...
%                         'between the two fields match.']);
                end
            else
               error('read_gagedata:coordStr',['The coordinate string '...
                   'combination ' char(39) sGageCurr.lon char(39) ...
                   ' and ' char(39) sGageCurr.lat char(39) ' are unknown.']); 
            end
        else %Name point based on linear index within model grid:
            [~, gageRow] = min(abs(sGageCurr.lat-lat));
            [~, gageCol] = min(abs(sGageCurr.lon-lon));
            if sGageCurr.lat > lat(1) - 0.5*abs(diff(lat(1:2))) || sGageCurr.lat < lat(end) + 0.5*abs(diff(lat(end-1:end))) 
                warning('CCHF_backbone:gagePtRow','The latitutde of the gauge point may be outside the area being modeled.');
            elseif  sGageCurr.lon < lon(1) - 0.5*abs(diff(lon(1:2))) || sGageCurr.lon > lon(end) + 0.5*abs(diff(lon(end-1:end))) 
                warning('CCHF_backbone:gagePtCol','The longitude of the gauge point may be outside the area being modeled.');
            end

            sObs.(['pt' num2str(round(sub2ind(szGrid, gageRow, gageCol)))]) = sGageCurr;
        end
    end
end
    

%     %If current observation structure, 'sGageCurr', has fieldnames lat and
%     %lon: add a substructure layer to the output structure;
%     %otherwise: directly copy current structure as substructure 
%     %('sGageCurr' composed of multiple substructures with 'lat', 'lon', etc.
%     if ~any(strcmpi(namesPtsWrt,'lon')) && any(strcmpi(fieldnames(sGageCurr.(namesPtsWrt{1})), 'lon'))     
%         for jj = 1 : numel(namesPtsWrt)
%             namesObs = fieldnames(sObs);
%             if ~regexpbl(namesObs,namesPtsWrt{jj})
%                 sObs.(namesPtsWrt{jj}) = sGageCurr.(namesPtsWrt{jj});
%             else %substructure fieldname already exists, create new unique fieldname
%                 cntrNewNm = 1;
%                 while any(strcmpi(namesObs,['pt' num2str(cntrNewNm)]))
%                     cntrNewNm = cntrNewNm + 1;
%                 end
%                 
%                 sObs.(['pt' num2str(cntrNewNm)]) = sGageCurr.(namesPtsWrt{jj});
%             end
%         end
%     else %current gagedata lat and lon fields
%         if ischar(sGageCurr.lat) && ischar(sGageCurr.lon) %lat and lon are strings and values are AVG
%             if strcmpi(sGageCurr.lat,'avg') && strcmpi(sGageCurr.lon,'avg')
%                 if ~any(strcmpi(fieldnames(sObs),'avg'))
%                    sObs.avg = sGageCurr;
%                 else
%                     for jj = 1 : numel(namesPtsWrt)
%                         sObs.avg.(namesPtsWrt{jj}) = sGageCurr.(namesPtsWrt{jj});
%                     end
%                     
%                     warning('read_gagedata:multFieldsAvg',...
%                         ['Multiple gagedata are area-averages. '...
%                         'Currently there is no check that the dates '...
%                         'between the two fields match.']);
%                 end
%             else
%                error('read_gagedata:coordStr',['The coordinate string '...
%                    'combination ' char(39) sGageCurr.lon char(39) ...
%                    ' and ' char(39) sGageCurr.lat char(39) ' are unknown.']); 
%             end
%         else
%             cntrNewNm = 1;
%             namesObs = fieldnames(sObs);
%             if ~isempty(namesObs)
%                 while any(strcmpi(namesObs,['pt' num2str(cntrNewNm)]))
%                     cntrNewNm = cntrNewNm + 1;
%                 end
%             end
%             sObs.(['pt' num2str(cntrNewNm)]) = sGageCurr;
%         end
%     end
% end

if all(~isnan(time))
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
        nms = fieldnames(sObs.(ptsCurr{ii}));

        %Remove non-data or date fields
        for kk = numel(nms) : -1 : 1
            if regexpbl(nms{kk},{'lat', 'lon', 'att', 'path'});
                nms(kk) = [];
            end
        end
%         nms(strcmpi(nms,'lat')) = [];
%         nms(strcmpi(nms,'latitude')) = [];
%         nms(strcmpi(nms,'lon')) = [];
%         nms(strcmpi(nms,'longitude')) = [];
%         nms(strcmpi(nms,'att')) = [];
%         nms(strcmpi(nms,'attributes')) = [];
%         nms(strcmpi(nms,'flag')) = [];

        %Find date fields:
        if any(strcmpi(nms,'date'))
            date{1} = sObs.(ptsCurr{ii}).date;
            if numel(date(:)) > 1
                date(2:end) = [];
            end
        elseif any(strcmpi(nms,'dateStart')) && any(strcmpi(nms,'dateEnd'))
            date{1} = sObs.(ptsCurr{ii}).dateStart;
            date{2} = sObs.(ptsCurr{ii}).dateEnd;
        else
            warning('read_gagedate:dateFind',['No date field found in sObs.' ptsCurr{ii} '.']);
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
        
        indKp = [min(indKp) max(indKp)];
        if ~isempty(indKp)
            for jj = 1 : numel(nms)
                if isfield(sObs.(ptsCurr{ii}), nms{jj}) && all(numel(sObs.(ptsCurr{ii}).(nms{jj})(:,1)) > indKp)
                    if numel(size(sObs.(ptsCurr{ii}).(nms{jj}))) == 2
                        sObs.(ptsCurr{ii}).(nms{jj}) ...
                            = sObs.(ptsCurr{ii}).(nms{jj})(indKp(1):indKp(2),:);
                    elseif numel(size(sObs.(ptsCurr{ii}).(nms{jj}))) == 3
                        sObs.(ptsCurr{ii}).(nms{jj}) ...
                            = sObs.(ptsCurr{ii}).(nms{jj})(indKp(1):indKp(2),:,:);
                    else
                        warning('read_gagedate:dim', ['The field ' char(39) nms{jj} ...
                            char(39) ' has ' num2str(numel(size(sObs.(ptsCurr{ii}).(nms{jj})))) ...
                            ' dimensions, which is unexpected.']);
                    end
                end
            end
        end
    end
end

if flagWrt == 1
    pathOut = path{1};
    [dirOut, ~, ~] = fileparts(pathOut);
    
    strTypes = blanks(0);
    for ii = 1 : numel(dataType)
        strTypes = [strTypes, dataType{ii}];
        if ii ~= 1 && ii ~= numel(dataType)
           strTypes = [strTypes, '-']; 
        end
    end

    pathOut = fullfile(dirOut,['CCHF_gage-' strTypes '.txt']);
    
    uiwait(msgbox(sprintf(['All the input gage data are being written'...
        ' to ' char(39) strrep(pathOut, '\', '/') char(39) '.' char(10) 'This file can '...
        'be used in place of all of the current gage files.']), ...
        '(Click OK to Proceed)','modal'));
    write_CCHF_gagedata(pathOut, sObs);
end