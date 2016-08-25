function sObs = read_CCHF_gagedata(path, lon, lat, varargin)


%mask can be provided as variable input argument.
if ~isempty(varargin(:)) && ~isempty(varargin{1})
    mask = varargin{1}; 
else
	mask = [];
end

szGrid = [numel(lat),numel(lon)]; 

sObs = struct;
cntrName = 0;

[~, nameInput, ~] = fileparts(path);


%Check name to see if data were written using 'write_CCHF_data'
if ~regexpbl(nameInput,'CCHF')
    warning('read_CCHF_data:inputName',['The data at ' char(39) regexprep(path,filesep,'-') char(39) ...
        ' does not appear to be in the standard CCHF data format.' ...
        char(10) 'This function will attempt to load the data but will likely crash.' ]);
end

%Open file:
fid = fopen(path, 'rt');

%First section of text is not important:
flagFirst = 0;
while flagFirst == 0
    lineCurr = fgets(fid);

    if ~regexpbl(lineCurr,'\w') %If the line has no alphabetic or numeric characters, it's a break.
        flagFirst = 1;
    end
end


%Next, read each variable-location pair:
%These pairs are seperated by a single blank line
flagVarPair = 0;
while flagVarPair == 0
    varCurr = blanks(0);
    lonCurr = [];
    latCurr = [];
    attCurr = blanks(0);

    %READ HEADER:
    flagSecHdr = 0;
    pathLd = '';
    while flagSecHdr == 0
        lineCurr = fgets(fid);

        if isnumeric(lineCurr) %Means the end of the file has been reached.
            %Close file and return to invoking function
            flagSecHdr = 1;
            indLab = [];
%             fclose(fid); 
%             
%             return
        else
            indLab = regexpi(lineCurr,':');
        end
        
        if ~isempty(indLab)
            switch strtrim(lineCurr(1:indLab-1))
                case 'Longitude'
                    lonCurr = str2double(strtrim(lineCurr(indLab+1:end)));
                    if isnan(lonCurr) && regexpbl(strtrim(lineCurr(indLab+1:end)),{'spatial','mean','avg','ave'})
                        lonCurr = 'avg';
                    elseif ~isnumeric(lonCurr)
                        warning('read_CCHF_data:Lon',['The current '...
                            'longitude value is NaN but it is not '...
                            'clear that the current field is a spatially-averaged']);
                    end
                case 'Latitude'
                    latCurr = str2double(strtrim(lineCurr(indLab+1:end)));
                    if isnan(latCurr) && regexpbl(strtrim(lineCurr(indLab+1:end)),{'spatial','mean','avg','ave'})
                        latCurr = 'avg';
                    elseif ~isnumeric(latCurr)
                        warning('read_CCHF_data:Lat',['The current '...
                            'latitude value is NaN but it is not '...
                            'clear that the current field is a spatially-averaged']);
                    end
                case 'Variable'
                    varCurr = strtrim(lineCurr(indLab+1:end));
                case 'Attributes'
                    attCurr = strtrim(lineCurr(indLab+1:end));
                case 'Path'
                    pathLd = strtrim(lineCurr(indLab+1:end));
                otherwise
                   flagSecHdr = 1;     
            end
        else
            flagSecHdr = 1;
            break  
        end
    end
    
    %Either load data from path or from lines within txt file:
    pathField = '';
    typeField = '';
    typeEval  = '';

    if ~isempty(pathLd)
        if regexpbl(pathLd, 'MOD')
            sGageCurr = import_MODIS_data(pathLd, lon, lat, mask, 'interp');
            
            nameStCurr = 'all';
            lonCurr = sGageCurr.lon;
            latCurr = sGageCurr.lat;
            dataCurr = sGageCurr.(varCurr);
            pathField = 'path_casi';
            typeField = 'type_casi';
            typeEval = sGageCurr.(typeField);
            
            %Determine if there is a start and end date or just one date:
            nmsTemp = fieldnames(sGageCurr);
            if any(strcmpi(nmsTemp,'dateStart')) && any(strcmpi(nmsTemp,'dateEnd')) 
                dateLab = {'Start','End'};
                datesCurr = [sGageCurr.dateStart, sGageCurr.dateEnd];
            else
                dateLab = {''};
                datesCurr = sGageCurr.date;
            end
            
            if regexpbl(nmsTemp,'flag')
                flag = sGageCurr.flag;
            else
                flag = []; 
            end
        else
            [~, fileLd, extLd] = fileparts(pathLad);
            warning('read_CCHF_gagedata',['A case has not been written to load ' ...
                'spatially distributed data of the type present in ' ...
                fileLd extLd '.']);
        end
    else
        
        %Close when end of file reached
        if isnumeric(lineCurr) %Means the end of the file has been reached.
            %Close file and return to invoking function
            fclose(fid); 

            return
        end
        
        %Remove blank spaces:
        lineCurr = strtrim(lineCurr);

        if strcmpi(lineCurr(end),',')
            nColCurr = numel(regexpi(lineCurr,','));
        else
            nColCurr = numel(regexpi(lineCurr,',')) + 1;
        end

        dataSpec = '';
        for ii = 1 : nColCurr 
            dataSpec = [dataSpec '%f,'];
        end

        %Test if following lines have same format:
        fTest = fid;
        lineTest = strtrim(fgets(fTest));

        if ~strcmpi(lineTest(end),',')
            dataSpec = dataSpec(1:end-1);
        end

        dataTest = cell2mat(textscan(lineTest,dataSpec));
        dataSpec = [dataSpec '\n'];


        A = [dataTest; reshape(fscanf(fid,dataSpec), nColCurr, [])'];

        dataCurr = A(:,1);
        if ~regexpbl(lineCurr,'flag')
            datesCurr = A(:,2:end);
            flag = [];
        else
            datesCurr = A(:,2:end-1);
            flag = A(:,end);
        end


        if regexpbl(lineCurr, {'start','end','begin','finish'})
            if regexpbl(lineCurr, {'start','end'},'and')
                dateLab = {'Start','End'};
            elseif regexpbl(lineCurr, 'start')
                dateLab = {'Start'};
            elseif regexpbl(lineCurr, 'end')
                dateLab = {'End'};
            else
                error('read_CCHF_data:unknownDateMarker','The current date marker has not been coded for.');
            end
        else
            dateLab = {''};
        end

        %PLACE ALL INFO FROM CURRENT iteration INTO SUBSTRUCTURE
        if regexpbl(lonCurr,'avg') && regexpbl(latCurr,'avg')
            nameStCurr = 'avg';
        else
            %Find current point names and see if the current data should be
            %added to one of these or if a new point should be created.
            namesExist = fieldnames(sObs);
            indFieldUse = [];
            for ii = 1 : numel(namesExist(:))
               if isnumeric(sObs.(namesExist{ii}).lon) && isnumeric(sObs.(namesExist{ii}).lat) && isequal(round2(sObs.(namesExist{ii}).lon,5), round2(lonCurr,5)) && isequal(round2(sObs.(namesExist{ii}).lat,5), round2(latCurr,5))
                  indFieldUse = ii; 
                  break
               end
            end

            if ~isempty(indFieldUse)
                nameStCurr = namesExist{indFieldUse};
            else
                [~, gageRow] = min(abs(latCurr-lat));
                [~, gageCol] = min(abs(lonCurr-lon));
    %             if latCurr > lat(1) + 0.5*abs(diff(lat(1:2))) || latCurr < lat(end) - 0.5*abs(diff(lat(end-1:end))) 
    %                 warning('CCHF_backbone:gagePtRow',['The latitutde of the gauge point (' num2str(latCurr) ') may be outside the area being modeled.']);
    %                 continue
    %             elseif  lonCurr < lon(1) - 0.5*abs(diff(lon(1:2))) || lonCurr > lon(end) + 0.5*abs(diff(lon(end-1:end))) 
    %                 warning('CCHF_backbone:gagePtCol',['The longitude of the gauge point (' num2str(lonCurr) ') may be outside the area being modeled.']);
    %                 continue
    %             end

                nameStCurr = ['pt' num2str(round(sub2ind(szGrid, gageRow, gageCol)))];
            end
        end
    end
    

    
    
    
    %Add substructure to output structure:
    sObs.(nameStCurr).('lon') = lonCurr; 
    sObs.(nameStCurr).('lat') = latCurr; 
    
    if ~isempty(flag)
       sObs.(nameStCurr).('flag') = flag; 
    end
    if ~isfield(sObs.(nameStCurr), varCurr)
        sObs.(nameStCurr).(varCurr) = dataCurr;
    else
        sObs.(nameStCurr).(varCurr)(end+1:end+numel(dataCurr)) = dataCurr;
    end
    if ~isempty(pathField)
        sObs.(nameStCurr).(pathField) = pathLd;
    end
    if ~isempty(typeField)
        sObs.(nameStCurr).(typeField) = typeEval;
    end
    
    
    sObs.(nameStCurr).(varCurr)(round(sObs.(nameStCurr).(varCurr)) == -9999) = NaN;
    
    
    for kk = 1 : numel(dateLab)
        if numel(dateLab) == 1
            ind = [1,numel(datesCurr(1,:))];
        elseif numel(dateLab) == 2
            ind = [numel(datesCurr(1,:))/numel(dateLab)*(kk-1)+1, numel(datesCurr(1,:))/numel(dateLab)*kk];
        else
            error('read_CCHF_data:ThreeDateMarkers','The current number of date markers has not been coded for.');
        end
        sObs.(nameStCurr).(['date' dateLab{kk}]) = datesCurr(:,ind(1):ind(2)); 
    end
    
    if ~isempty(attCurr)
        sObs.(nameStCurr).('attributes') = attCurr;
    end
    
    %Close when end of file reached
    if isnumeric(lineCurr) %Means the end of the file has been reached.
        %Close file and return to invoking function
        fclose(fid); 

        return
    end
end


fclose(fid);