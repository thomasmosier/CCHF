function write_CCHF_gagedata(path, sObs)



namesPtsWrt = fieldnames(sObs);
sNum = '%.6g';

[dir, nameInput, ~] = fileparts(path);

% strComb = 'CCHF_format_combined';
% if ~regexpbl(nameInput, strComb)
%     nameInput = [nameInput, '_', strComb];
% end

path = fullfile(dir, [nameInput '.txt']);

%    fid = fopen(path,'w+');
fid = fopen(path,'wt');

%Write intro string
fprintf(fid,'%s\n', ['CCHF formatted gagedata from ' char(39) nameInput char(39) '.']);
fprintf(fid,'%s\n', 'Each point-variable pair of data will be written to this file but be separated by an empty line.');

%Write all substructures into same file, seperated by a blank line:
for ii = 1 : numel(namesPtsWrt(:))
    %IDENTIFY HEADER INFORMATION:
    namesCurr = fieldnames(sObs.(namesPtsWrt{ii}));
    
    %Do not know why 'fields' or 'indGage' are fields, but they occur sometimes, so remove them:
    namesCurr(strcmpi(namesCurr,'fields')) = [];
    namesCurr(strcmpi(namesCurr,'indGage')) = [];
    
    if regexpbl(namesCurr,'lon')
        lonCurr = sObs.(namesPtsWrt{ii}).lon;
        if isnan(lonCurr)
           lonCurr = 'spatial mean'; 
        end
        namesCurr(strcmpi(namesCurr,'lon')) = [];
    else
        error('write_CCHF_data:noLon',[namesPtsWrt{ii} ' does not have'...
            ' a longitude field and therefore cannot be written to file.']);
    end

    if regexpbl(namesCurr,'lat')
        latCurr = sObs.(namesPtsWrt{ii}).lat;
        if isnan(latCurr)
           latCurr = 'spatial mean'; 
        end
        namesCurr(strcmpi(namesCurr,'lat')) = [];
    else
        error('write_CCHF_data:noLat',[namesPtsWrt{ii} ' does not have'...
            ' a latitude field and therefore cannot be written to file.']);
    end

    if regexpbl(namesCurr,{'date','time'})
        if regexpbl(namesCurr,'date')
        	strDateSrch = 'date';
        else regexpbl(namesCurr,'time')
        	strDateSrch = 'time';
        end
        indDateNms = regexpi(namesCurr, strDateSrch);
        dateNms = namesCurr(~cellfun('isempty',indDateNms));  
        
        date = cell(numel(dateNms),1);
        for kk = 1 : numel(dateNms)
            date{kk} = sObs.(namesPtsWrt{ii}).(dateNms{kk});
            namesCurr(strcmpi(namesCurr,dateNms{kk})) = [];
        end
    else
        error('write_CCHF_data:noDate',[namesPtsWrt{ii} ' does not have'...
            ' a date field and therefore cannot be written to file.']);
    end
    
    if regexpbl(namesCurr,'attributes')
        attCurr = sObs.(namesPtsWrt{ii}).attributes;
        namesCurr(strcmpi(namesCurr,'attributes')) = [];
    else
        attCurr = blanks(0);
    end
    
    if regexpbl(namesCurr,'flag')
        flagCurr = sObs.(namesPtsWrt{ii}).flag;
        namesCurr(strcmpi(namesCurr,'flag')) = [];
    else
        flagCurr = blanks(0);
    end
    
    %If there is a path field, remove it:
    for kk = numel(namesCurr(:)) : -1 : 1
        if regexpbl(namesCurr{kk},'path')
            namesCurr(kk) = [];
        elseif regexpbl(namesCurr{kk},'type')
            namesCurr(kk) = [];
        end
    end


    %Loop over each perceived variable in the current sub-structure:
    for jj = 1 : numel(namesCurr(:))
        %Don't write current field if it has fewer elements the the date
        if numel(sObs.(namesPtsWrt{ii}).(namesCurr{jj})) < numel(date{1}(:,1))
            if size(sObs.(namesPtsWrt{ii}).(namesCurr{jj})) == 2
                warning('write_CCHF_gagedata:missingData', ...
                    ['Gage data corresponding to ' namesPtsWrt{ii} '.' ...
                    namesCurr{jj} ' is not being written because there ' ...
                    'appear to be missing entries.']);
                continue
            end
        end
        
        fprintf(fid,'\n'); %Write blank line between each section.
        
        %%WRITE HEADER:
        fprintf(fid,'%s\t', 'Variable:');
        fprintf(fid,'%s\n', namesCurr{jj});

        fprintf(fid,'%s\t', 'Longitude:');
        if isnumeric(lonCurr)
            fprintf(fid,[sNum '\n'], lonCurr);
        else
            fprintf(fid,'%s\n', lonCurr);
        end

        fprintf(fid,'%s\t', 'Latitude:');
        if isnumeric(latCurr)
            fprintf(fid,[sNum '\n'], latCurr);
        else
            fprintf(fid,'%s\n', latCurr);
        end

        if ~isempty(attCurr)
            fprintf(fid,'%s\t', 'Attributes:');
            fprintf(fid,'%s\n', attCurr);
        end
        
        %If data have 3 dim (e.g. MODIS) write path name rather than data
        if numel(size(sObs.(namesPtsWrt{ii}).(namesCurr{jj}))) == 3
            indPath = [];
            nmsTemp = fieldnames(sObs.(namesPtsWrt{ii}));
            for ll = 1 : numel(nmsTemp(:))
                if regexpbl(nmsTemp{ll}, 'path')
                    indPath = ll;
                end
            end
            if ~isempty(indPath)
                fprintf(fid,'%s\t', 'Path:');
                fprintf(fid,'%s\n', sObs.(namesPtsWrt{ii}).(nmsTemp{indPath}));
            else
                warning('write_CCHF_gagedata:noPath',['The path for ' ...
                    namesPtsWrt{ii} '.' namesCurr{jj} ' does not exist. ' ...
                    'Therefore it will not be written to the CCHF '...
                    'formatted observation data file.']);
            end
        else %Write data
            strDataHd = 'Value, ';
            strDataFrmt = [sNum ', '];
            dateWrite = [];
            for kk = 1 : numel(dateNms)
                if kk == 2
                    strDataHd = [strDataHd ', '];
                    strDataFrmt = [strDataFrmt ' '];
                end
                if numel(dateNms) == 1
                    strDateX = '';
                elseif numel(dateNms) == 2
                    if regexpbl(dateNms{kk},{'start','begin'})
                        strDateX = '_start';
                    elseif regexpbl(dateNms{kk},{'end','finish'})
                        strDateX = '_end';
                    else
                        error('write_CCHF_date:dateQual',[char(39) ...
                            dateNms{kk} char(39) ' is an unknown date qualifier.  Program for this.'])
                    end
                else
                    error('write_CCHF_date:MultDates','There are three date fields.  This has not been programmed for.');
                end

                switch numel(date{kk}(1,:))
                    case 1
                        strDataHd = [strDataHd, 'Year' strDateX];
                        strDataFrmt = [strDataFrmt, '%d,'];
                    case 2
                        strDataHd = [strDataHd, 'Year' strDateX ...
                            ', Month' strDateX];
                        strDataFrmt = [strDataFrmt, '%d, %d,'];
                    case 3
                        strDataHd = [strDataHd, 'Year' strDateX ...
                            ', Month' strDateX, ', Day', strDateX];
                        strDataFrmt = [strDataFrmt, '%d, %d, %d,'];
                    case 4
                        strDataHd = [strDataHd, 'Year' strDateX ...
                            ', Month' strDateX, ', Day', strDateX, ...
                            ', Hour', strDateX];
                        strDataFrmt = [strDataFrmt, '%d, %d, %d, %d,'];
                    case 5
                        strDataHd = [strDataHd, 'Year' strDateX ...
                            ', Month' strDateX, ', Day', strDateX, ...
                            ', Hour', strDateX, ', Minute', strDateX];
                        strDataFrmt = [strDataFrmt, '%d, %d, %d, %d, %d,'];
                    otherwise
                        error('write_CCHF_gagedata:unknownTimeFormat', ...
                            ['The case where the date vector has precision to ' ...
                            num2str(numel(date{kk}(1,:))) ' temporal '...
                            'subdivisions has not been programmed for.']);
                end

                dateWrite = [dateWrite, date{kk}];

                if ~isempty(flagCurr)
                    strDataHd = [strDataHd ', flag'];
                    strDataFrmt = [strDataFrmt ' %d,'];
                end
            end

            %%WRITE DATA AND DATE:
            %Write line describing columns:
            fprintf(fid,'%s\n', strDataHd);
            %Write data and dates:
            for kk = 1 : numel(dateWrite(:,1))
                %Write data:
                if ~isempty(flagCurr)
                    fprintf(fid,strDataFrmt, [sObs.(namesPtsWrt{ii}).(namesCurr{jj})(kk), dateWrite(kk,:), flagCurr(kk)]);
                else
                    fprintf(fid,strDataFrmt, [sObs.(namesPtsWrt{ii}).(namesCurr{jj})(kk), dateWrite(kk,:)]);
                end
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);