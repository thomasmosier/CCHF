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

function write_CCHF_gagedata(path, sObs)


varLon = 'longitude';
varLat = 'latitude';

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
    
    if regexpbl(namesCurr, varLon)
        lonCurr = sObs.(namesPtsWrt{ii}).(varLon);
        if isnan(lonCurr)
           lonCurr = 'spatial mean'; 
        end
        namesCurr(strcmpi(namesCurr,varLon)) = [];
    else
        error('write_CCHF_data:noLon',[namesPtsWrt{ii} ' does not have'...
            ' a longitude field and therefore cannot be written to file.']);
    end

    if regexpbl(namesCurr,varLat)
        latCurr = sObs.(namesPtsWrt{ii}).(varLat);
        if isnan(latCurr)
           latCurr = 'spatial mean'; 
        end
        namesCurr(strcmpi(namesCurr,varLat)) = [];
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
        indFlag = [];
        for kk = 1 : numel(namesCurr)
            if regexpbl(namesCurr{kk}, 'flag')
                indFlag(end+1) = kk;
            end
        end
        if numel(indFlag) == 1
            flagCurr = sObs.(namesPtsWrt{ii}).(namesCurr{indFlag});
            namesCurr(strcmpi(namesCurr,'flag')) = [];
        elseif numel(indFlag) > 1
            warning('writeCchfGageData:multFlags',['The pt ' namesPtsWrt{ii} ' has ' num2str(numel(indFlag)) ' flag fields. A maximum of one has been programmed for.']);
        end
    else
        flagCurr = blanks(0);
    end
    
    %If there is a path field, remove it:
    for kk = numel(namesCurr(:)) : -1 : 1
        if regexpbl(namesCurr{kk},'path')
            namesCurr(kk) = [];
        elseif regexpbl(namesCurr{kk},'site')
            namesCurr(kk) = [];
        elseif regexpbl(namesCurr{kk},'type')
            namesCurr(kk) = [];
        end
    end


    %Loop over each perceived variable in the current sub-structure:
    for jj = 1 : numel(namesCurr(:))
        dateNmsCurr = dateNms;
        
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
        
        %If data have 3 dim (e.g. MODIS or geodetic mass balance grid)
        %write path name rather than array
        if numel(size(sObs.(namesPtsWrt{ii}).(namesCurr{jj}))) == 3 || all(size(sObs.(namesPtsWrt{ii}).(namesCurr{jj})) ~= 1)
            indPath = [];
            nmsTemp = fieldnames(sObs.(namesPtsWrt{ii}));
            for ll = 1 : numel(nmsTemp(:))
                if regexpbl(nmsTemp{ll}, [namesCurr{jj} '_path'])
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
            
            if regexpbl(dateNmsCurr,{'start','begin','end','finish'})
                for zz = numel(dateNmsCurr(:)) : -1 : 1
                    if ~regexpbl(dateNmsCurr{zz}, namesCurr{jj})
                        dateNmsCurr(zz) = [];
                    end
                end  
            end
            
            if numel(dateNmsCurr(:)) > 2
                error('write_CCHF_date:MultDates','There are at least three date fields. This has not been programmed for.');
%             elseif isempty(dateNmsCurr)
%                 error('write_CCHF_date:noDates','There are no date fields. This has not been programmed for.');
            end
            
            for kk = 1 : numel(dateNmsCurr)
                if kk == 2
                    strDataHd = [strDataHd ', '];
                    strDataFrmt = [strDataFrmt ' '];
                end
                if numel(dateNmsCurr) == 1
                    strDateX = '';
                elseif numel(dateNmsCurr) == 2
                    if regexpbl(dateNmsCurr{kk},{'start','begin'})
                        strDateX = '_start';
                    elseif regexpbl(dateNmsCurr{kk},{'end','finish'})
                        strDateX = '_end';
                    else
                        error('write_CCHF_date:dateQual',[char(39) ...
                            dateNmsCurr{kk} char(39) ' is an unknown date qualifier.  Program for this.'])
                    end
                else
                    error('write_CCHF_date:MultDates','There are at least three date fields. This has not been programmed for.');
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
            if ~isempty(dateWrite)
                for kk = 1 : numel(dateWrite(:,1))
                    %Write data:
                    if ~isempty(flagCurr)
                        fprintf(fid,strDataFrmt, [sObs.(namesPtsWrt{ii}).(namesCurr{jj})(kk), dateWrite(kk,:), flagCurr(kk)]);
                    else
                        fprintf(fid,strDataFrmt, [sObs.(namesPtsWrt{ii}).(namesCurr{jj})(kk), dateWrite(kk,:)]);
                    end
                    fprintf(fid,'\n');
                end
            elseif ~isempty(sObs.(namesPtsWrt{ii}).(namesCurr{jj}))
                for kk = 1 : numel(sObs.(namesPtsWrt{ii}).(namesCurr{jj})(:,1))
                    %Write data:
                    fprintf(fid,strDataFrmt, sObs.(namesPtsWrt{ii}).(namesCurr{jj})(kk));
                    fprintf(fid,'\n');
                end
            end
        end
    end
end

fclose(fid);