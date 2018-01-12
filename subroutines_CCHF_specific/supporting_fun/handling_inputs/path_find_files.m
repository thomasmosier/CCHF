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

function sPath = path_find_files(sPath, sMeta)

disp('Finding files to load during model runs. This may take several minutes (depends on number of files in directory).');

gcmRefP = [-1,-1,-1]; dateStartP = [-1,-1,-1]; dateEndP = [-1,-1,-1];
%Loop over each of the variables to load
for ll = 1 : numel(sMeta.varLd)
    %Skip if path not present
    if ~isfield(sPath,sMeta.varLd{ll}) || isempty(sPath.(sMeta.varLd{ll}))
        continue
    end
    
    %Find files in present path:
    fileNcTemp = dir(fullfile(sPath.(sMeta.varLd{ll}),'*.nc'));
    if isempty(fileNcTemp)
        fileNcTemp = dir(fullfile(sPath.(sMeta.varLd{ll}),'*.asc'));
    end
    if isempty(fileNcTemp)
        fileNcTemp = dir(fullfile(sPath.(sMeta.varLd{ll}),'*.txt'));
    end
    if isempty(fileNcTemp)
        warning('load_ts:uknownDataType',['Currently this function can '...
            'only work with netCDF, ASC, and TXT files.']);
        return
    end
    
    %Extract field:
    fileNcTemp = extractfield(fileNcTemp, 'name');

    %Test if function can read file name (based on CMIP5 convention):
    testNc = CMIP5_time(fileNcTemp{1});

    %Files are in format that can be read by 'CMIP5_time':
    if sum(isnan(testNc(1,:))) <= 1
            indTest = numel(sMeta.dateRun(1,:));
        testDate = CMIP5_time(fileNcTemp{1});
            testDate = testDate(:,~isnan(testDate(1,:)));
            indTest = min(numel(testDate(1,:)),indTest); 

        dateStart = nan(numel(fileNcTemp),indTest);
        dateEnd   = nan(numel(fileNcTemp),indTest);

        for kk = 1 : numel(fileNcTemp)
            testDate = CMIP5_time(fileNcTemp{kk});
            dateStart(kk,:) = testDate(1,1:indTest);
            dateEnd(kk,:) = testDate(2,1:indTest);
        end

        if regexpbl(fileNcTemp{1}, '.nc')
            attTime = ncinfo(fullfile(sPath.(sMeta.varLd{ll}), fileNcTemp{1}), 'time');
                attTime = squeeze(struct2cell(attTime.Attributes))';

            [gcmRef, gcmUnits] = NC_time_units(attTime);
            cal =  NC_cal(attTime);

            if ll == 1 || ~isequal(gcmRef, gcmRefP)
                daysModRun  = days_since(gcmRef, sMeta.dateRun, cal);
            end
            if ll == 1 || ~isequal(dateStart, dateStartP) || ~isequal(dateEnd, dateEndP)
                daysStart = days_since(gcmRef, dateStart, cal);
                daysEnd   = days_since(gcmRef, dateEnd, cal);
                if regexpbl(gcmUnits,'hour')
                    daysStart = daysStart / 24;
                    daysEnd = daysEnd / 24;
                end
            end
        else
            gcmRef = [0, 0, 0];
            
            testNc = CMIP5_time(fileNcTemp{1});
            indRem = find(isnan(testNc(1,:)));
            if ~isempty(indRem)
                gcmRef(   :, indRem) = [];
            end
            cal = 'gregorian';

            daysStart   = days_since(gcmRef, dateStart, cal);
            daysEnd     = days_since(gcmRef,   dateEnd, cal);
            daysModRun  = days_since(gcmRef, sMeta.dateRun, cal);
        end
        
        sPath.([sMeta.varLd{ll} 'File']) = cell(numel(sMeta.dateRun(:,1)),1);
        for kk = 1 : numel(daysModRun)
            indUse = intersect(find(daysModRun(kk) >= daysStart), find(daysModRun(kk) <= daysEnd));
            if numel(indUse) == 1
                sPath.([sMeta.varLd{ll} 'File']){kk} = fullfile(sPath.(sMeta.varLd{ll}), fileNcTemp{indUse});
            elseif numel(indUse) == 0
                sPath.([sMeta.varLd{ll} 'File']){kk} = blanks(0);
                error('path_file_find:nofile',['No climate file found for ' ...
                    num2str(sMeta.dateRun(kk,1)) '-' num2str(sMeta.dateRun(kk,2)) ...
                    '. May need to change start or end times.']);
            else
                error('path_file_find:multfile',['Multiple climate files found for ' ...
                    num2str(sMeta.dateRun(kk,1)) '-' num2str(sMeta.dateRun(kk,2)) '.' ]);
            end
        end
        
    else
        %Try underscore convention:
        indUnd = regexpi(fileNcTemp{1}, '_');
      	testDigit = zeros(numel(indUnd) + 1, 1);
        for ii = 1 : numel(testDigit)
            if ii == 1 
                if ~isnan(str2double(fileNcTemp{1}(1:indUnd(1)-1)))
                    testDigit(ii) = 1; 
                end
            elseif ii == numel(testDigit) 
                indExt = regexpi(fileNcTemp{1}, '\.');
                if ~isnan(str2double(fileNcTemp{1}(indUnd(end)+1:indExt(end)-1)))
                    testDigit(ii) = 1; 
                end
            else
                if ~isnan(str2double(fileNcTemp{1}(indUnd(ii-1)+1:indUnd(ii)-1)))
                    testDigit(ii) = 1;
                end
            end
        end
        
        %Give up if fewer than two numbers:
        if sum(testDigit) < 2
            warning('path_find_files:unknown naming convention.', ...
                ['Rest of model may not work because the ' char(39) ...
                'path_find_files' char(39) ' function is aborting early.']);
            continue
        else
            dateFiles = nan(numel(fileNcTemp),sum(testDigit));
            
            %loop over files:
            for kk = 1 : numel(fileNcTemp)
                indUnd = regexpi(fileNcTemp{kk}, '_');
                indExt = regexpi(fileNcTemp{kk}, '\.');
                if sum(testDigit) == 2
                    dateFiles(kk,:) = ...
                        [str2double(fileNcTemp{kk}(indUnd(end-1)+1:indUnd(end)-1)), ...
                        str2double(fileNcTemp{kk}(indUnd(end)+1:indExt(end)-1))];
                elseif sum(testDigit) == 3
                    dateFiles(kk,:) = ...
                        [str2double(fileNcTemp{kk}(indUnd(end-2)+1:indUnd(end-1)-1)), ...
                        str2double(fileNcTemp{kk}(indUnd(end-1)+1:indUnd(end)-1)), ...
                        str2double(fileNcTemp{kk}(indUnd(end)+1:indExt(end)-1))];
                else
                    warning('path_find_files:subdaily',['Subdaily '...
                        'time-series inputs in ASCII format have not been programmed for.']);
                end
            end
            
            
            sPath.([sMeta.varLd{ll} 'File']) = cell(numel(sMeta.dateRun(:,1)),1);
            for kk = 1 : numel( sMeta.dateRun(:,1))
                indUse = find(ismember(dateFiles, sMeta.dateRun(kk,:),'rows') == 1);
                
                if numel(indUse) == 1
                    sPath.([sMeta.varLd{ll} 'File']){kk} = fullfile(sPath.(sMeta.varLd{ll}), fileNcTemp{indUse});
                elseif numel(indUse) == 0
                    sPath.([sMeta.varLd{ll} 'File']){kk} = blanks(0);
                    error('path_file_find:nofile',['No climate file found for ' ...
                        num2str(sMeta.dateRun(kk,1)) '-' num2str(sMeta.dateRun(kk,2)) ...
                        '. May need to change start or end times.']);
                else
                    error('path_file_find:multfile',['Multiple climate files found for ' ...
                        num2str(sMeta.dateRun(kk,1)) '-' num2str(sMeta.dateRun(kk,2)) '.' ]);
                end
            end
            %Put special indicator stating resolution of file:
            if sum(testDigit) == 2
                sPath.([sMeta.varLd{ll} 'Res']) = 'month';
            elseif  sum(testDigit) == 3
                sPath.([sMeta.varLd{ll} 'Res']) = 'day';
            end
        end
    end
    
%     %Check to see if it appears that daily time-series files are grouped
%     %together (e.g. in monhtly files)
%     nTotFiles = numel(sPath.([sMeta.varLd{ll} 'File']));
%     nUniqFiles = numel(unique(sPath.([sMeta.varLd{ll} 'File'])));
%     if nTotFiles / nUniqFiles > 10
%         warning('pathFileFind:groupedFiles', ['The ratio of time-steps to ' ...
%             sMeta.varLd{ll} ' files is approx. ' num2str(round(nTotFiles / nUniqFiles)) ...
%             '. It is recommended that individual files are used for each ' ...
%             'time step because this will greatly increase model evaluation speed.']);
%     end
end