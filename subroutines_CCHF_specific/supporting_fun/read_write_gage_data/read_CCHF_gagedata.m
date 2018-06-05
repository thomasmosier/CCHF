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

function [sObs, varargout] = read_CCHF_gagedata(path, lon, lat, varargin)


%mask can be provided as variable input argument.
if ~isempty(varargin(:)) && ~isempty(varargin{1})
    mask = varargin{1}; 
else
	mask = [];
end

% szGrid = [numel(lat),numel(lon)]; 
% cntrName = 0;

sObs = struct;

[foldInput, nameInput, ~] = fileparts(path);


%Check name to see if data were written using 'write_CCHF_data'
if ~regexpbl(nameInput,'CCHF')
    warning('read_CCHF_data:inputName',['The data at ' char(39) regexprep(path,filesep,'-') char(39) ...
        ' does not appear to be in the standard CCHF data format.' ...
        char(10) 'This function will attempt to load the data but will likely crash.' ]);
end

pathObsSv = fullfile(foldInput, [nameInput '.mat']);

if exist(pathObsSv, 'file')
    disp('CCHF format observation data loaded from previous model run.');
    load(pathObsSv);
else
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
    varAll = cell(0,1);
    flagVarPair = 0;
    while flagVarPair == 0
        nameStCurr = 'NAN';
        varCurr = blanks(0);
        lonCurr = [];
        latCurr = [];
        attCurr = blanks(0);
        flgData = [];

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
                        if isnan(lonCurr) && regexpbl(strtrim(lineCurr(indLab+1:end)),{'mean','avg','average'})
                            lonCurr = 'avg';
                        elseif isnan(lonCurr) && regexpbl(strtrim(lineCurr(indLab+1:end)),{'all', 'model grid'})
                            lonCurr = 'all';
                        elseif ~isnumeric(lonCurr)
                            warning('read_CCHF_data:Lon',['The current '...
                                'longitude value is NaN but it is not '...
                                'clear that the current field is a spatially-averaged']);
                        end
                    case 'Latitude'
                        latCurr = str2double(strtrim(lineCurr(indLab+1:end)));
                        if isnan(latCurr) && regexpbl(strtrim(lineCurr(indLab+1:end)),{'mean','avg','average'})
                            latCurr = 'avg';
                        elseif isnan(latCurr) && regexpbl(strtrim(lineCurr(indLab+1:end)), {'all', 'model grid'})
                            latCurr = 'all';
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

        if ~isempty(varCurr)
            varAll{end+1} = varCurr;
        end
        
%         if regexpbl(varCurr, 'flag')
%             keyboard
%             while ~regexpbl(fgets(fid),'\w') && ~isnumeric(lineCurr)
%             end
%         end

        %Either load data from path or from lines within txt file:
        pathField = '';
        typeField = '';
        typeEval  = '';

        if ~isempty(pathLd)
            if regexpbl(pathLd, 'MOD')
%                 disp('interpolating MODIS')
                sGageCurr = import_MODIS_data(pathLd, lon, lat, mask, 'interp');

                nameStCurr = 'all';
                lonCurr = sGageCurr.lon;
                latCurr = sGageCurr.lat;
                dataCurr = sGageCurr.(varCurr);
                pathField = [varCurr '_path'];
                typeField = [varCurr '_type'];
                typeEval = sGageCurr.(typeField);

                %Determine if there is a start and end date or just one date:
                nmsTemp = fieldnames(sGageCurr);
                strDateStart = [varCurr '_dateStart'];
                strDateEnd = [varCurr '_dateEnd'];
                if any(strcmpi(nmsTemp, strDateStart)) && any(strcmpi(nmsTemp, strDateEnd)) 
                    dateLab = {'Start','End'};
                    datesCurr = [sGageCurr.(strDateStart), sGageCurr.(strDateEnd)];
                else
                    dateLab = {''};
                    datesCurr = sGageCurr.([varCurr '_date']);
                end

                varFlag1 = [varCurr '_flag'];
                varFlag2 = [varCurr 'flag'];
                if regexpbl(nmsTemp, varFlag1)
                    flgData = sGageCurr.(varFlag1);
                elseif regexpbl(nmsTemp, varFlag2)
                    flgData = sGageCurr.(varFlag2);
                else
                    flgData = []; 
                end
            elseif regexpbl(pathLd, {'geodetic','AST-KH9','AST-SRTM'})
%                 disp('AST_Kh9')
                sGageCurr = import_geodetic_data(pathLd, lon, lat, mask, 'interp');

                nameStCurr = 'all';
                lonCurr = sGageCurr.lon;
                latCurr = sGageCurr.lat;
                dataCurr = sGageCurr.(varCurr);
                pathField = [varCurr '_path'];
    %             typeField = [varCurr '_type'];
    %             typeEval = sGageCurr.(typeField);

                %Determine if there is a start and end date or just one date:
                nmsTemp = fieldnames(sGageCurr);
                strDateStart = [varCurr '_dateStart'];
                strDateEnd = [varCurr '_dateEnd'];
                if any(strcmpi(nmsTemp, strDateStart)) && any(strcmpi(nmsTemp, strDateEnd)) 
                    dateLab = {'Start','End'};
                    datesCurr = [sGageCurr.(strDateStart), sGageCurr.(strDateEnd)];
                else
                    dateLab = {''};
                    datesCurr = sGageCurr.([varCurr '_date']);
                end

            else
                [~, fileLd, extLd] = fileparts(pathLd);
                warning('read_CCHF_gagedata:spatDist',['A case has not been written to load ' ...
                    'spatially distributed data of the type present in ' ...
                    fileLd extLd '.']);
            end
        else
            %Close when end of file reached
            if isnumeric(lineCurr) %Means the end of the file has been reached.
%                 %Close file and return to invoking function
%                 fclose(fid); 
% 
%                 if nargout == 2
%                     varargout{1} = varAll;
%                 end
%                 disp('returning')
                break
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
                flgData = [];
            else
                datesCurr = A(:,2:end-1);
                flgData = A(:,end);
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
            elseif regexpbl(lonCurr,'all') && regexpbl(latCurr,'all')
                nameStCurr = 'all';
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
                    pt = gage_pt_on_grid(lonCurr, latCurr, lon, lat);

                    nameStCurr = ['pt' num2str(pt)];
                end
            end
        end


        %Add substructure to output structure:
        sObs.(nameStCurr).('lon') = lonCurr; 
        sObs.(nameStCurr).('lat') = latCurr; 

        if ~isempty(flgData)
            if regexpbl(nmsTemp, varFlag1)
                sObs.(nameStCurr).(varFlag1) = flgData;
            elseif regexpbl(nmsTemp, varFlag2)
                sObs.(nameStCurr).(varFlag2) = flgData;
            else
                sObs.(nameStCurr).(varFlag1) = flgData; 
            end
            
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
            sObs.(nameStCurr).([varCurr '_date' dateLab{kk}]) = datesCurr(:,ind(1):ind(2)); 
        end

        if ~isempty(attCurr)
            sObs.(nameStCurr).([varCurr 'attributes']) = attCurr;
        end

        %Close when end of file reached
        if isnumeric(lineCurr) %Means the end of the file has been reached.
            break
%             if nargout == 2
%                 varargout{1} = varAll;
%             end
%             %Close file and return to invoking function
%             fclose(fid); 
% 
%             return
        end
    end
    
    fclose(fid);
    
    disp(['CCHF format observation data is being saved to: ' pathObsSv]);
    save(pathObsSv, 'sObs', 'varAll');
end


%Ensure no date fields for "flag":
ptsLp = fieldnames(sObs);
for ii = 1 : numel(ptsLp)
    fldsCurr = fieldnames(sObs.(ptsLp{ii}));
    for jj = numel(fldsCurr) : -1 : 1
        if ~regexpbl(fldsCurr{jj}, 'flag')
            fldsCurr(jj) = '';
        end
    end
    
    if ~isempty(fldsCurr)
        for jj = numel(fldsCurr) : -1 : 1
            if ~isequal(fldsCurr{jj}(end-3:end), 'flag')
                
                sObs.(ptsLp{ii}) = rmfield(sObs.(ptsLp{ii}), fldsCurr{jj});
            end
        end
    end
end



if nargout == 2
    varargout{1} = varAll;
end