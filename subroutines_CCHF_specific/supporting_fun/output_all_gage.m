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

function [output, varargout] = output_all_gage(output, sObs, lon, lat)
%This function ensures all output fields present in observation are in the 
%output cell-array 

% disp('inside output_all_gage')
namesObs = fieldnames(sObs);

%Find the cells that each observation point corresponds to and
%average observation data if multiple points occur within same model cell.
        
%Define lat and lon edges and model grid size:
latEdg = [lat(1) + abs(0.5*diff(lat(1:2))), lat(end) - abs(0.5*diff(lat(end-1:end)))]; %[N, S]
lonEdg = [lon(1) - abs(0.5*diff(lon(1:2))), lon(end) + abs(0.5*diff(lon(end-1:end)))]; %[W, E]
szModel = [numel(lat), numel(lon)];

%FIND MODELLING GRID POINT FOR EACH OBSERVATION POINT
%Initialize array:
pts = nan(numel(namesObs),1); 
flagOutBnds = -888;
flagAvg = -777;

%Loop over all observation points
for ii = 1 : numel(namesObs)
    namesCurr = fieldnames(sObs.(namesObs{ii}));
    if ismember('lat',namesCurr)
        latCurr = sObs.(namesObs{ii}).('lat');
    elseif ismember('latitude',namesCurr)
        latCurr = sObs.(namesObs{ii}).('latitude');  
    end
    if ismember('lon',namesCurr)
        lonCurr = sObs.(namesObs{ii}).('lon');
    elseif ismember('longitude',namesCurr)
        lonCurr = sObs.(namesObs{ii}).('longitude');
    end

    if isnumeric(latCurr) && isnumeric(lonCurr) 
        if latCurr > latEdg(1) || latCurr < latEdg(2) || lonCurr < lonEdg(1) || lonCurr > lonEdg(2)
            pts(ii) = flagOutBnds;
        else
            %Find point in model grid that current point corresponds to:
            [~, indLatCurr] = min(abs(lat - latCurr));
            [~, indLonCurr] = min(abs(lon - lonCurr));

            pts(ii) = sub2ind(szModel, indLatCurr, indLonCurr);
        end
    elseif ischar(latCurr) && ischar(lonCurr)
        pts(ii) = flagAvg;
    else
        warning('output_all_gage:ptType','The current point is of an unknown type. This may cause issues.');
    end
end


%FOR ALL UNIQUE MODELING GRID POINTS WITH OBS, IDENTIFY OBSERVATIONS
%OF SAME TYPE.  FOR ALL REPEATED OBS, COMBINE DATA INTO SINGLE
%SUB-STRUCTURE AND THEN AVERAGE IF DATES REPEATED

%Loop over all unique grid points:
ptsUniq = unique(pts);

for ii = 1 : numel(ptsUniq) 
    %Identify pts matching current unique pt
    indPtsCurr = find(pts == ptsUniq(ii));

    %Cases to treat: (1) pt outside grid and (2) multiple
    %observations for same modelling grid point
    if ptsUniq(ii) == flagOutBnds %If pt outside of modeling grid, remove
        for jj = 1 : numel(indPtsCurr)
            sObs = rmfield(sObs, namesObs{indPtsCurr(jj)});
        end
    elseif numel(indPtsCurr) > 1 %If multiple points have same values, average all points
        %IDENTIFY ALL DATA VARIABLE NAMES PRESENT FOR CURRENT MODELLING GRID POINT:
        %initialize cell for variable names
        nmVarPts = cell(numel(indPtsCurr),1);

        %Find all variable names of observation points in current unique modeling grid pt:
        for jj = 1 : numel(indPtsCurr)
            %Remove non-data fields: 
            namesCurrPt = fieldnames(sObs.(namesObs{indPtsCurr(jj)}));
            for kk = numel(namesCurrPt) : -1 : 1
                if regexpbl(namesCurrPt{kk},{'lat','lon','date','time','att'})
                    namesCurrPt(kk) = [];
                end
            end

            if numel(namesCurrPt) == 1
                nmVarPts{jj} = namesCurrPt{1};
            else
                %This shouldn't happen if all working as envisioned
                error('output_all_gage:multDataField',['The current '...
                    'point has more than one data field but should only have one.']);
            end
        end 

        %Find unique variable names:
        nmDataUniq = unique(nmVarPts);

        %FOR EACH UNIQUE VARIABLE IN CURRENT GRID POINT, FIND IF 
        %THERE ARE MULTIPLE OBSERVATIONS
        %Loop over unique variable names:
        for jj = 1 : numel(nmDataUniq)
            %Find all observations with same variable name for
            %current modeling point:
            indVarCurr = find(strcmpi(nmDataUniq{jj},nmVarPts) == 1);

            %IF MULTIPLE OBSERVATIONS HAVE SAME VARIABLE NAME,
            %MERGE AND AVERAGE
            if numel(indVarCurr) > 1
                %Move all data to first point in current series
                for kk = 2 : numel(indVarCurr) %Loop over 2nd through last observation points with the same data name and in same modelling grid point
                    %Find number of points to transfer:
                    nTrans = numel(sObs.(namesObs{indPtsCurr(indVarCurr(kk))}).(nmDataUniq{jj}));

                    %TRANSFER DATE/TIME INFO TO 1ST UNIQUE POINT
                    %WITH OBSERVATION
                    %Get field names of current point (need for
                    %date / time)
                    nmDateTrans = fieldnames(sObs.(namesObs{indPtsCurr(indVarCurr(kk))}));
                    for ll = 1 : numel(nmDateTrans)
                        if regexpbl(nmDateTrans{ll},{'date','time'})
                            sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateTrans{ll})(end+1:end+nTrans,:) ...
                                = sObs.(namesObs{indPtsCurr(indVarCurr(kk))}).(nmDateTrans{ll})(:,:);
                        end

                        %Get latitude and longitude positions to average:
                        if regexpbl(nmDateTrans{ll},{'lat'})
                            sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateTrans{ll})(end+1,1) ...
                                = sObs.(namesObs{indPtsCurr(indVarCurr(kk))}).(nmDateTrans{ll})(:);
                        end
                        if regexpbl(nmDateTrans{ll},{'lon'})
                            sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateTrans{ll})(end+1,1) ...
                                = sObs.(namesObs{indPtsCurr(indVarCurr(kk))}).(nmDateTrans{ll})(:);
                        end
                    end

                    %TRANSFER DATA INFO TO 1ST UNIQUE POINT
                    %WITH OBSERVATION
                    sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDataUniq{jj})(end+1:end+nTrans,:) ...
                        = sObs.(namesObs{indPtsCurr(indVarCurr(kk))}).(nmDataUniq{jj})(:,:);

                    %REMOVE NON-FIRST SUBSTRUCTURE ENTRIES:
                    sObs = rmfield(sObs, namesObs{indPtsCurr(indVarCurr(kk))});
                end


                %SORT OBS DATES FROM FIRST ENTRY FOR EACH UNIQUE
                %VARIABLE BY TIME AND AVERAGE ENTRIES WITH
                %OVERLAPPING TIMES:
                nmDateSort = fieldnames(sObs.(namesObs{indPtsCurr(indVarCurr(1))}));
                indUniqDatesTemp = cell(numel(nmDateSort), 1);
                flagSort = 0; cntr = 0;
                for ll = 1 : numel(nmDateSort)
                    if regexpbl(nmDateSort{ll},{'date','time'})
                        if flagSort == 0
                            [sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll}), sortOrd] ...
                                = sortrows(sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll}));
                            flagSort = 1;
                        else
                            sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll}) ...
                                = sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll})(sortOrd,:);
                        end

                        cntr = cntr + 1;
                        %Find which dates are unique for each date entry:
                        [~, indUniqDatesTemp{cntr}, ~] ...
                            = unique(sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll}),'rows');
                    end
                end
                indUniqDatesTemp(cntr+1:end) = [];

                %SORT DATA USING ORDER DETERMINED BY DATES ABOVE:
                sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDataUniq{jj}) ...
                    = sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDataUniq{jj})(sortOrd);

                if numel(indUniqDatesTemp) == 1
                    %Find indices of repeated dates:
                    indRepDates = setdiff((1:numel(sortOrd)), indUniqDatesTemp{1});
                    indUniqDates = indUniqDatesTemp{1};
                elseif numel(indUniqDatesTemp) == 2
                    %Find indices of repeated dates:
                    indRepDates = cell(2,1);
                    indRepDates{1} = setdiff((1:numel(sortOrd)), indUniqDatesTemp{1});
                    indRepDates{2} = setdiff((1:numel(sortOrd)), indUniqDatesTemp{2});
                    indRepDates = intersect(indRepDates{1},indRepDates{2});
                    indUniqDates = setdiff((1:numel(sortOrd)), indRepDates);
                else
                    error('output_all_gage:nOut',[num2str(numel(indUniqDatesTemp)) ...
                        ' date indices were found, which is an unexpected number.']);
                end


                %Find each set of repeated dates:
                indDataAvg = cell(numel(indUniqDates),1);
                [ indDataAvg{:} ] = deal(nan(1,1));
                for ll = 1 : numel(indUniqDates)
                    indDataAvg{ll}(1) = indUniqDates(ll);

                    cntrUniq = 1;
                    while ~ismember(indUniqDates(ll) + cntrUniq, indUniqDates) && indUniqDates(ll) + cntrUniq <= numel(sortOrd)
                       indDataAvg{ll}(end+1) = indUniqDates(ll) + cntrUniq;
                      cntrUniq = cntrUniq + 1;
                    end
                end

                %Average data for repeated indices:
                dataAvg = nan(numel(indDataAvg),1);
                for ll = 1 : numel(indDataAvg)
                   dataAvg(ll) ...
                       = nanmean(sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDataUniq{jj})(indDataAvg{ll}));
                end
                sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDataUniq{jj}) = dataAvg;

                %Average lat and lon values for repeated indices:
                for ll = 1 : numel(nmDateSort)
                    if regexpbl(nmDateSort{ll},{'lat','lon'})
                        sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll}) ...
                           = nanmean(sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll}));
                    end
                end

                %Delete all repeated dates:
                if ~isempty(indRepDates)
                    for ll = 1 : numel(nmDateSort)
                        if regexpbl(nmDateSort{ll},{'date','time'})
                            sObs.(namesObs{indPtsCurr(indVarCurr(1))}).(nmDateSort{ll})(indRepDates) = [];
                        end
                    end
                end
            end
        end
    end
end
     
% 
% %If DEM nan, remove points from sObs:
% namesObs = fieldnames(sObs);
% for ii = 1 : numel(namesObs)
%     if isnumeric(sObs.(namesObs{ii}).lat) && isnumeric(sObs.(namesObs{ii}).lon)
%         %Find point in model grid that current point corresponds to:
%         [~, indLatCurr] = min(abs(lat - sObs.(namesObs{ii}).lat));
%         [~, indLonCurr] = min(abs(lon - sObs.(namesObs{ii}).lon));
% 
%         if isnan(varargin{1}.dem(indLatCurr,indLonCurr))
%            sObs = rmfield(sObs, namesObs{ii});
%         end
%     end
% end

%Write output
varargout{1} = sObs;


%Ensure all output fields present in observation are in the 
%output cell-array:
namesObs = fieldnames(sObs);
%Loop over all sub structures of observation data.  These correspond to different points.
for ii = 1 : numel(namesObs)
    namesCurr = fieldnames(sObs.(namesObs{ii}));
    
    if ismember('lat',namesCurr)
        latCurr = sObs.(namesObs{ii}).('lat');
        namesCurr(strcmpi(namesCurr,'lat')) = [];
    elseif ismember('latitude',namesCurr)
        latCurr = sObs.(namesObs{ii}).('latitude');
        namesCurr(strcmpi(namesCurr,'latitude')) = [];   
    end
    if ismember('lon',namesCurr)
        lonCurr = sObs.(namesObs{ii}).('lon');
        namesCurr(strcmpi(namesCurr,'lon')) = [];
    elseif ismember('longitude',namesCurr)
        lonCurr = sObs.(namesObs{ii}).('longitude');
        namesCurr(strcmpi(namesCurr,'longitude')) = [];   
    end

    %Remove date/time and attribute fields:
    for jj = numel(namesCurr) : -1 : 1
        if regexpbl(namesCurr{jj}, {'date', 'time', 'attributes', 'att', 'flag', 'lat', 'lon', 'path', 'type'})
            namesCurr(jj) = [];
        end
    end
    
    for jj = 1 : numel(namesCurr)
        switch namesCurr{jj}
%             case 'casi'
%                 fieldObsCurr = {'casi','snw'};
%             case 'stake'
%                 fieldObsCurr = {'icdwe','sndwe'};
%             case 'geodetic'
%                 fieldObsCurr = 'geodetic';
            case 'flowrate'
                fieldObsCurr = 'flow';
            otherwise
                fieldObsCurr = namesCurr(jj);
        end
        
        %Identify fields that currently exist and others that must be
        %created:
        if ~isempty(output)
            if iscell(output) && ~isempty(output(:))
                if iscell(output(:,1)) && numel(output(:,1)) <= 1
                    output(:,1) = output{:,1};
%                 else
%                     error('outputAllGage:outputTypeUnknown',['The output array is type ' ...
%                         class(output(:,1)) ', which has not been programmed for.']);
                end
                
                %Find fields present in both obs and output:
                [fieldOutCurr, ~, indOut] = intersect(fieldObsCurr, output(:,1));
                %Find fields present in obs but not present in output:
                [fieldNotOutCurr, ~] = setdiff(fieldObsCurr, output(:,1));
            else
                fieldOutCurr = '';
                fieldNotOutCurr = {fieldObsCurr};
                output = cell(0,2);
            end
        else
            fieldOutCurr = '';
            fieldNotOutCurr = {fieldObsCurr};
            output = cell(0,2);
        end
       

        %Loop over fields not currently in output: Add fields and
        %output locations:
        if ~isempty(fieldNotOutCurr)
            for kk = 1 : numel(fieldNotOutCurr)
                if isnumeric(lonCurr) && isnumeric(latCurr)
                    output(end+1,1:2) = {char(fieldNotOutCurr{kk}), [lonCurr, latCurr]};
                elseif strcmpi(lonCurr,'avg') || strcmpi(lonCurr,'all') || strcmpi(lonCurr,'writegrid') || strcmpi(lonCurr,'annualgrid') 
                    output(end+1,1:2) = {fieldNotOutCurr{kk}, {lonCurr}};
                end
            end
        end
        
        %Loop over fields currently in obs and output: Ensure coordinates
        %present in both:
        if ~isempty(fieldOutCurr)
            for kk = 1 : numel(fieldOutCurr)
                %This ensures that any lon-lat ordered pairs are rows
                %rather than columns
                if iscell(output{indOut(kk),2}) && numel(output{indOut(kk),2}(1,:)) > numel(output{indOut(kk),2}(:,1))
                    output{indOut(kk),2} = output{indOut(kk),2}';
                end
                
               if isnumeric(output{indOut(kk),2}) && ~ismember([lonCurr,latCurr], output{indOut(kk),2},'rows')
                    if isnumeric(lonCurr) && isnumeric(latCurr)
                        output{indOut(kk),2} = [output{indOut(kk),2}; [lonCurr,latCurr]];
                    elseif strcmpi(lonCurr,'avg')
                        output{indOut(kk),2} = [output{indOut(kk),2}; {'avg'}];
                    elseif strcmpi(lonCurr,'all')
                        output{indOut(kk),2} = [output{indOut(kk),2}; {'all'}];
%                     elseif strcmpi(lonCurr,'writegrid')
%                         output{indOut(kk),2} = [output{indOut(kk),2}; {'writegrid'}];
%                     elseif strcmpi(lonCurr, 'annualgrid')
%                         output{indOut(kk),2} = [output{indOut(kk),2}; {'annualgrid'}];
                    end
               elseif ischar(output{indOut(kk),2})
                    if isnumeric(lonCurr) && isnumeric(latCurr)
                        output{indOut(kk),2} = {output{indOut(kk),2}; [lonCurr,latCurr]};
                    elseif strcmpi(lonCurr,'avg')
                        output{indOut(kk),2} = {output{indOut(kk),2}; 'avg'};
                    elseif strcmpi(lonCurr,'all')
                        output{indOut(kk),2} = [output{indOut(kk),2}; {'all'}];
%                     elseif strcmpi(lonCurr,'writegrid')
%                         output{indOut(kk),2} = [output{indOut(kk),2}; {'writegrid'}];
%                     elseif strcmpi(lonCurr,'annualgrid')
%                         output{indOut(kk),2} = [output{indOut(kk),2}; {'annualgrid'}];
                    end
               elseif iscell(output{indOut(kk),2})
                    if isnumeric(lonCurr) && isnumeric(latCurr)
                        sameCrd = cellfun(@(x) any(ismember([lonCurr,latCurr], x)), output{indOut(kk),2});
                    elseif strcmpi(lonCurr,'avg')
                        sameCrd = cellfun(@(x) any(strcmpi('avg', x)), output{indOut(kk),2});
                    elseif strcmpi(lonCurr,'all')
                        sameCrd = cellfun(@(x) any(strcmpi('all', x)), output{indOut(kk),2});
%                     elseif strcmpi(lonCurr,'writegrid')
%                         sameCrd = cellfun(@(x) any(strcmpi('writegrid', x)), output{indOut(kk),2});
%                     elseif strcmpi(lonCurr,'annualgrid')
%                         sameCrd = cellfun(@(x) any(strcmpi('annualgrid', x)), output{indOut(kk),2});
                    end

                   %Check if point already exists in array:
                   if ~any(sameCrd)
                        if isnumeric(lonCurr) && isnumeric(latCurr)
                            output{indOut(kk),2} = [output{indOut(kk),2}; [lonCurr,latCurr]];
                        elseif strcmpi(lonCurr,'avg')
                            output{indOut(kk),2} = [output{indOut(kk),2}; 'avg'];
                        elseif strcmpi(lonCurr,'all')
                            output{indOut(kk),2} = [output{indOut(kk),2}; 'all'];
%                         elseif strcmpi(lonCurr,'writegrid')
%                             output{indOut(kk),2} = [output{indOut(kk),2}; 'writegrid'];
%                         elseif strcmpi(lonCurr,'annualgrid')
%                             output{indOut(kk),2} = [output{indOut(kk),2}; 'annualgrid'];
                        else
                            error('output_all_gage:unknownCrdMarker', ...
                                ['The longitude marker ' lonCurr ' is not known.']);
                        end
                   end
               end
            end
        end
    end
end 



%%Add column to output cell array that contains linear index 
%(this will be used in funtion 'print_model_state' to speed up processing time):
[ output{:,3} ]= deal(cell(1,1));
for ii = 1 : numel(output(:,1))
    if iscell(output{ii,2})
        output{ii, 3} = cell(numel(output{ii,2}),1);
        
        for jj = 1 : numel(output{ii,2})
            if isnumeric(output{ii,2}{jj})
                [~, gageRow] = min(abs(output{ii,2}{jj}(2)-lat));
                [~, gageCol] = min(abs(output{ii,2}{jj}(1)-lon));
                if output{ii,2}{jj}(2) > lat(1) + 0.5*abs(diff(lat(1:2))) || output{ii,2}{jj}(2) < lat(end) - 0.5*abs(diff(lat(end-1:end))) 
                    warning('SETI_backbone:gagePtRow',['The latitutde of the ' output{ii,1} ' gauge point may be outside the area being modeled.']);
                elseif  output{ii,2}{jj}(1) < lon(1) - 0.5*abs(diff(lon(1:2))) || output{ii,2}{jj}(1) > lon(end) + 0.5*abs(diff(lon(end-1:end))) 
                    warning('SETI_backbone:gagePtCol',['The longitude of the ' output{ii,1} ' gauge point may be outside the area being modeled.']);
                end

                output{ii,3}{jj} = round(sub2ind(szModel, gageRow, gageCol));
            elseif ischar(output{ii,2}{jj})
                if regexpbl(output{ii,2}{jj},{'avg','mean'})
                    output{ii,3}{jj} = 'avg';
                elseif strcmpi(output{ii,2}{jj},{'all'})
                    output{ii,3}{jj} = 'all';
                elseif strcmpi(output{ii,2}{jj},{'writegrid'})
                    output{ii,3}{jj} = 'writegrid';
                elseif strcmpi(output{ii,2}{jj},{'annualgrid'})
                    output{ii,3}{jj} = 'annualgrid';
                end
            end
        end
    elseif isnumeric(output{ii,2})
        output{ii, 3} = cell(numel(output{ii,2}(:,1)),1);

        for jj = 1 : numel(output{ii,2}(:,1))
            [~, gageRow] = min(abs(output{ii,2}(jj,2)-lat));
            [~, gageCol] = min(abs(output{ii,2}(jj,1)-lon));
            if output{ii,2}(jj,2) > lat(1) + 0.5*abs(diff(lat(1:2))) || output{ii,2}(jj,2) < lat(end) - 0.5*abs(diff(lat(end-1:end))) 
                warning('SETI_backbone:gagePtRow','The latitutde of the gauge point may be outside the area being modeled.');
            elseif  output{ii,2}(jj,1) < lon(1) - 0.5*abs(diff(lon(1:2))) || output{ii,2}(jj,1) > lon(end) + 0.5*abs(diff(lon(end-1:end))) 
                warning('SETI_backbone:gagePtCol','The longitude of the gauge point may be outside the area being modeled.');
            end

            output{ii,3}{jj} = ['pt' num2str(round(sub2ind(szModel, gageRow, gageCol)))];
        end
    elseif ischar(output{ii,2})
        if regexpbl(output{ii,2},{'avg','mean'})
            output{ii,3}{1} = 'avg';
        elseif strcmpi(output{ii,2},{'all'})
            output{ii,3}{1} = 'all';
        elseif strcmpi(output{ii,2},{'writegrid'})
            output{ii,3}{1} = 'writegrid';
        elseif strcmpi(output{ii,2},{'annualgrid'})
            output{ii,3}{1} = 'annualgrid';
        end
    else
        error('output_all_gage:unknownOutputType','Unknown type');
    end
end


