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

function varargout = mod_v_obs_v2(sObs, sMod, methRank, varargin)
%Match all observation data to model output for observations taken at one
%point in time or over a duration (must have fields '.[field]_dateStart' and
%'.[field]_dateStart') and calculate performance statistic.  Currently
%all observations of a same type are lumped together before being compared
%(e.g. if many stake measurements present, these will be lumped together
%before calculating a performance statistic).



flagPlot = 0;
flagScatter = 0;
flagGrid = 0;
typeCombine = '';
flagOneStat = 0;
statCombine = 'linear';
outputPath = '';
rmMod = '';
lon = nan;
lat = nan;
wrtGridTyp = [];
if ~isempty(varargin)
    for ii = 1 : numel(varargin) 
        if ischar(varargin{ii})
            switch varargin{ii}
                case 'rmMod'
                    rmMod = 'rmMod';
                case 'plot'
                    flagPlot = 1;

                    [dirCheck,~,~] = fileparts(varargin{ii+1});
                    if isdir(dirCheck)
                        outputPath = varargin{ii+1};
                    else
                        warning('mod_v_obs:dirOutput',['The argument after '...
                            char(39) 'plot' char(39) ' is not an output '...
                            'directory. Therefore no plots will be written '...
                            'to file. Check input argument order']);
                    end
                case 'combineType'
                    typeCombine = 'type';
                case 'lumpType'
                    typeCombine = 'lump';
                case 'oneStat'
                    flagOneStat = 1;
                case 'square'
                    statCombine = 'square';
                case 'scatter'
                   flagScatter = 1; 
                case 'grid'
                    flagGrid = 1;
                case 'lon'
                    if isnumeric(varargin{ii+1})
                        lon = varargin{ii+1};
                    else
                        error('mod_v_obs:lonNotNum','The optional argument for longitude is non-numeric');
                    end
                case 'lat'
                    if isnumeric(varargin{ii+1})
                        lat = varargin{ii+1};
                    else
                        error('mod_v_obs:latNotNum','The optional argument for latitude is non-numeric');
                    end
                case 'wrtGridTyp'
                    wrtGridTyp = varargin{ii+1};
            end
        end
    end
end

% if flagOneStat && typeCombine
%     typeCombine = 0;
%     warning('mod_v_obs:lump',['Both ' char(39) 'lump by type' char(39) ...
%         ' and ' char(39) 'lump all' char(39) ' were entered. ' char(39)...
%         'lump all' char(39) ' will be used.']);
% end

if flagScatter && ~flagPlot
    flagScatter = 0;
    warning('mod_v_obs:scatterNoPlot','The scatter option is selected but no plot is requested.');
end

namesObs = fieldnames(sObs);
namesMod = fieldnames(sMod);



%%FIND ALL DATA FIELDS PRESENT IN OBSERVATION DATA:
obsTypes = cell(0,1);
%Add all fields in observation data
for ii = 1 : numel(namesObs(:)) %Loop over all observation points
    nmDataCurrObs = fieldnames(sObs.(namesObs{ii}));
    obsTypes(end+1:end+numel(nmDataCurrObs),1) = nmDataCurrObs;
end
%Remove all fields dealing with coordinates, time, or attributes:
fieldsRem = {'lat','lon','att','date','time','flag'};
for ii = numel(obsTypes) : -1 : 1
   if regexpbl(obsTypes{ii}, fieldsRem) 
       obsTypes(ii) = [];
   elseif regexpbl(obsTypes{ii}, {'_path', 'path_'})
        obsTypes(ii) = [];
    elseif regexpbl(obsTypes{ii}, {'_type','type_'})
        obsTypes(ii) = [];
   end
end

obsTypes = unique(obsTypes);

%Create output array for observations-model comparisons:
dataAll = cell(0,5);
%Format is: 
    %{ii,1} = {Measurement Type; reference date dateRef; [lon, lat]; eval metric};
    %{ii,2} = dates/days (Model)
    %{ii,3} = data (Model)
    %{ii,4} = days (Obs)
    %{ii,5} = data (Obs)

    
cntrObs = 0;
%Find matches between observation and model output:
for ii = 1 : numel(namesObs(:)) %Loop over all data in the Observation data
    indCurr = strcmpi(namesObs{ii},namesMod);
    
    if all(indCurr == 0)
        warning('mod_vs_obs:missingField',['Observations at ' ...
            namesObs{ii} ' are being skipped because output are not '...
            'available in the model output.']);
        continue
    end
        
    
    fieldsCurrMod = fieldnames(sMod.(namesMod{indCurr}));
    nmDataCurrObs = fieldnames(sObs.(namesObs{ii}));
        %Remove non-data fields
        nmDataCurrObs(strcmpi(nmDataCurrObs,'lat')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'latitude')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'lon')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'longitude')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'att')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'attributes')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'flag')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'date')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'dateStart')) = [];
        nmDataCurrObs(strcmpi(nmDataCurrObs,'dateEnd')) = [];
        %remove any field starting with 'path_' in it:
        for mm = numel(nmDataCurrObs(:)) : -1 : 1
            if regexpbl(nmDataCurrObs{mm}, {'path_', '_path', 'type_', '_type', '_date', 'date_', '_density', 'density_', 'flag_', '_flag'})
                nmDataCurrObs(mm) = [];
            end
        end

    %Loop over fields in observation data (expected that one
    %modelled output exists for each observation field
    for kk = 1 : numel(nmDataCurrObs(:))
        cntrObs = cntrObs + 1; 
        
        %Initialize metadata cell array:
        dataAll{cntrObs,1} = cell(4,1);
    
        %Find indices for observation fields:
        indField = [];
        if regexpbl(nmDataCurrObs{kk},'stake')
            indField = [find(strcmpi(fieldsCurrMod, 'icdwe') == 1); ...
                find(strcmpi(fieldsCurrMod, 'sndwe') == 1)];
            if numel(indField) < 2
                error('mod_v_obs:stakeFields',['Stake data are '...
                    'being used for model assessment.  '...
                    'Either the snowChange or iceChange '...
                    'fields, which are both required, are not '...
                    'present in the model output.']);
            end
            
            evalType = 'sum';
            intrpType = 'sum';
        elseif regexpbl(nmDataCurrObs{kk},{'snow','radar'},'and') || regexpbl(nmDataCurrObs{kk},'swe')
            indField = find(strcmpi(fieldsCurrMod, 'swe') == 1);
            
            evalType = 'sum';
            intrpType = 'sum';
        elseif regexpbl(nmDataCurrObs{kk}, {'casi', 'icx', 'snx', 'sca'})
            %This is tricky because the data to compare with depends on the
            %metric. If metric is ''Parajka' (following Parajka & Blöschl,
            %2008), then compare MODIS SCA with modelled SWE. Otherwise,
            %compare observed SCA with modelled SCA.
            indField = find(strcmpi(fieldsCurrMod, nmDataCurrObs{kk}) == 1);

            if isfield(sObs.(namesObs{ii}),[nmDataCurrObs{kk} '_type'])
                evalType = sObs.(namesObs{ii}).([nmDataCurrObs{kk} '_type']);
            else
                evalType = 'max';
                warning('mod_v_obs:unknownEvalType',['The evaluation '...
                    'type is assumed to be ' evalType ' for aggregating ' ...
                    nmDataCurrObs{kk} ' observations.']);
            end
            
            intrpType = 'mean';
        elseif regexpbl(nmDataCurrObs{kk}, {'geodetic', 'mb', 'massbal'})
            indField = find(strcmpi(fieldsCurrMod, 'geodetic') == 1);
            evalType = 'sum';
            intrpType = 'sum';
        elseif regexpbl(nmDataCurrObs{kk}, {'flow'})
            indField = find(strcmpi(fieldsCurrMod, 'flow') == 1);
            evalType = 'mean';
            intrpType = 'mean';
        elseif regexpbl(fieldsCurrMod, nmDataCurrObs{kk})
            error('modVObs:unknownFieldPresent', ['The gage field ' nmDataCurrObs{kk} ' is present in the modelled data but has not been programmed for. Please add.']);
        else
            error('modVObs:unknownField', ['The gage field ' nmDataCurrObs{kk} ' has not been programmed for.']);
        end
        %Throw error if type not found
        if isempty(indField)
           error('modVObs:noIndField', ['The model indice for gage field ' nmDataCurrObs{kk} ' was not located.']);
        end
        

        %Retrieve dates and temporal resolution of modelled data:
        dataAll{cntrObs,2} = sMod.(namesMod{indCurr}).date; %Modelled date is always .date
        tResMod = numel(dataAll{cntrObs,2}(1,:));
        
        
        %Retrieve data from model:
        if ~isempty(indField) %Modelled field exist for current observation field
            if numel(indField) == 2
                dataAll{cntrObs,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(1)}) ...
                        + sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(2)});
            elseif numel(indField) == 1
                dataAll{cntrObs,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField});
            end
            
            %SPECIAL TREATMENT FOR STAKE DATA:
%             if regexpbl(nmDataCurrObs{kk},'stake') 
%                 if all(sMod.(namesMod{indCurr}).('iceBl')) == 1
%                     dataAll{cntrObs,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(1)}) ...
%                         + sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(2)});
%                 else
%                     dataAll{cntrObs,3} = nan(numel(sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(1)})),1);
%                 end
%             else
%                 dataAll{cntrObs,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField});
%             end
        else
            warning('mod_v_obs:missingModOutput',['Field ' nmDataCurrObs{kk} ...
                ' at observation ' namesObs{ii} ' is being skipped '...
                'because no modelled indices exist at this point.']);
            continue
        end 

        %Retrieve data from observation:
        if  numel(nmDataCurrObs(kk)) == 1
            dataAll{cntrObs,5} = sObs.(namesObs{ii}).(nmDataCurrObs{kk});
            
            %If there is a 'flag' field, value of 1 means don't use data in
            %analysis.
            if isfield(sObs.(namesObs{ii}), 'flag')
                %Different flag treatment for 3D vs. 2D arrays
                if numel(size(dataAll{cntrObs,5})) == 3
                    if numel(sObs.(namesObs{ii}).flag) == numel(dataAll{cntrObs,5}(:,1,1))
                        for mm = 1 : numel(dataAll{cntrObs,5}(:,1,1))
                            if sObs.(namesObs{ii}).flag(mm) == 1
%                                 if any2d(~isnan(squeeze(dataAll{cntrObs,5}(mm,:,:)))) 
%                                     disp(num2str(mm))
%                                 end
                                dataAll{cntrObs,5}(mm,:,:) = nan;
                            end
                        end
                    else
                        error('mod_v_obs:error3DFlagSz','The data quality flag has an unexpected size.');
                    end
                elseif numel(size(dataAll{cntrObs,5})) == 2
                    if numel(sObs.(namesObs{ii}).flag) == numel(dataAll{cntrObs,5})
                        dataAll{cntrObs,5}(sObs.(namesObs{ii}).flag == 1) = nan;
                    else
                        error('mod_v_obs:error2DFlagSz','The data quality flag has an unexpected size.');
                    end
                end
            end
        elseif numel(nmDataCurrObs(kk)) > 1
            error('mod_v_gage:extraObsFields','Extra data fields exist for current observation variable.');
        else
            error('mod_v_gage:noObsFields','No data fields exist for current observation variable.');
        end

        
        if isfield(sObs.(namesObs{ii}), [nmDataCurrObs{kk} '_date']) %Observations occur at single point in time
            dataAll{cntrObs,4} = sObs.(namesObs{ii}).([nmDataCurrObs{kk} '_date']);
            tResObs = numel(dataAll{cntrObs,4}(1,:));

            %If temporal resolution doesn't match, convert (e.g. 
            %if observation data monthly and modelled data daily):
            if tResObs ~= tResMod 
                if tResObs == 3 && tResMod == 2
                    datesAvg = unique(dataAll{cntrObs,4}(:,1:2),'rows');
                    dataAvg = nan(numel(datesAvg(:,1)),1);
                    for ll = 1 : numel(datesAvg(:,1))
                        indAvg = ismember(dataAll{cntrObs,4}(:,1:2), datesAvg(ll,1:2),'rows');
                        dataAvg(ll) = mean(dataAll{cntrObs,5}(indAvg));
                    end

                    dataAll{cntrObs,5} =  dataAvg; %Substitute averages for observations
                    dataAll{cntrObs,4} = datesAvg;
                elseif tResObs == 2 && tResMod == 3
                    datesAvg = unique(dataAll{cntrObs,2}(:,1:2),'rows');
                    dataAvg = nan(numel(datesAvg(:,1)),1);
                    for ll = 1 : numel(datesAvg(:,1))
                        indAvg = ismember(dataAll{cntrObs,2}(:,1:2), datesAvg(ll,1:2),'rows');
                        dataAvg(ll) = mean(dataAll{cntrObs,3}(indAvg));
                    end

                    dataAll{cntrObs,3} = dataAvg; %Substitute averages for modelled
                    dataAll{cntrObs,2} = datesAvg;
                else
                    error('mod_v_obs:tempResDiff',['Each '...
                        'time-series element for the observed '...
                        'data has ' ...
                        num2str(numel(dataAll{cntrObs,4}(1,:))) ...
                        ' componenents while the modelled data '...
                        'has ' num2str(numel(dataAll{cntrObs,2}(1,:)))...
                        '.  This has not been programmed for.']);
                end
            end

            %Find dates common to both modelled and observed data
            [~, dateUseObs, dateUseMod] = intersect(dataAll{cntrObs,4}, dataAll{cntrObs,2},'rows');
            
%             if (sum(~isnan(dateUseObs)) < 0.7*sum(~isnan(dataAll{cntrObs,5}))) || (numel(dateUseObs) < 0.7*numel(dataAll{cntrObs,5}))
%                 warning('mod_v_obs:unusedObs',...
%                     [num2str(sum(~isnan(dateUseObs))) ' of the ' ...
%                     num2str(sum(~isnan(dataAll{cntrObs,5}))) ...
%                     ' observations for ' char(nmDataCurrObs) ' are being '...
%                     'used.  Consider adapting the model '...
%                     'run time accordingly.']);
%             end
            
            dataAll{cntrObs,2} = dataAll{cntrObs,2}( dateUseMod,:);
            dataAll{cntrObs,3} = dataAll{cntrObs,3}( dateUseMod,:);
            dataAll{cntrObs,4} = dataAll{cntrObs,4}( dateUseObs,:);
            dataAll{cntrObs,5} = dataAll{cntrObs,5}( dateUseObs,:);
                       
            if isempty(dateUseObs) || isempty(dateUseMod)
                warning('mod_v_obs:noOverlapObs',...
                    ['There is no overlap between modelled output '...
                    'and observations for ' char(nmDataCurrObs) ...
                    '. Consider adapting the model '...
                    'run time accordingly.']);
                continue
            end

        elseif isfield(sObs.(namesObs{ii}), [nmDataCurrObs{kk} '_dateStart']) && isfield(sObs.(namesObs{ii}), [nmDataCurrObs{kk} '_dateEnd'])
            dataAll{cntrObs,4} = cell(2,1);
            dataAll{cntrObs,4}{1} = sObs.(namesObs{ii}).([nmDataCurrObs{kk} '_dateStart']); 
            dataAll{cntrObs,4}{2} = sObs.(namesObs{ii}).([nmDataCurrObs{kk} '_dateEnd']);
            tResObs = numel(dataAll{cntrObs,4}{1}(1,:));


            %Match observation and modelled dates differently
            %depending on resolutions:             
            if tResObs ~= tResMod 
                if tResObs == 3 && tResMod == 2
                    %Interpolate modelled data to daily resolution
                    %and extract for same dates present in
                    %observations

                    datesTempModInterp = date_vec_fill([dataAll{cntrObs,2}(1,:),1],[dataAll{cntrObs,2}(end,:),eomday(dataAll{cntrObs,2}(end,1), dataAll{cntrObs,2}(end,2))],'gregorian');

                    warning('off', 'all'); %Warning caused if all modelled data are NaN
                    dataAll{cntrObs,3} = interp_month2day(dataAll{cntrObs,2}, dataAll{cntrObs,3}, datesTempModInterp, 'pchip', 1);
                    warning('on', 'all');
                    dataAll{cntrObs,2} = datesTempModInterp;
                else
                    error('mod_v_obs:tempResDiff',['Each '...
                        'time-series element for the observed '...
                        'data has ' ...
                        num2str(numel(dataAll{cntrObs,4}(1,:))) ...
                        ' componenents while the modelled data '...
                        'has ' num2str(numel(dataAll{cntrObs,2}(1,:)))...
                        '.  This has not been programmed '...
                        'for.  The observation data has '...
                        'start and end dates']);
                end
            end
            
            
            %Extract modelled data for time-series present in observations:
            if numel(size(dataAll{cntrObs,5})) == 2 && any(size(dataAll{cntrObs,5}) == 1) %Pt data
                dataTempMod2Obs = nan(numel(dataAll{cntrObs,5}),1);
                for ll = 1 : numel(dataAll{cntrObs,5})
                    [~, indUpMod(1)] = ismember(dataAll{cntrObs,4}{1}(ll,:), dataAll{cntrObs,2}, 'rows');
                    [~, indUpMod(2)] = ismember(dataAll{cntrObs,4}{2}(ll,:), dataAll{cntrObs,2}, 'rows');

                    if all(indUpMod ~= 0) && indUpMod(2) > indUpMod(1) %If start and end dates of observation present in model, copy data
                        if regexpbl(evalType, 'sum')
                            dataTempMod2Obs(ll) = nansum(dataAll{cntrObs,3}(indUpMod(1):indUpMod(2)));
                        elseif regexpbl(evalType, 'mean')
                            dataTempMod2Obs(ll) = nanmean(dataAll{cntrObs,3}(indUpMod(1):indUpMod(2)));
                        elseif regexpbl(evalType, 'max')
                            dataTempMod2Obs(ll) = nanmax(dataAll{cntrObs,3}(indUpMod(1):indUpMod(2)));
                        else
                            error('mod_v_obs:unknownEvalType', [evalType ...
                                ' is an unknown evaluation type and a '...
                                'case has not been programmed for it.']);
                        end
                    else
                        dataAll{cntrObs,5}(ll) = nan;
                    end
                end
            elseif numel(size(dataAll{cntrObs,5})) == 3 || ~any(size(dataAll{cntrObs,5}) == 1)%gridded data
                if numel(size(dataAll{cntrObs,5})) == 3
                    nTsMod = numel(dataAll{cntrObs,5}(:,1,1));
                    dataTempMod2Obs = nan(size(dataAll{cntrObs,5}),'single');             
                else
                    nTsMod = 1;
                    dataTempMod2Obs = nan([1,size(dataAll{cntrObs,5})],'single'); 
                end
                
                for ll = 1 : nTsMod
                    [~, indUpMod(1)] = ismember(dataAll{cntrObs,4}{1}(ll,:), dataAll{cntrObs,2}, 'rows');
                    [~, indUpMod(2)] = ismember(dataAll{cntrObs,4}{2}(ll,:), dataAll{cntrObs,2}, 'rows');

                    %Special treatment for gridded geodetic measurements:
                    if regexpbl(nmDataCurrObs{kk}, {'geodetic', 'ice_change'})
                        %Geodetic observations are gridded and span two points
                        %in time (typically many years). Observation time 
                        %interval may not overlap perfectly with modelled 
                        %time-series. Therefore, formulate units as 
                        %change per 365 days
                        if indUpMod(1) == 0
                            indUpMod(1) = 1;
                        end
                        if indUpMod(2) == 0
                            indUpMod(2) = numel(dataAll{cntrObs,2}(:,1));
                        end
                        
                        nDaysMod = days_since(dataAll{cntrObs,2}(indUpMod(1),:), dataAll{cntrObs,2}(indUpMod(2),:), 'gregorian');
                        dataTempMod2Obs(ll,:,:) = (365/nDaysMod)*squeeze(sum(dataAll{cntrObs,3}(indUpMod(1):indUpMod(2),:,:), 1));
                        
                        nDaysObs = days_since(dataAll{cntrObs,4}{1}(ll,:), dataAll{cntrObs,4}{2}(ll,:), 'gregorian');
                        if numel(size(dataAll{cntrObs,5}(ll,:,:))) == 3
                            dataAll{cntrObs,5}(ll,:,:) = (365/nDaysObs)*squeeze(dataAll{cntrObs,5}(ll,:,:));
                        else
                            mbTemp = dataAll{cntrObs,5}(:,:);
                            dataAll{cntrObs,5} = nan([1, size(mbTemp)]);
                            dataAll{cntrObs,5}(ll,:,:) = (365/nDaysObs)*mbTemp;
                        end
                    elseif all(indUpMod ~= 0) && indUpMod(2) > indUpMod(1) %If start and end dates of observation present in model, copy data
                        if regexpbl(evalType, 'sum')
                            dataTempMod2Obs(ll,:,:) = squeeze(sum(dataAll{cntrObs,3}(indUpMod(1):indUpMod(2),:,:), 1));
                        elseif regexpbl(evalType, 'mean')
                            dataTempMod2Obs(ll,:,:) = squeeze(mean(dataAll{cntrObs,3}(indUpMod(1):indUpMod(2),:,:), 1));
                        elseif regexpbl(evalType, 'max')
                            dataTempMod2Obs(ll,:,:) = squeeze(max(dataAll{cntrObs,3}(indUpMod(1):indUpMod(2),:,:), [], 1));
                        else
                            error('mod_v_obs:unknownEvalType', [evalType ...
                                ' is an unknown evaluation type and a '...
                                'case has not been programmed for it.']);
                        end
                    else
                        dataAll{cntrObs,5}(ll,:,:) = nan;
                    end
                end
            else
                error('mod_v_obs:errDim', ['The current variable has ' ...
                    num2str(numel(size(numel(dataAll{cntrObs,5})))) ...
                    ' dimensions, which has not been programmed for.']);
            end

            %Copy aggregated model data to data array:
            dataAll{cntrObs,3} = dataTempMod2Obs;
            
            %Remove all empty elements for point observation data:
            if numel(size(dataAll{cntrObs,5})) == 2
                indRemCurr = find(isnan(dataAll{cntrObs,5}));
                if ~isempty(indRemCurr)
                    dataAll{cntrObs,5}(indRemCurr) = []; %Observation data field
                    dataAll{cntrObs,3}(indRemCurr) = []; %Modeled data field
                    dataAll{cntrObs,4}{1}(indRemCurr,:) = []; %Observation start date
                    dataAll{cntrObs,4}{2}(indRemCurr,:) = []; %Observation end date
                end
            end
            
            %Set model output dates to be same as observation (because
            %model output has been aggregated to observation time elements)
            dataAll{cntrObs,2} = dataAll{cntrObs,4}; %Model start date
        else
            error('mod_v_obs:noDate',['No time field appears to be present for ' sObs.(namesObs{ii}) ' observation data.'])
        end


        %Save time signature as reference date and "days_since" instead of
        %date:
        dateRef = nan(tResMod, 1);
        if iscell(dataAll{cntrObs,2}) && iscell(dataAll{cntrObs,4})
            if ~isempty(dataAll{cntrObs,2}{1}) && numel(dataAll{cntrObs,2}{1}(1,:)) > 1 && ~isempty(dataAll{cntrObs,4}{1}) && numel(dataAll{cntrObs,4}{1}(1,:)) > 1
                dateRef = dataAll{cntrObs,2}{1}(1,:);
                for ll = 1 : numel(dataAll{cntrObs,2})
                    dataAll{cntrObs,2}{ll} = days_since(dateRef, dataAll{cntrObs,2}{ll}, 'gregorian');
                    dataAll{cntrObs,4}{ll} = days_since(dateRef, dataAll{cntrObs,4}{ll}, 'gregorian');
                end
            end
        elseif iscell(dataAll{cntrObs,2})
            if ~isempty(dataAll{cntrObs,2}{1}) && numel(dataAll{cntrObs,2}{1}(1,:)) > 1
                dateRef = dataAll{cntrObs,2}{1}(1,:);
                for ll = 1 : numel(dataAll{cntrObs,2})
                    dataAll{cntrObs,2}{ll} = days_since(dateRef, dataAll{cntrObs,2}{ll}, 'gregorian');
                end
            end
            dataAll{cntrObs,4} = days_since(dateRef, dataAll{cntrObs,4}, 'gregorian');
        elseif iscell(dataAll{cntrObs,4})
            if ~isempty(dataAll{cntrObs,4}{1}) && numel(dataAll{cntrObs,4}{1}(1,:)) > 1
                dateRef = dataAll{cntrObs,4}{1}(1,:);
                for ll = 1 : numel(dataAll{cntrObs,4})
                    dataAll{cntrObs,4}{ll} = days_since(dateRef, dataAll{cntrObs,4}{ll}, 'gregorian');
                end
            end
            dataAll{cntrObs,2} = days_since(dateRef, dataAll{cntrObs,2}, 'gregorian');
        else
            if ~isempty(dataAll{cntrObs,2}) && numel(dataAll{cntrObs,2}(1,:)) > 1
                dateRef = dataAll{cntrObs,2}(1,:);
                dataAll{cntrObs,2} = days_since(dateRef, dataAll{cntrObs,2}, 'gregorian');
                dataAll{cntrObs,4} = days_since(dateRef, dataAll{cntrObs,4}, 'gregorian');
            end
        end
        
        if all(isnan(dateRef))
           warning('mod_v_obs:noRefDate',...
                    ['The reference date for ' char(nmDataCurrObs) ...
                    ' could not be found. This variable is therefore '...
                    'being skipped.']);
                continue 
        end

        
        %Write metadata for current observation-model:
            % {Measurement Type; reference date dateRef; [lon, lat]; evalulation metric ; evalulation type (sum, max, mean, ...)};
        %This could be defined within the ii loop, but by having it in the
        %kk loop it will be empty if no datafield found.
        dataAll{cntrObs,1}(1:3) = {nmDataCurrObs{kk}; dateRef; [sObs.(namesObs{ii}).lon, sObs.(namesObs{ii}).lat]};
        dataAll{cntrObs,1}{5} = intrpType;
    end    
end


%Remove empty series from combined array:
for kk = numel(dataAll(:,1)) : -1 : 1
    if isempty(dataAll{kk,1}) 
        if numel(size(dataAll{kk,3})) == 3 && (all(all2d(isnan(dataAll{kk,3}))) || all(all2d(isnan(dataAll{kk,5}))))
            dataAll(kk,:) = [];
        elseif all2d(isnan(dataAll{kk,3})) || all2d(isnan(dataAll{kk,5}))
            dataAll(kk,:) = [];
        end
    elseif isempty(dataAll{kk,3}) && isempty(dataAll{kk,5})
        dataAll(kk,:) = [];
    end
end


%Parse evaluation metrics to use:
if ischar(methRank)
   methRank = {methRank};
elseif ~iscell(methRank)
    error('modVObs:unknownRankingClass', ['The fitness metric is of class ' class(methRank) ', which has not been programmed for.'])
end
nMetric = numel(methRank(:));
for ii = 1 : numel(dataAll(:,1))
    
    obsCurr = dataAll{ii,1}{1};
    evalCurr = cell(nMetric, 1);
    [evalCurr{:}] = deal(blanks(0));
    for jj = 1 : nMetric
        %'Parajka' can only be used for snow covered area
        if regexpbl(obsCurr, {'casi', 'icx', 'snx', 'sca'})
            if regexpbl(methRank{jj}, 'Parajka')
                evalCurr{jj} = 'Parajka';
            else
                evalCurr{jj} = methRank{jj};
            end
        elseif regexpbl(obsCurr, {'geodetic', 'mb', 'massbal'})
            if regexpbl(methRank{jj}, 'mbe')
                evalCurr{jj} = 'mbe';
            elseif regexpbl(methRank{jj}, 'mae')
                evalCurr{jj} = 'mae';
            elseif regexpbl(methRank{jj}, 'mape')
                evalCurr{jj} = 'mae';
            else
                evalCurr{jj} = methRank{jj};
            end
        else
            if regexpbl(methRank{jj}, 'Parajka') %Remove 'Parajka' metric reference
                indMeth = regexpi(methRank{jj},'Parajka');
                if indMeth(1) == 1
                    indMeth = [8, numel(methRank{jj})];
                else
                    indMeth = [1, indMeth - 1];
                end
                evalCurr{jj} = methRank{jj}(indMeth(1):indMeth(2));
            else
                evalCurr{jj} = methRank{jj};
            end
        end
        %Remove underscores from start or end of method rank string:
        if regexpbl(evalCurr{jj}(1),'_')
            evalCurr{jj} = evalCurr{jj}(2:end);
        end
        if regexpbl(evalCurr{jj}(end),'_')
            evalCurr{jj} = evalCurr{jj}(1:end-1);
        end
    end
    
    dataAll{ii,1}{4} = evalCurr;
end


%%CALCULATE FITNESS SCORES
nMetric = numel(methRank);
scoreTemp = cell(nMetric, 1);
nmScore   = cell(nMetric, 1);
if regexpbl(typeCombine,'type') %Evaluate individually for each seperate observation and then combined based on unique types of observations (i.e. combine observation performance for same type)
    [scoreTemp{:}] = deal( nan(numel(obsTypes), 1));
    [  nmScore{:}] = deal(cell(numel(obsTypes), 1));
    for kk = 1 : nMetric
        for ii = 1 : numel(obsTypes)
            indCurrType = nan(numel(dataAll(:,1)),1);
            for indCurr = 1 : numel(dataAll(:,1))
                indCurrType(indCurr) = strcmpi(obsTypes{ii}, dataAll{indCurr,1}{1});
            end
            indCurrType = find(indCurrType == 1);

            fitTemp = nan(size(indCurrType));
            for indCurr = 1 : numel(indCurrType)
                if ~isempty(dataAll{indCurrType(indCurr),3}) && ~isempty(dataAll{indCurrType(indCurr),5}) 
                    fitTemp(indCurr) = ...
                        fitness( dataAll{indCurrType(indCurr),5}, dataAll{indCurrType(indCurr),3}, dataAll{indCurrType(indCurr),1}{4}{kk}, rmMod);
                end
            end

            %Combine scored for current type:
            if ~isempty(fitTemp)
                switch statCombine
                    case 'linear'
                        scoreTemp{kk}(ii) = nanmean(fitTemp);
                    case 'square'
                        scoreTemp{kk}(ii) = sqrt(nansum(fitTemp.^2));
                end

                nmScore{kk}{ii} = char(obsTypes{ii});
            end
        end
    end
elseif regexpbl(typeCombine,'lump') %Combine all data of same observation prior to evaluating
    [scoreTemp{:}] = deal( nan(numel(obsTypes), 1));
    [  nmScore{:}] = deal(cell(numel(obsTypes), 1));
    for kk = 1 : nMetric
        for ii = 1 : numel(obsTypes)
            indCurrType = nan(numel(dataAll(:,1)),1);
            for jj = 1 : numel(dataAll(:,1))
                indCurrType(jj) = strcmpi(obsTypes{ii}, dataAll{jj,1}{1});
            end
            indCurrType = find(indCurrType == 1);

            dataFitTempMod = [];
            dataFitTempObs = [];
            for jj = 1 : numel(indCurrType)
                dataFitTempMod = [dataFitTempMod; dataAll{indCurrType(jj),3}];
                dataFitTempObs = [dataFitTempObs; dataAll{indCurrType(jj),5}];
            end

            if ~isempty(dataFitTempObs) && ~isempty(dataFitTempMod) 
                scoreTemp{kk}(ii) = fitness(dataFitTempObs, dataFitTempMod, dataAll{indCurrType(jj),1}{4}{kk}, rmMod);
                nmScore{kk}{ii} = char(obsTypes{ii});
            else
                scoreTemp{kk}(ii) = nan;
            end
        end
    end
else
    [scoreTemp{:}] = deal( nan(numel(dataAll(:,1)),1));
    [  nmScore{:}] = deal(cell(numel(dataAll(:,1)), 1));
    for kk = 1 : nMetric
        for ii = 1 : numel(dataAll(:,1))
            if ~isempty(dataAll{ii,5}) && ~isempty(dataAll{ii,3}) 
                scoreTemp{kk}(ii) = fitness(dataAll{ii,5}, dataAll{ii,3}, dataAll{ii,1}{4}{kk}, rmMod);
                nmScore{kk}{ii} = [char(dataAll{ii,1}{1}), ' (' num2str(round2(dataAll{ii,1}{3}(1),3)) ', ' num2str(round2(dataAll{ii,1}{3}(2),3)) ')'];
            end
        end
    end
end

%Remove empty fields / scores
for kk = 1 : nMetric
    for ii = numel(nmScore{kk}) : -1 : 1
       if isempty(nmScore{kk}{ii})
           nmScore{kk}(ii) = [];
           scoreTemp{kk} = scoreTemp{kk}((1:numel(scoreTemp)) ~= ii); 
       end
    end
end
    
if flagOneStat == 1
    %%COMBINE FITNESS SCORES:
    scoreOut = cell(nMetric, 1);
    for kk = 1 : nMetric
        switch statCombine
            case 'linear'
                scoreOut{kk} = nanmean(scoreTemp{kk});
            case 'square'
                scoreOut{kk} = sqrt(nansum(scoreTemp{kk}.^2));
        end 

        nmTemp = blanks(0);
        for ii = numel(nmScore{kk})
            nmTemp = [nmTemp '-' nmScore{kk}{ii}];
        end

        nmScore{kk} = {'Combined' nmTemp};
   end
else
    scoreOut = scoreTemp;
end


%%PLOT OBSERVATION AND MODEL DATA:
if flagPlot == 1 && numel(dataAll(:,1)) ~= 0
    %%SET PROPERTIES:
    %Set plot size:
    szFrame = [8 5];
    %Set font size:
    ftSz = 14;
        ftSzAx = ftSz - 2;
        ftSzLgd = ftSz - 4;
    %Set line width:
    lnWd = 2;
        lnWdA = lnWd - 0.5;
    %Set font:
    strFont = 'Arial';

    %Use these colors instead?
    %From color Brewer):
    clrBrewer = [228,26,28; ...
        55,126,184; ...
        77,175,74; ...
        152,78,163]/255;
    custGray = [0.5,0.5,0.5];

    for ii = 1 : numel(obsTypes)
        %Guess units:
        switch obsTypes{ii}
            case 'stake'
                unitCurr = 'm';
                varCurr = 'Change in Snow and Ice Water Equivalent';
            case 'snowradar'
                unitCurr = 'm';
                varCurr = 'Snow Water Equivalent';
            case 'flow'
                unitCurr = 'm^3/s';
                varCurr = 'Discharge';
            case 'flowrate'
                unitCurr = 'm^3/s';
                varCurr = 'Discharge';
            case 'sca'
                unitCurr = '% Area';
                varCurr = 'Snow Covered Area';
            case 'casi'
                unitCurr = '% Area';
                varCurr = 'Snow and Ice Covered Area';
            case 'geodetic'
                unitCurr = 'm/yr';
                varCurr = 'Change in Snow and Ice Water Equivalent';
            otherwise
                unitCurr = '?';
                varCurr = obsTypes{ii};
        end

        if ~isempty(typeCombine)
            indCurrType = nan(numel(dataAll(:,1)),1);
            for indCurr = 1 : numel(dataAll(:,1))
                indCurrType(indCurr) = strcmpi(obsTypes{ii}, dataAll{indCurr,1}{1});
            end
            if all(indCurrType == 0)
               continue 
            end
            indCurrType = {find(indCurrType == 1)};
            nPlot = 1;
        else
            nPlot = numel(dataAll(:,1));
            indCurrType = cell(nPlot,1);
            for ll = 1 : nPlot
                indCurrType{ll} = ll;
            end
        end
        

        %Loop over each plot:
        for ll = 1 : nPlot
            szCurr = size(dataAll{indCurrType{ll}(1),5});
            if numel(dataAll{indCurrType{ll}(1),2}) == 2
                nDtCurr = numel(dataAll{indCurrType{ll}(1),2}{1}(:,1));
                dateTyp = 2;
            else
                nDtCurr = numel(dataAll{indCurrType{ll}(1),2}(:,1));
                dateTyp = 1;
            end

            %Make gridded plots if data are gridded
            if ((szCurr(1) == nDtCurr && numel(szCurr) == 3) || (szCurr(1) ~= nDtCurr && numel(szCurr) == 2)) && flagGrid
                if isdir(outputPath)
                    dirWrite = outputPath;
                    fileWrite = '';
                else
                    [dirWrite, fileWrite, ~] = fileparts(outputPath);
                end
                
                %Use imagesc to plot data:
                custClrMap = ... 
                    [0,0,1; ... %blue
                    0.2,0.2,1; ...
                    0.4,0.4,1; ...
                    0.6,0.6,1; ...
                    0.8,0.8,1; ...
                    1, 1, 1; ... %white
                    1, 0.8, 0.8; ...
                    1, 0.6, 0.6; ...
                    1, 0.4, 0.4; ...
                    1, 0.2, 0.2; ...
                    1, 0, 0]; %red
%                     0.7, 0.7, 0.7]; %Gray 
                
                %Initialize arrays:
                gridObsCurr = nan([nDtCurr, szCurr(end-1:end)], 'single');
                gridModCurr = nan([nDtCurr, szCurr(end-1:end)], 'single');
                dataErrCurr = nan([nDtCurr, szCurr(end-1:end)], 'single');
                dataMapeCurr = nan([nDtCurr, szCurr(end-1:end)], 'single');
                
                %loop over all iterations of the observations
                for kk = 1 : nDtCurr
                    %Format data for plotting (Mod - Obs):
                    %Set values to nan where change in observation is exactly 0: 
                    if numel(szCurr) == 3
                        gridObsCurr(kk,:,:) = squeeze(dataAll{indCurrType{ll}(1),5}(kk,:,:));
                        gridModCurr(kk,:,:) = squeeze(dataAll{indCurrType{ll}(1),3}(kk,:,:));
                        dataErrCurr(kk,:,:) = gridModCurr(kk,:,:) - gridObsCurr(kk,:,:);
                        
                        dataMapeTemp = 100*squeeze((dataAll{indCurrType{ll}(1),3}(kk,:,:) - gridObsCurr(kk,:,:))./gridObsCurr(kk,:,:));
                        dataMapeTemp(isnan(squeeze(gridObsCurr(kk,:,:)))) = nan;
                        dataMapeTemp(isnan(squeeze(dataAll{indCurrType{ll}(1),3}(kk,:,:)))) = nan;
                        
                        dataMapeCurr(kk,:,:) = dataMapeTemp;
                        
%                         sz2d = size(gridObsCurr);
                    else
                        gridObsCurr(kk,:,:) = dataAll{indCurrType{ll}(1),5};
                        gridModCurr(kk,:,:) = dataAll{indCurrType{ll}(1),3};
                        dataErrCurr(kk,:,:) = gridModCurr - gridObsCurr;
                        dataMapeCurr(kk,:,:) = 100*(gridModCurr - gridObsCurr)./gridObsCurr;
                       
%                         sz2d = size(dataAll{indCurrType{ll}(1),5});
                    end
%                     %Set to nan:
%                     if strcmpi(obsTypes{ii}, 'geodetic')
%                         indNan2d = find(gridObsCurr == 0);
%                         [nanRow, nanCol] = ind2sub(sz2d, indNan2d);
%                         sz3d = size(dataErrCurr);
%                         indNan3d = kk + (nanRow-1)*sz3d(1) + (nanCol-1)*sz3d(2)*sz3d(1);
%                         dataErrCurr(indNan3d) = nan;
%                         datMapeCurr(indNan3d) = nan;
%                     end
                    
                    if nDtCurr <= 10 %Only make output figures if there are 10 or fewer observations
                        %Plot observed grid:
                        hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'Color',[1 1 1]);

                        %Remove any "nan border" caused by different
                        %geographic domains:
                        [gridObsPlot, latPlot, lonPlot] = rm_nan_border(squeeze(gridObsCurr(kk,:,:)), lat, lon);
                        [gridModPlot, latPlot, lonPlot] = rm_nan_border(squeeze(gridModCurr(kk,:,:)), lat, lon);
                        
                        cLim = max([max2d(abs(gridModPlot)), max2d(abs(gridObsPlot))]);
                        
%                         gridObsCurr(isnan(gridObsCurr)) = cLim + 0.1*cLim;
                        if all(~isnan(lon)) && all(~isnan(lat))
                            imagesc(lonPlot, latPlot, gridObsPlot);
                            xlim([min(lonPlot), max(lonPlot)]);
                            ylim([min(latPlot), max(latPlot)]);
                        else   
                            imagesc(gridObsPlot);
                        end
                        set(gca,'Color',[0.8 0.8 0.8]);
                        shading flat; 
                        colormap(custClrMap);
                        caxis([-cLim cLim]);
                        hCBar = colorbar;  
                        hTitle = title(['Observed ' varCurr ' (' unitCurr ')']);

                        %Set plot properties:
                        set(gca, ...
                            'FontName'   , strFont);
                        set([hTitle, hCBar], ...
                            'FontName'   , strFont, ...
                            'fontSize', ftSzAx);
                        set(gca, ...
                            'Box'         , 'on', ...
                            'TickDir'     , 'in'     , ...
                            'TickLength'  , [.02 .02] , ...
                            'LineWidth'   , lnWdA , ...
                            'fontSize', ftSzAx);   

                        if ~isempty(fileWrite)
                            pathWrite = fullfile(dirWrite, [fileWrite '_' num2str(ii)]);
                        else
                            pathWrite = fullfile(dirWrite, [obsTypes{ii} '_obs_grid']);
                        end

                        if nDtCurr > 1
                            pathWrite = [char(pathWrite), '_' num2str(kk)];
                        end

                        savefig(hFig, [pathWrite '.fig']);
                        print(hFig, [pathWrite '.eps'],'-depsc2');
                        print(hFig, [pathWrite '.png'],'-dpng','-r600');
                        
                        
                        %Plot modelled grid:
                        hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'Color',[1 1 1]);

%                         gridModCurr(isnan(gridModCurr)) = cLim + 0.1*cLim;
                        if all(~isnan(lon)) && all(~isnan(lat))
                            imagesc(lonPlot, latPlot, gridModPlot);
                            xlim([min(lonPlot), max(lonPlot)]);
                            ylim([min(latPlot), max(latPlot)]);
                        else   
                            imagesc(gridModPlot);
                        end
                        set(gca,'Color',[0.8 0.8 0.8]);
                        shading flat; 
                        colormap(custClrMap);
                        caxis([-cLim cLim]);
                        hCBar = colorbar;  
                        hTitle = title(['Modelled ' varCurr ' (' unitCurr ')']);

                        %Set plot properties:
                        set(gca, ...
                            'FontName'   , strFont);
                        set([hTitle, hCBar], ...
                            'FontName'   , strFont, ...
                            'fontSize', ftSzAx);
                        set(gca, ...
                            'Box'         , 'on', ...
                            'TickDir'     , 'in'     , ...
                            'TickLength'  , [.02 .02] , ...
                            'LineWidth'   , lnWdA , ...
                            'fontSize', ftSzAx);   

                        if ~isempty(fileWrite)
                            pathWrite = fullfile(dirWrite, [fileWrite '_' num2str(ii)]);
                        else
                            pathWrite = fullfile(dirWrite, [obsTypes{ii} '_mod_grid']);
                        end

                        if nDtCurr > 1
                            pathWrite = [char(pathWrite), '_' num2str(kk)];
                        end

                        savefig(hFig, [pathWrite '.fig']);
                        print(hFig, [pathWrite '.eps'],'-depsc2');
                        print(hFig, [pathWrite '.png'],'-dpng','-r600');
                        
                        
                        %PLOT LEVEL ERROR
                        hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'Color',[1 1 1]);

                        %Remove any "nan border" caused by different
                        %geographic domains:
                        [dataPlot, latPlot, lonPlot] = rm_nan_border(squeeze(dataErrCurr(kk,:,:)), lat, lon);
                        
                        cLim = max2d(abs(dataPlot));
%                         dataPlot(isnan(dataPlot)) = cLim + 0.1*cLim;
                        if all(~isnan(lon)) && all(~isnan(lat))
                            imagesc(lonPlot, latPlot, dataPlot);
                            xlim([min(lonPlot), max(lonPlot)]);
                            ylim([min(latPlot), max(latPlot)]);
                        else   
                            imagesc(dataPlot);
                        end
                        set(gca,'Color',[0.8 0.8 0.8]);
                        shading flat; 
                        colormap(custClrMap);
                        caxis([-cLim cLim]);
                        hCBar = colorbar;  
                        hTitle = title(['Modelled Minus Observed ' varCurr ' (' unitCurr ')']);

                        %Set plot properties:
                        set(gca, ...
                            'FontName'   , strFont);
                        set([hTitle, hCBar], ...
                            'FontName'   , strFont, ...
                            'fontSize', ftSzAx);
                        set(gca, ...
                            'Box'         , 'on', ...
                            'TickDir'     , 'in'     , ...
                            'TickLength'  , [.02 .02] , ...
                            'LineWidth'   , lnWdA , ...
                            'fontSize', ftSzAx);   

                        if ~isempty(fileWrite)
                            pathWrite = fullfile(dirWrite, [fileWrite '_' num2str(ii)]);
                        else
                            pathWrite = fullfile(dirWrite, [obsTypes{ii} '_mod_v_obs_lvl_err']);
                        end

                        if nDtCurr > 1
                            pathWrite = [char(pathWrite), '_' num2str(kk)];
                        end
                        
                        savefig(hFig, [pathWrite '.fig']);
                        print(hFig, [pathWrite '.eps'],'-depsc2');
                        print(hFig, [pathWrite '.png'],'-dpng','-r600');
                        

                        %PLOT PERCENT ERROR:
                        hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'Color',[1 1 1]);

                        %Remove any "nan border" caused by different
                        %geographic domains:
                        [dataPlot, latPlot, lonPlot] = rm_nan_border(squeeze(dataMapeCurr(kk,:,:)), lat, lon);
                        dataPlot(isinf(dataPlot)) = nan;
                        
                        
                        cLim = max2d(abs(dataPlot));
                        cLim(cLim > 200) = 100;
                        
%                         dataPlot(isnan(dataPlot)) = cLim + 0.1*cLim;
                        if all(~isnan(lon)) && all(~isnan(lat))
                            imagesc(lonPlot, latPlot, dataPlot);
                            xlim([min(lonPlot), max(lonPlot)]);
                            ylim([min(latPlot), max(latPlot)]);
                        else   
                            imagesc(dataPlot);
                        end
                        set(gca,'Color',[0.8 0.8 0.8]);
                        shading flat; 
                        colormap(custClrMap);
                        caxis(gca, [-cLim cLim]);
                        hCBar = colorbar(gca);  
                        hTitle = title(['Percent Error ' varCurr ' (%; thresh = 100%)']);

                        %Set plot properties:
                        set(gca, ...
                            'FontName'   , strFont);
                        set([hTitle, hCBar], ...
                            'FontName'   , strFont, ...
                            'fontSize', ftSzAx);
                        set(gca, ...
                            'Box'         , 'on', ...
                            'TickDir'     , 'in'     , ...
                            'TickLength'  , [.02 .02] , ...
                            'LineWidth'   , lnWdA , ...
                            'fontSize', ftSzAx);   

                        if ~isempty(fileWrite)
                            pathWrite = fullfile(dirWrite, [fileWrite '_' num2str(ii)]);
                        else
                            pathWrite = fullfile(dirWrite, [obsTypes{ii} '_mod_v_obs_per_err']);
                        end

                        if nDtCurr > 1
                            pathWrite = [char(pathWrite), '_' num2str(kk)];
                        end

                        savefig(hFig, [pathWrite '.fig']);
                        print(hFig, [pathWrite '.eps'],'-depsc2');
                        print(hFig, [pathWrite '.png'],'-dpng','-r600');
                        
                    end
                end %End of date loop
                

                %%Write evaluation grids to file:
                dateRef = dataAll{indCurrType{ll}(1),1}{2};
                %Get dates and Estimate time bounds:
                if  dateTyp == 2
                    daysBndsOut = [dataAll{indCurrType{ll}(1),2}{1}(:), dataAll{indCurrType{ll}(1),2}{2}(:)];
                    
                    daysUse = nanmean(daysBndsOut, 2);
                    dateOut = days_2_date(daysUse, dateRef, 'gregorian');
                    dateWrtOut = [{days_2_date(daysBndsOut(:,1), dateRef, 'gregorian')}, {days_2_date(daysBndsOut(:,2), dateRef, 'gregorian')}];
                elseif dateTyp == 1
                    daysOut = dataAll{indCurrType{ll}(1),2};
                    
                    if nDtCurr > 2
                        dayStepTemp = abs(diff(daysOut(:,1)));
                        dayStep = mode(dayStepTemp);
                        if dayStep ~= dayStepTemp(1)
                           warning('mod_v_obs:dayStepWrong',['The  time step for ' ...
                               obsTypes{ii} ' was estimated as ' num2str(dayStep) ...
                               ' but this appears to be wrong.']); 
                        end
                        daysBndsOut = [daysOut - 0.5*dayStep; daysOut(end,:) + 0.5*dayStep];
                    else
                        daysBndsOut = [dataAll{indCurrType{ll}(1),2}{1}(1,:); dataAll{indCurrType{ll}(1),2}{2}(1,:)];
                    end
                    
                    dateOut = days_2_date(daysOut, dateRef, 'gregorian');
                    dateWrtOut = {dateOut};
                else
                    error('modVObs:unknownDateType',['The date type with ' num2str(dateTyp) ' has not been programmed for.'])
                end
                
                %Set parameters for writing files:
                nDec = 2;
                hdrOut = ESRI_hdr(lon, lat, 'corner');
                
                if regexpbl(wrtGridTyp, {'nc', 'netcdf'})
                    extWrt = 'nc';
                    
                    dateBndsOut = days_2_date(daysBndsOut, dateRef, 'gregorian');
                elseif regexpbl(wrtGridTyp, 'asc')
                    extWrt = 'asc';
                else
                    error('modVObs:unknownWrtTyp', ['The write type ' wrtGridTyp ' has not been programmed for.'])
                end

                %Write model output:
                pathOut = fullfile(dirWrite, [obsTypes{ii} '_grids'], [obsTypes{ii} '_mod.' extWrt]);
                [foldCurr, fileCurr, extCurr] = fileparts(pathOut);
                if ~isempty(wrtGridTyp)
                    if regexpbl(wrtGridTyp, {'nc', 'netcdf'})
                        if exist(pathOut, 'file')
                           delete(pathOut); 
                        end
                        print_grid_NC_v2(pathOut, gridModCurr, obsTypes{ii}, lon, lat, dateOut, dateBndsOut, unitCurr);
                    elseif regexpbl(wrtGridTyp, 'asc')
                        fileNmCurr = file_nm(fileCurr, 'ts', dateWrtOut);
                        
                        for kk = 1 : nDtCurr
                            pathOutCurr = fullfile(foldCurr, [fileNmCurr{kk}, extCurr]);
                            if exist(pathOutCurr, 'file')
                               delete(pathOutCurr); 
                            end
                            write_ESRI_v4(squeeze(gridModCurr(kk,:,:)), hdrOut, pathOutCurr, nDec); 
                        end
                    end
                end
                %Write time-cumulative grid:
                evalCurr = dataAll{indCurrType{ll},1}{5};
                switch evalCurr
                    case 'mean'
                        gridWrt = squeeze(nanmean(gridModCurr, 1));
                    case 'sum'
                        gridWrt = squeeze(nansum(gridModCurr, 1));
                    case 'min'
                        gridWrt = squeeze(nanmin(gridModCurr, 1));
                    case 'max'
                        gridWrt = squeeze(nanmax(gridModCurr, 1));
                    otherwise
                        error('modVObs:unknownEvalType', ['The evaluation type ' dataAll{indCurrType{ll},1}{5} ' is not known.']);
                end
                pathOutAvg = fullfile(dirWrite, [fileCurr, '_' evalCurr, '.asc']);
                write_ESRI_v4(gridWrt, hdrOut, pathOutAvg, nDec);
                
                %Write observed output:
                pathOut = fullfile(dirWrite, [obsTypes{ii} '_grids'], [obsTypes{ii} '_obs.' extWrt]);
                [foldCurr, fileCurr, extCurr] = fileparts(pathOut);
                if ~isempty(wrtGridTyp)
                    if regexpbl(wrtGridTyp, {'nc', 'netcdf'})
                        if exist(pathOut, 'file')
                           delete(pathOut); 
                        end
                        print_grid_NC_v2(pathOut, gridObsCurr, obsTypes{ii}, lon, lat, dateOut, dateBndsOut, unitCurr);
                    elseif regexpbl(wrtGridTyp, 'asc')
                        fileNmCurr = file_nm(fileCurr, 'ts', dateWrtOut);
                        
                        for kk = 1 : nDtCurr
                            pathOutCurr = fullfile(foldCurr, [fileNmCurr{kk}, extCurr]);
                            if exist(pathOutCurr, 'file')
                               delete(pathOutCurr); 
                            end
                            write_ESRI_v4(squeeze(gridObsCurr(kk,:,:)), hdrOut, pathOutCurr, nDec); 
                        end
                    end
                end
                %Write time-cumulative grid:
                switch evalCurr
                    case 'mean'
                        gridWrt = squeeze(nanmean(gridObsCurr, 1));
                    case 'sum'
                        gridWrt = squeeze(nansum(gridObsCurr, 1));
                    case 'min'
                        gridWrt = squeeze(nanmin(gridObsCurr, 1));
                    case 'max'
                        gridWrt = squeeze(nanmax(gridObsCurr, 1));
                    otherwise
                        error('modVObs:unknownEvalType', ['The evaluation type ' dataAll{indCurrType{ll},1}{5} ' is not known.']);
                end
                pathOutAvg = fullfile(dirWrite, [fileCurr, '_' evalCurr, '.asc']);
                write_ESRI_v4(gridWrt, hdrOut, pathOutAvg, nDec);
                
                %Write level error:
                pathOut = fullfile(dirWrite, [obsTypes{ii} '_grids'], [obsTypes{ii} '_mod_minus_obs.' extWrt]);
                [foldCurr, fileCurr, extCurr] = fileparts(pathOut);
                if ~isempty(wrtGridTyp)
                    if regexpbl(wrtGridTyp, {'nc', 'netcdf'})
                        if exist(pathOut, 'file')
                           delete(pathOut); 
                        end
                        print_grid_NC_v2(pathOut, dataErrCurr, obsTypes{ii}, lon, lat, dateOut, dateBndsOut, unitCurr);
                    elseif regexpbl(wrtGridTyp, 'asc')
                        fileNmCurr = file_nm(fileCurr, 'ts', dateWrtOut);
                        
                        for kk = 1 : nDtCurr
                            pathOutCurr = fullfile(foldCurr, [fileNmCurr{kk}, extCurr]);
                            if exist(pathOutCurr, 'file')
                               delete(pathOutCurr); 
                            end
                            write_ESRI_v4(squeeze(dataErrCurr(kk,:,:)), hdrOut, pathOutCurr, nDec); 
                        end
                    end
                end
                %Write average:
                pathOutAvg = fullfile(dirWrite, [fileCurr, '_avg', '.asc']);
                write_ESRI_v4(squeeze(nanmean(dataErrCurr, 1)), hdrOut, pathOutAvg, nDec);
                
                %Write Percent Error:
                pathOut = fullfile(dirWrite, [obsTypes{ii} '_grids'], [obsTypes{ii} '_mod_percent_error.' extWrt]);
                [foldCurr, fileCurr, extCurr] = fileparts(pathOut);
                if ~isempty(wrtGridTyp)
                    if regexpbl(wrtGridTyp, {'nc', 'netcdf'})
                        if exist(pathOut, 'file')
                           delete(pathOut); 
                        end
                        print_grid_NC_v2(pathOut, dataMapeCurr, obsTypes{ii}, lon, lat, dateOut, dateBndsOut, unitCurr);
                    elseif regexpbl(wrtGridTyp, 'asc')
                        fileNmCurr = file_nm(fileCurr, 'ts', dateWrtOut);
                                          
                        for kk = 1 : nDtCurr
                            pathOutCurr = fullfile(foldCurr, [fileNmCurr{kk}, extCurr]);
                            if exist(pathOutCurr, 'file')
                               delete(pathOutCurr); 
                            end
                            write_ESRI_v4(squeeze(dataMapeCurr(kk,:,:)), hdrOut, pathOutCurr, nDec); 
                        end
                    end
                end
                %Write average:
                pathOutAvg = fullfile(dirWrite, [fileCurr, '_avg', '.asc']);
                write_ESRI_v4(squeeze(nanmean(dataMapeCurr, 1)), hdrOut, pathOutAvg, nDec);
            end %End of grid analysis
            

            %Time-series/Scatter plots (even of gridded data)
            nSeries = numel(indCurrType{ll});

            hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame]);
            hold on
            if flagScatter
                hTs = nan(nSeries,1);
                xMin = 0;
                xMax = 0;
                strSeries = cell(nSeries,1);
            else
                hTs = nan(2*nSeries,1);
                xMin = 10000;
                xMax = 0;
                strSeries = cell(2*nSeries,1);
            end

            %Set colors:
            if nSeries <= 4
                colorsUse = clrBrewer;
            else
                colorsUse = distinguishable_colors( nSeries );
            end

            %Plot each series for the given data type:
            for indCurr = 1 : nSeries
                if flagScatter %Display data as scatter plot, otherwise make time-series
                    if numel(size(dataAll{indCurrType{ll}(indCurr),3})) == 2
                        plotModData  = dataAll{indCurrType{ll}(indCurr),3};
                    elseif numel(size(dataAll{indCurrType{ll}(indCurr),3})) == 3
                        plotModData = nan(numel(dataAll{indCurrType{ll}(indCurr),3}(:,1,1)),1);
                        for mm = 1 : numel(dataAll{indCurrType{ll}(indCurr),3}(:,1,1))
                            plotModData(mm)  = mean2d(squeeze(dataAll{indCurrType{ll}(indCurr),3}(mm,:,:)));
                        end
                    end
                    if numel(size(dataAll{indCurrType{ll}(indCurr),5})) == 2
                        plotObsData  = dataAll{indCurrType{ll}(indCurr),5};
                    elseif numel(size(dataAll{indCurrType{ll}(indCurr),5})) == 3
                        plotObsData = nan(numel(dataAll{indCurrType{ll}(indCurr),5}(:,1,1)),1);
                        for mm = 1 : numel(dataAll{indCurrType{ll}(indCurr),5}(:,1,1))
                            plotObsData(mm)  = mean2d(squeeze(dataAll{indCurrType{ll}(indCurr),5}(mm,:,:)));
                        end
                    end

                    hTs(indCurr) = scatter(plotObsData,plotModData, ...
                        2*ftSz, 'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',colorsUse(mod(indCurr-1, nSeries)+1,:));
                    xMin = min(min(real([dataAll{indCurrType{ll}(indCurr),3}(:); dataAll{indCurrType{ll}(indCurr),5}(:)])), xMin);
                    xMax = max(max(real([dataAll{indCurrType{ll}(indCurr),3}(:); dataAll{indCurrType{ll}(indCurr),5}(:)])), xMax);

                    %Create legend entries:
                    lgdCurr = dataAll{indCurrType{ll}(indCurr),1}{3};
                    if ischar(lgdCurr)
                        if regexpbl(lgdCurr, {'all', 'avg'})
                            strSeries{indCurr} = 'Grid Mean';
                        else
                            error('modVObs:gridLabel', [lgdCurr ' is not a recognized observation location labeling option.']);
                        end
                    elseif isnumeric(lgdCurr)
                        strSeries{indCurr} = ['Pt (' num2str(round2(lgdCurr(1),3)) ', ' num2str(round2(lgdCurr(2),3)) ')'];
                    else
                        error('modVObs:gridLabelClass', [class(lgdCurr) ' is not a recognized format for the observation location label.']);
                    end
                    
                else %time-series without linear time-step increases (such as glacier stakes)
                    cntrTs = 2*(indCurr-1) + 1 : 2*indCurr;

                    if iscell(dataAll{indCurrType{ll}(indCurr),2})
                        plotDates = [];
                        plotModData  = [];
                        plotObsData  = [];
                        for kk = 1 : numel(dataAll{indCurrType{ll}(indCurr),2}{1})
                            tempModDates = (dataAll{indCurrType{ll}(indCurr),2}{1}(kk) : dataAll{indCurrType{ll}(indCurr),2}{2}(kk));
                            plotDates = [plotDates, tempModDates];
                            if numel(size(dataAll{indCurrType{ll}(indCurr),3})) == 2
                                plotModData  = [plotModData, dataAll{indCurrType{ll}(indCurr),3}(kk)*ones(1, numel(tempModDates))];
                            elseif numel(size(dataAll{indCurrType{ll}(indCurr),3})) == 3
                                plotModData  = [plotModData, mean2d(squeeze(dataAll{indCurrType{ll}(indCurr),3}(kk,:,:)))*ones(1, numel(tempModDates))];
                            end
                            if numel(size(dataAll{indCurrType{ll}(indCurr),5})) == 2
                                plotObsData  = [plotObsData, dataAll{indCurrType{ll}(indCurr),5}(kk)*ones(1, numel(tempModDates))];
                            elseif numel(size(dataAll{indCurrType{ll}(indCurr),5})) == 3
                                plotObsData  = [plotObsData, mean2d(squeeze(dataAll{indCurrType{ll}(indCurr),5}(kk,:,:)))*ones(1, numel(tempModDates))];
                            end

                        end

                        hTs(cntrTs)= plot(plotDates, plotModData, 'x', plotDates, plotObsData, 'o'); 

                        %Find min and max dates:
                        xMin = min(min(plotDates), xMin);
                        xMax = max(max(plotDates), xMax);
                    else %Regular time-series
                        hTs(cntrTs)= plot(dataAll{indCurrType{ll}(indCurr),2}, dataAll{indCurrType{ll}(indCurr),3}, '-x', dataAll{indCurrType{ll}(indCurr),4}, dataAll{indCurrType{ll}(indCurr),5}, '-o'); 

                        %Find min and max dates:
                        xMin = min(dataAll{indCurrType{ll}(indCurr),2}(1), xMin);
                        xMax = max(dataAll{indCurrType{ll}(indCurr),2}(end), xMax);
                    end

                    %Set properties of current series:
                    set(hTs(cntrTs(1)),'Color',([0.5,0.5,0.5] + 0.5*colorsUse(mod(indCurr-1, nSeries)+1,:)),'LineWidth',lnWdA);
                    set(hTs(cntrTs(2)), 'Color',colorsUse(mod(indCurr-1, nSeries)+1,:));

                    %Create legend entries:
                    lgdCurr = dataAll{indCurrType{ll}(indCurr),1}{3};
                    if ischar(lgdCurr)
                        if regexpbl(lgdCurr, {'all', 'avg'})
                            strCurrPt = 'Grid Mean';
                        else
                            error('modVObs:gridLabel', [lgdCurr ' is not a recognized observation location labeling option.']);
                        end
                    elseif isnumeric(lgdCurr)
                        strCurrPt = ['Pt (' num2str(round2(lgdCurr(1),3)) ', ' num2str(round2(lgdCurr(2),3)) ')'];
                    else
                        error('modVObs:gridLabelClass', [class(lgdCurr) ' is not a recognized format for the observation location label.']);
                    end
                    
                    %strCurrPt = ['Pt (' num2str(round2(dataAll{indCurrType{ll}(indCurr),1}{3}(1),3)) ', ' num2str(round2(dataAll{indCurrType{ll}(indCurr),1}{3}(2),3)) ')'];
                    strSeries(cntrTs(1):cntrTs(2)) = {['Mod @ ' strCurrPt]; ['Obs @ ' strCurrPt]};
                end
            end

            if xMax == Inf || (exist('yMax','var') && yMax == Inf)
                h = gcf;
                axesObjs = get(h, 'Children');  %axes handles
                dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

                xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
                ydata = get(dataObjs, 'YData');

                xMax = 0;
                for zz = 1 : numel(xdata)
                    indNInf = find(xdata{zz} ~= Inf);
                    if ~isempty(indNInf)
                       xMax = max([xMax,xdata{zz}(indNInf)]); 
                    end
                end

                yMax = 0;
                for zz = 1 : numel(ydata)
                    indNInf = find(ydata{zz} ~= Inf);
                    if ~isempty(indNInf)
                       yMax = max([yMax,ydata{zz}(indNInf)]); 
                    end
                end
            end

            if flagScatter
                xlim(real([xMin,xMax]));
                ylim(real([xMin,xMax]));

                %Add reference line(s)
                hRef(1) = line([xMin xMax],[xMin xMax]);

                %Add text tag:
                yMin = ylim;
                    yMax = yMin(2);
                    yMin = yMin(1);

                if xMin < 0
                    hRef(2) = line([xMin xMax],[0 0]);
                    hRef(3) = line([0 0], [yMin yMax]);

                    hText(1) = text(double(0.015*(xMax-xMin) + xMin), double(0.95*(yMax-yMin) + yMin),'Model Overpredicts');
                    hText(2) = text(double( 0.73*(xMax-xMin) + xMin), double(0.68*(yMax-yMin) + yMin),'Model Underpredicts');

                    hText(3) = text(double(0.015*(xMax-xMin) + xMin), double(0.32*(yMax-yMin) + yMin),'Model Underpredicts');
                    hText(4) = text(double( 0.73*(xMax-xMin) + xMin), double(0.07*(yMax-yMin) + yMin),'Model Overpredicts');
                else
                    hText(1) = text(double(0.05*(xMax-xMin) + xMin), double(0.95*(yMax-yMin) + yMin),'Model Overpredicts');
                    hText(2) = text(double( 0.7*(xMax-xMin) + xMin), double(0.07*(yMax-yMin) + yMin),'Model Underpredicts');
                end


                %Set Labels:
                hXLab = xlabel(['Observed (' unitCurr ')']);
                hYLab = ylabel(['Modelled (' unitCurr ')']);

                %Set labels
                if xMax - xMin + 1 <= 10
                    ptsTicks = (xMin:1:xMax);

%                     nXTicks = xMax - xMin + 1;
%                     ptsTicks =  round(linspace(xMin,xMax,nXTicks));
                elseif xMax - xMin + 1 <= 20                    
                    ptsTicks = (xMin:2:xMax);

%                     nXTicks = xMax - xMin + 1;
%                     ptsTicks =  round(linspace(xMin,xMax,nXTicks));
                elseif xMax - xMin + 1 <= 40
                    ptsTicks = (xMin:4:xMax);
                elseif xMax - xMin + 1 <= 60
                    ptsTicks = (xMin:6:xMax);
                elseif xMax - xMin + 1 <= 80
                    ptsTicks = (xMin:8:xMax);
                elseif xMax - xMin + 1 <= 100
                    ptsTicks = (xMin:10:xMax);
                else
                    nXTicks = 12;
                    ptsTicks =  round(linspace(xMin,xMax,nXTicks));
                end

                set(gca, 'XTickLabel', round2(ptsTicks,1), 'XTick', ptsTicks);
                set(gca, 'YTickLabel', round2(ptsTicks,1), 'YTick', ptsTicks);


            else
                hRef = line([xMin xMax],[0 0]);

                %Set Labels:
                hXLab = xlabel('Date');
                hYLab = ylabel([varCurr ' (' unitCurr ')']);

                %EDIT X-TICK LABELS:
                xMin = max(xMin, 1);
                %Create dates for x-axis:
                strDates = date2str(unique(days_2_date((xMin:xMax),dataAll{indCurrType{ll}(indCurr),1}{2}(1,:),'gregorian'),'rows'),'m/d/y');
                %Set labels
                if ischar(strDates)
                    nXTicks = 1;
                else
                    nXTicks = min(numel(strDates),10);
                end
                ptsDates = round(linspace(1,numel(strDates),nXTicks));
                ptsTicks =  round(linspace(xMin,xMax,nXTicks));
                set(gca, 'XTickLabel',strDates(ptsDates), 'XTick', ptsTicks);
                set(gca, 'XTickLabelRotation', 45);
                %Set axis and data fram properties:
                xlim([xMin-5,xMax+5]); 
            end
            hold off

            %Create Legend:
            hLgd = legend(hTs,strSeries,'location','bestoutside');


            %Set plot properties:
            set(gca, ...
                'FontName'   , strFont);
            set([hXLab, hYLab], ...
                'FontName'   , strFont);
            set(hLgd, ...
                'FontSize'   , ftSzLgd, ...
                'LineWidth', lnWdA);
            set([hXLab, hYLab]  , ...
                'FontSize'   , ftSz, ...
                'Color', 'black');
            set(hRef(:),'LineWidth', lnWd, ...
                'LineStyle','--', ...
                'color',custGray);
            set(gca, ...
                'Box'         , 'off', ...
                'TickDir'     , 'in'     , ...
                'TickLength'  , [.02 .02] , ...
                'LineWidth'   , lnWdA , ...
                'fontSize', ftSzAx);   
            if ~flagScatter
                set(gca, 'YMinorTick'  , 'on');
            else
                set(hText, ...
                    'FontSize'   , ftSzLgd, ...
                    'LineWidth', lnWdA);
            end

            if isdir(outputPath)
                dirWrite = outputPath;
                fileWrite = '';
            else
                [dirWrite, fileWrite, ~]= fileparts(outputPath);
            end
            if ~isempty(fileWrite)
                pathWrite = fullfile(dirWrite, [fileWrite '_' num2str(ii)]);
            else
                pathWrite = fullfile(dirWrite, ['mod_v_obs_' num2str(ii)]);
            end

            if flagScatter
                pathWrite = [char(pathWrite), '_scatter'];
            else
                pathWrite = [char(pathWrite), '_ts'];
            end

            savefig(hFig, [pathWrite '.fig']);
            print(hFig, [pathWrite '.eps'],'-depsc2');
            print(hFig, [pathWrite '.png'],'-dpng','-r600');
        end
    end
end

if nargout >= 1
    varargout{1} = scoreOut; 
	if nargout == 2
        for kk = 1 : nMetric
            %Rename specific obsure output names:
            for ii = 1 : numel(nmScore{kk}(:))
                if regexpbl(nmScore{kk}{ii}, 'casi')
                   nmScore{kk}{ii} = 'cryosphere covered area'; 
                end
            end
        end
    	varargout{2} = nmScore;
	end
end
