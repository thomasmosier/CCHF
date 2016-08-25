function varargout = mod_v_obs_v2(sObs, sMod, methRank, varargin)
%Match all observation data to model output for observations taken at one
%point in time or over a duration (must have fields '.dateStart' and
%'.dateEnd') and calculate performance statistic.  Currently
%all observations of a same type are lumped together before being compared
%(e.g. if many stake measurements present, these will be lumped together
%before calculating a performance statistic).



flagPlot = 0;
flagScatter = 0;
typeCombine = '';
flagOneStat = 0;
statCombine = 'linear';
outputPath = '';
rmMod = '';
if ~isempty(varargin)
    for ii = 1 : numel(varargin) 
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
   elseif regexpbl(obsTypes{ii}, 'path_')
        obsTypes(ii) = [];
    elseif regexpbl(obsTypes{ii}, 'type_')
        obsTypes(ii) = [];
   end
end

obsTypes = unique(obsTypes);

%Create output array for observations-model comparisons:
dataAll = cell(numel(namesObs(:)),5);
%Format is: 
    %{ii,1} = {Measurement Type; reference date dateRef; [lon, lat]; eval metric};
    %{ii,2} = dates/days (Model)
    %{ii,3} = data (Model)
    %{ii,4} = days (Obs)
    %{ii,5} = data (Obs)


%Find matches between observation and model output:
for ii = 1 : numel(namesObs(:)) %Loop over all data in the Observation data
    indCurr = strcmpi(namesObs{ii},namesMod);
    
    %Initialize metadata cell array:
    dataAll{ii,1} = cell(4,1);
    
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
            if regexpbl(nmDataCurrObs{mm}, 'path_')
                nmDataCurrObs(mm) = [];
            elseif regexpbl(nmDataCurrObs{mm}, 'type_')
                nmDataCurrObs(mm) = [];
            end
        end

    %Loop over fields in observation data (expected that one
    %modelled output exists for each observation field
    for kk = 1 : numel(nmDataCurrObs(:))
        
        %Find indices for observation fields:
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
        elseif regexpbl(nmDataCurrObs{kk},{'snow','radar'},'and') || regexpbl(nmDataCurrObs{kk},'swe')
            indField = find(strcmpi(fieldsCurrMod, 'swe') == 1);
            
            evalType = 'sum';
        elseif regexpbl(nmDataCurrObs{kk}, {'casi', 'icx', 'snx', 'sca'})
            %This is tricky because the data to compare with depends on the
            %metric. If metric is ''Parajka' (following Parajka & Blöschl,
            %2008), then compare MODIS SCA with modelled SWE. Otherwise,
            %compare observed SCA with modelled SCA.
            if regexpbl(methRank, 'Parajka')
                indField = find(strcmpi(fieldsCurrMod, 'snw') == 1);
            else
                indField = find(strcmpi(fieldsCurrMod, nmDataCurrObs{kk}) == 1);
            end
            
            if isfield(sObs.(namesObs{ii}),['type_' nmDataCurrObs{kk}])
                evalType = sObs.(namesObs{ii}).(['type_' nmDataCurrObs{kk}]);
            else
                evalType = 'max';
                warning('mod_v_obs:unknownEvalType',['The evaluation '...
                    'type is assumed to be ' evalType ' for aggregating ' ...
                    nmDataCurrObs{kk} ' observations.']);
            end
        else
            indField = find(strcmpi(fieldsCurrMod, nmDataCurrObs{kk}) == 1);
        end
        
        %Select appropriate evaluation metric ('Parajka' can only be used for snow covered area)
        if regexpbl(nmDataCurrObs{kk}, {'casi', 'icx', 'snx', 'sca'})
            if regexpbl(methRank, 'Parajka')
                dataAll{ii,1}{4} = 'Parajka';
            else
                dataAll{ii,1}{4} = methRank;
            end
        else
            if regexpbl(methRank, 'Parajka') %Remove 'Parajka' metric reference
                indMeth = regexpi(methRank,'Parajka');
                if indMeth(1) == 1
                    indMeth = [8, numel(methRank)];
                else
                    indMeth = [1, indMeth - 1];
                end
                dataAll{ii,1}{4} = methRank(indMeth(1):indMeth(2));
            else
                dataAll{ii,1}{4} = methRank;
            end
        end
        %Remove underscores from start or end of method rank string:
        if regexpbl(dataAll{ii,1}{4}(1),'_')
            dataAll{ii,1}{4} = dataAll{ii,1}{4}(2:end);
        end
        if regexpbl(dataAll{ii,1}{4}(end),'_')
            dataAll{ii,1}{4} = dataAll{ii,1}{4}(1:end-1);
        end
        

        %Retrieve dates and temporal resolution of modelled data:
        dataAll{ii,2} = sMod.(namesMod{indCurr}).date;
        tResMod = numel(dataAll{ii,2}(1,:));
        
        
        %Retrieve data from model:
        if ~isempty(indField) %Modelled field do not exist for current observation field
            if numel(indField) == 2
                dataAll{ii,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(1)}) ...
                        + sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(2)});
            elseif numel(indField) == 1
                dataAll{ii,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField});
            end
            
            %SPECIAL TREATMENT FOR STAKE DATA:
%             if regexpbl(nmDataCurrObs{kk},'stake') 
%                 if all(sMod.(namesMod{indCurr}).('iceBl')) == 1
%                     dataAll{ii,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(1)}) ...
%                         + sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(2)});
%                 else
%                     dataAll{ii,3} = nan(numel(sMod.(namesMod{indCurr}).(fieldsCurrMod{indField(1)})),1);
%                 end
%             else
%                 dataAll{ii,3} = sMod.(namesMod{indCurr}).(fieldsCurrMod{indField});
%             end
        else
            warning('mod_v_obs:missingModOutput',['Field ' nmDataCurrObs{1} ...
                ' at observation ' namesObs{ii} ' is being skipped '...
                'because no modelled indices exist at this point.']);
            continue
        end 

        %Retrieve data from observation:
        if  numel(nmDataCurrObs) == 1
            dataAll{ii,5} = sObs.(namesObs{ii}).(nmDataCurrObs{kk});
            
            %If there is a 'flag' field, value of 1 means don't use data in
            %analysis.
            if isfield(sObs.(namesObs{ii}), 'flag')
                %Different flag treatment for 3D vs. 2D arrays
                if numel(size(dataAll{ii,5})) == 3
                    if numel(sObs.(namesObs{ii}).flag) == numel(dataAll{ii,5}(:,1,1))
                        for mm = 1 : numel(dataAll{ii,5}(:,1,1))
                            if sObs.(namesObs{ii}).flag(mm) == 1
%                                 if any2d(~isnan(squeeze(dataAll{ii,5}(mm,:,:)))) 
%                                     disp(num2str(mm))
%                                 end
                                dataAll{ii,5}(mm,:,:) = nan;
                            end
                        end
                    else
                        error('mod_v_obs:error3DFlagSz','The data quality flag has an unexpected size.');
                    end
                elseif numel(size(dataAll{ii,5})) == 2
                    if numel(sObs.(namesObs{ii}).flag) == numel(dataAll{ii,5})
                        dataAll{ii,5}(sObs.(namesObs{ii}).flag == 1) = nan;
                        
                    else
                        error('mod_v_obs:error2DFlagSz','The data quality flag has an unexpected size.');
                    end
                end
            end
        elseif numel(nmDataCurrObs) > 1
            error('mod_v_gage:extraObsFields','Extra data fields exist for current observation variable.');
        else
            error('mod_v_gage:noObsFields','No data fields exist for current observation variable.');
        end

        
        if isfield(sObs.(namesObs{ii}), 'date') %Observations occur at single point in time
            dataAll{ii,4} = sObs.(namesObs{ii}).date;
            tResObs = numel(dataAll{ii,4}(1,:));

            %If temporal resolution doesn't match, convert (e.g. 
            %if observation data monthly and modelled data daily):
            if tResObs ~= tResMod 
                if tResObs == 3 && tResMod == 2
                    datesAvg = unique(dataAll{ii,4}(:,1:2),'rows');
                    dataAvg = nan(numel(datesAvg(:,1)),1);
                    for ll = 1 : numel(datesAvg(:,1))
                        indAvg = ismember(dataAll{ii,4}(:,1:2), datesAvg(ll,1:2),'rows');
                        dataAvg(ll) = mean(dataAll{ii,5}(indAvg));
                    end

                    dataAll{ii,5}  = dataAvg; %Substitute averages for observations
                    dataAll{ii,4} = datesAvg;
                elseif tResObs == 2 && tResMod == 3
                    datesAvg = unique(dataAll{ii,2}(:,1:2),'rows');
                    dataAvg = nan(numel(datesAvg(:,1)),1);
                    for ll = 1 : numel(datesAvg(:,1))
                        indAvg = ismember(dataAll{ii,2}(:,1:2), datesAvg(ll,1:2),'rows');
                        dataAvg(ll) = mean(dataAll{ii,3}(indAvg));
                    end

                    dataAll{ii,3}  =  dataAvg; %Substitute averages for modelled
                    dataAll{ii,2} = datesAvg;
                else
                    error('mod_v_obs:tempResDiff',['Each '...
                        'time-series element for the observed '...
                        'data has ' ...
                        num2str(numel(dataAll{ii,4}(1,:))) ...
                        ' componenents while the modelled data '...
                        'has ' num2str(numel(dataAll{ii,2}(1,:)))...
                        '.  This has not been programmed for.']);
                end
            end

            %Find dates common to both modelled and observed data
            [~, dateUseObs, dateUseMod] = intersect(dataAll{ii,4}, dataAll{ii,2},'rows');

            dataAll{ii,2} = dataAll{ii,2}( dateUseMod,:);
            dataAll{ii,3} = dataAll{ii,3}( dateUseMod,:);
            dataAll{ii,4} = dataAll{ii,4}( dateUseObs,:);
            dataAll{ii,5} = dataAll{ii,5}( dateUseObs,:);

            if sum(~isnan(dateUseObs)) < 0.7*sum(~isnan(dataAll{ii,5}))
                warning('mod_v_obs:unusedObs',...
                    [num2str(sum(~isnan(dateUseObs))) ' of the ' ...
                    num2str(sum(~isnan(dataAll{ii,5}))) ...
                    ' observation time-points are being '...
                    'used.  Consider adapting the model '...
                    'run time accordingly.']);
            end
        elseif isfield(sObs.(namesObs{ii}), 'dateStart') && isfield(sObs.(namesObs{ii}), 'dateEnd')
            dataAll{ii,4} = cell(2,1);
            dataAll{ii,4}{1} = sObs.(namesObs{ii}).dateStart;
            dataAll{ii,4}{2} = sObs.(namesObs{ii}).dateEnd;
            tResObs = numel(dataAll{ii,4}{1}(1,:));


            %Match observation and modelled dates differently
            %depending on resolutions:             
            if tResObs ~= tResMod 
                if tResObs == 3 && tResMod == 2
                    %Interpolate modelled data to daily resolution
                    %and extract for same dates present in
                    %observations

                    datesTempModInterp = date_vec_fill([dataAll{ii,2}(1,:),1],[dataAll{ii,2}(end,:),eomday(dataAll{ii,2}(end,1), dataAll{ii,2}(end,2))],'gregorian');

                    warning('off', 'all'); %Warning caused if all modelled data are NaN
                    dataAll{ii,3} = interp_month2day(dataAll{ii,2}, dataAll{ii,3}, datesTempModInterp, 'pchip', 1);
                    warning('on', 'all');
                    dataAll{ii,2} = datesTempModInterp;
                else
                    error('mod_v_obs:tempResDiff',['Each '...
                        'time-series element for the observed '...
                        'data has ' ...
                        num2str(numel(dataAll{ii,4}(1,:))) ...
                        ' componenents while the modelled data '...
                        'has ' num2str(numel(dataAll{ii,2}(1,:)))...
                        '.  This has not been programmed '...
                        'for.  The observation data has '...
                        'start and end dates']);
                end
            end
            
            
            %Extract modelled data for time-series present in observations:
            if numel(size(dataAll{ii,5})) == 2 %Pt data
                dataTempMod2Obs = nan(numel(dataAll{ii,5}),1);
                for ll = 1 : numel(dataAll{ii,5})
                    [~, indUpMod(1)] = ismember(dataAll{ii,4}{1}(ll,:), dataAll{ii,2}, 'rows');
                    [~, indUpMod(2)] = ismember(dataAll{ii,4}{2}(ll,:), dataAll{ii,2}, 'rows');

                    if all(indUpMod ~= 0) && indUpMod(2) > indUpMod(1) %If start and end dates of observation present in model, copy data
                        if regexpbl(evalType, 'sum')
                            dataTempMod2Obs(ll) = nansum(dataAll{ii,3}(indUpMod(1):indUpMod(2)));
                        elseif regexpbl(evalType, 'mean')
                            dataTempMod2Obs(ll) = nanmean(dataAll{ii,3}(indUpMod(1):indUpMod(2)));
                        elseif regexpbl(evalType, 'max')
                            dataTempMod2Obs(ll) = nanmax(dataAll{ii,3}(indUpMod(1):indUpMod(2)));
                        else
                            error('mod_v_obs:unknownEvalType', [evalType ...
                                ' is an unknown evaluation type and a '...
                                'case has not been programmed for it.']);
                        end
                    else
                        dataAll{ii,5}(ll) = nan;
                    end
                end
            elseif numel(size(dataAll{ii,5})) == 3 %gridded data
                dataTempMod2Obs = nan(size(dataAll{ii,5}),'single');

                for ll = 1 : numel(dataAll{ii,5}(:,1,1))
                    [~, indUpMod(1)] = ismember(dataAll{ii,4}{1}(ll,:), dataAll{ii,2}, 'rows');
                    [~, indUpMod(2)] = ismember(dataAll{ii,4}{2}(ll,:), dataAll{ii,2}, 'rows');

                    if all(indUpMod ~= 0) && indUpMod(2) > indUpMod(1) %If start and end dates of observation present in model, copy data
                        if regexpbl(evalType, 'sum')
                            dataTempMod2Obs(ll,:,:) = squeeze(sum(dataAll{ii,3}(indUpMod(1):indUpMod(2),:,:),[], 1));
                        elseif regexpbl(evalType, 'mean')
                            dataTempMod2Obs(ll,:,:) = squeeze(mean(dataAll{ii,3}(indUpMod(1):indUpMod(2),:,:),[], 1));
                        elseif regexpbl(evalType, 'max')
                            dataTempMod2Obs(ll,:,:) = squeeze(max(dataAll{ii,3}(indUpMod(1):indUpMod(2),:,:),[], 1));
                        else
                            error('mod_v_obs:unknownEvalType', [evalType ...
                                ' is an unknown evaluation type and a '...
                                'case has not been programmed for it.']);
                        end
                    else
                        dataAll{ii,5}(ll,:,:) = nan;
                    end
                end
            else
                error('mod_v_obs:errDim', ['The current variable has ' ...
                    num2str(numel(size(numel(dataAll{ii,5})))) ...
                    ' dimensions, which has not been programmed for.']);
            end

            %Copy aggregated model data to data array:
            dataAll{ii,3} = dataTempMod2Obs;
            
            %Remove all empty elements for point observation data:
            if numel(size(dataAll{ii,5})) == 2
                indRemCurr = find(isnan(dataAll{ii,5}));
                if ~isempty(indRemCurr)
                    dataAll{ii,5}(indRemCurr) = []; %Observation data field
                    dataAll{ii,3}(indRemCurr) = []; %Modeled data field
                    dataAll{ii,4}{1}(indRemCurr,:) = []; %Observation start date
                    dataAll{ii,4}{2}(indRemCurr,:) = []; %Observation end date
                end
            end
            
            %Set model output dates to be same as observation (because
            %model output has been aggregated to observation time elements)
            dataAll{ii,2} = dataAll{ii,4}; %Model start date
        else
            error('mod_v_obs:noDate',['No time field appears to be present for ' sObs.(namesObs{ii}) ' observation data.'])
        end


        %Save time signature as reference date and "days_since" instead of
        %date:
        if iscell(dataAll{ii,2}) && iscell(dataAll{ii,4})
            if ~isempty(dataAll{ii,2}{1}) && numel(dataAll{ii,2}{1}(1,:)) > 1 && ~isempty(dataAll{ii,4}{1}) && numel(dataAll{ii,4}{1}(1,:)) > 1
                dateRef = dataAll{ii,2}{1}(1,:);
                for ll = 1 : numel(dataAll{ii,2})
                    dataAll{ii,2}{ll} = days_since(dateRef, dataAll{ii,2}{ll}, 'gregorian');
                    dataAll{ii,4}{ll} = days_since(dateRef, dataAll{ii,4}{ll}, 'gregorian');
                end
            end
        elseif iscell(dataAll{ii,2})
            if ~isempty(dataAll{ii,2}{1}) && numel(dataAll{ii,2}{1}(1,:)) > 1
                dateRef = dataAll{ii,2}{1}(1,:);
                for ll = 1 : numel(dataAll{ii,2})
                    dataAll{ii,2}{ll} = days_since(dateRef, dataAll{ii,2}{ll}, 'gregorian');
                end
            end
            dataAll{ii,4} = days_since(dateRef, dataAll{ii,4}, 'gregorian');
        elseif iscell(dataAll{ii,4})
            if ~isempty(dataAll{ii,4}{1}) && numel(dataAll{ii,4}{1}(1,:)) > 1
                dateRef = dataAll{ii,4}{1}(1,:);
                for ll = 1 : numel(dataAll{ii,4})
                    dataAll{ii,4}{ll} = days_since(dateRef, dataAll{ii,4}{ll}, 'gregorian');
                end
            end
            dataAll{ii,2} = days_since(dateRef, dataAll{ii,2}, 'gregorian');
        else
            if ~isempty(dataAll{ii,2}) && numel(dataAll{ii,2}(1,:)) > 1
                dateRef = dataAll{ii,2}(1,:);
                dataAll{ii,2} = days_since(dateRef, dataAll{ii,2}, 'gregorian');
                dataAll{ii,4} = days_since(dateRef, dataAll{ii,4}, 'gregorian');
            end
        end

        
        %Write metadata for current observation-model:
            % {Measurement Type; reference date dateRef; [lon, lat]};
        %This could be defined within the ii loop, but by having it in the
        %kk loop it will be empty if no datafield found.
        dataAll{ii,1}(1:3) = {nmDataCurrObs{kk}; dateRef; [sObs.(namesObs{ii}).lon, sObs.(namesObs{ii}).lat]};
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
    end
end



%%CALCULATE FITNESS SCORES
if regexpbl(typeCombine,'type')
    scoreTemp = nan(numel(obsTypes), 1);
    nmScore = cell(numel(obsTypes), 1);

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
                    fitness( dataAll{indCurrType(indCurr),5}, dataAll{indCurrType(indCurr),3}, dataAll{indCurrType(indCurr),1}{4}, rmMod);
            end
        end
        
        %Combine scored for current type:
        if ~isempty(fitTemp)
            switch statCombine
                case 'linear'
                    scoreTemp(ii) = nanmean(fitTemp);
                case 'square'
                    scoreTemp(ii) = sqrt(nansum(fitTemp.^2));
            end

            nmScore{ii} = char(obsTypes{ii});
        end
    end
elseif regexpbl(typeCombine,'lump')
    scoreTemp = nan(numel(obsTypes), 1);
    nmScore = cell(numel(obsTypes), 1);
    for ii = 1 : numel(obsTypes)
        indCurrType = nan(numel(dataAllCmpr(:,1)),1);
        for jj = 1 : numel(dataAllCmpr(:,1))
            indCurrType(jj) = strcmpi(obsTypes{ii}, dataAllCmpr{jj,1}{1});
        end
        indCurrType = find(indCurrType == 1);
        
        dataFitTempMod = [];
        dataFitTempObs = [];
        for jj = 1 : numel(indCurrType)
            dataFitTempMod = [dataFitTempMod; dataAllCmpr{indCurrType(jj),3}];
            dataFitTempObs = [dataFitTempObs; dataAllCmpr{indCurrType(jj),5}];
        end
        
        if ~isempty(dataFitTempObs) && ~isempty(dataFitTempMod) 
            scoreTemp(ii) = fitness(dataFitTempObs, dataFitTempMod, dataAll{indCurrType(jj),1}{4}, rmMod);
            nmScore{ii} = char(obsTypes{ii});
        else
            scoreTemp(ii) = nan;
        end
    end
else
    scoreTemp = nan(numel(dataAll(:,1)),1);
    nmScore = cell(numel(dataAll(:,1)), 1);
    for ii = 1 : numel(dataAll(:,1))
        if ~isempty(dataAll{ii,5}) && ~isempty(dataAll{ii,3}) 
            scoreTemp(ii) = fitness(dataAll{ii,5}, dataAll{ii,3}, dataAll{ii,1}{4}, rmMod);
            nmScore{ii} = [char(dataAll{ii,1}{1}), ' (' num2str(round2(dataAll{ii,1}{3}(1),3)) ', ' num2str(round2(dataAll{ii,1}{3}(2),3)) ')'];
        end
    end
end

%Remove empty fields / scores
for ii = numel(nmScore) : -1 : 1
   if isempty(nmScore{ii})
       nmScore(ii) = [];
       scoreTemp = scoreTemp((1:numel(scoreTemp)) ~= ii); 
   end
end
    
if flagOneStat == 1
   %%COMBINE FITNESS SCORES:
    switch statCombine
        case 'linear'
            scoreOut = nanmean(scoreTemp);
        case 'square'
            scoreOut = sqrt(nansum(scoreTemp.^2));
    end 
    
    nmTemp = blanks(0);
    for ii = numel(nmScore)
        nmTemp = [nmTemp '-' nmScore{ii}];
    end
    
    nmScore = {'Combined' nmTemp};
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
            %Skip is metric is 
            if regexpbl(dataAll{indCurrType{ll}(1),1}{4}, 'Parajka')
               continue 
            end
            
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
                    strSeries{indCurr} = ['Pt (' num2str(round2(dataAll{indCurrType{ll}(indCurr),1}{3}(1),3)) ', ' num2str(round2(dataAll{indCurrType{ll}(indCurr),1}{3}(2),3)) ')'];
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
                                plotObsData  = [plotObsData, mean2d(squeeze(dataAll{indCurrType{ll}(indCurr),5}(kk)))*ones(1, numel(tempModDates))];
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
                    strCurrPt = ['Pt (' num2str(round2(dataAll{indCurrType{ll}(indCurr),1}{3}(1),3)) ', ' num2str(round2(dataAll{indCurrType{ll}(indCurr),1}{3}(2),3)) ')'];
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
                'TickDir'     , 'out'     , ...
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
        %Rename specific obsure output names:
        for ii = 1 : numel(nmScore(:))
            if regexpbl(nmScore{ii}, 'casi')
               nmScore{ii} = 'cryosphere covered area'; 
            end
        end
    	varargout{2} = nmScore;
	end
end
