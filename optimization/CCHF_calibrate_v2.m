function [sModOut, outParam] = CCHF_calibrate_v2(sObs, sOpt, sPath, sHydro, sMeta)


%Load option values:
sMeta.mode = 'calibrate';
nPop = find_att(sOpt.att,'population');
nGen = find_att(sOpt.att,'generations');

prmBnds = cell2mat(sMeta.coefAtt(:,2:3));
fitTest = sOpt.fitTest; %Needed for parfor.

%If multiple stages of calibration, (1) determine which parameters fit each
%category and (2) which observation data to use
%Initialize:
stagePrmNm = cell(sOpt.nStage,1);
stageBestPrm = cell(sOpt.nStage,1); %Values of calibrated parameters from each stage
stagePrmInd = cell(sOpt.nStage,1);
stageDInd = cell(sOpt.nStage,1);
stageNInd = cell(sOpt.nStage,1);
stagePInd = cell(sOpt.nStage,1);
stageObsStruc = cell(sOpt.nStage,1);

disp(['Calibration is initiating using the ' char(39) sOpt.type ...
    char(39) ' method.']);


%Populate a series of cell arrays that contain information such as which 
%parameters to calibrate during each stage:
if sOpt.nStage == 1 %All calibration being done in single stage
    stagePrmNm{1} = sMeta.coefAtt(:,1); %Contains names of parameters calibrated in each stage.
    stageObsStruc{1} = sObs; %Contains observations used in each calibration stage.
    stagePrmInd{1} = (1:numel(stagePrmNm{1}))'; %Indices (from sMeta.coef) of parameters calibrated in each stage.
    stageDInd{1} = []; %Indices of all parameters to be calibrated in future stages
    stageNInd{1} = []; %Indices of all parameters not used during current calibration stage
    stagePInd{1} = []; %Indices of all parameters calibrated in previous stage (cumulative)
else %Calibration conducted in multiple stages
    for ss = 1 : sOpt.nStage
        %Find parameters for each stage:
        stageTemp = cell(0,1);
        for jj = 1 : numel(sOpt.stageVar(ss,:))
            indVarUse = cellfun(@(x) ~isempty(regexpi(x,sOpt.stageVar{ss,jj})), sMeta.coef(:,6));
            stageTemp = [stageTemp; sMeta.coef(indVarUse,1)];
        end
        stagePrmNm{ss} = stageTemp;
        
        %Find all parameter indices not in current stage:
        for jj = 1 : numel(sMeta.coef(:,1))
            if any(strcmpi(sMeta.coef(jj,1),stagePrmNm{ss}(:)))
                stagePrmInd{ss} = [stagePrmInd{ss}; jj];
            else
                stageNInd{ss} = [stageNInd{ss}; jj];
            end
        end
        
        %Find all indices calibrated in previous stages:
        if ss == 1
           stagePInd{ss} = []; 
        else
            for ll = 1 : ss - 1
                stagePInd{ss} = [stagePInd{ss}; stagePrmInd{ll}];
                %This line allows parameters to be included in both
                %calibration stages (maybe not good idea?)
                stagePInd{ss} = setdiff(stagePInd{ss},stagePrmInd{ll+1});
            end
        end
        
        stageDInd{ss} = setdiff(stageNInd{ss}, stagePInd{ss});
        
        %Find observation data for each stage:
        ptsObs = fields(sObs);
        blUse = zeros(numel(ptsObs),1);
        for jj = 1 : numel(ptsObs)
            fieldsCurr = fields(sObs.(ptsObs{jj}));
            %Remove metadata fields:
            fieldsCurr(strcmpi(fieldsCurr,'lat')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'latitude')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'lon')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'longitude')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'att')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'attributes')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'flag')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'date')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'dateStart')) = [];
            fieldsCurr(strcmpi(fieldsCurr,'dateEnd')) = [];
            
            %Four categories are 'input', 'cryo', 'land', 'routing'
            if any(strcmpi(fieldsCurr{1},{'flow','flowrate','streamflow'}))
                typeCurr = 'routing';
            elseif any(strcmpi(fieldsCurr{1},{'stake','swe','snowradar','sca','casi','snowcover'}))
                typeCurr = 'cryo';
            elseif any(strcmpi(fieldsCurr{1},{'tmp','tas','tmx','tmn','pre','pr'}))
                typeCurr = 'input';
            else
                warning('CCHF_calibrate:unkownObs',[char(39) ...
                    fieldsCurr{1} char(39) ' observation data has not '...
                    'been programmed for in the calibration routine. '...
                    'Please add.']);
            end

            if any(strcmpi(typeCurr, sOpt.stageVar(ss,:)))   
                blUse(jj) = 1;
            end  
        end
        
        indUse = find(blUse == 1);
        for jj = 1 : numel(indUse)
            stageObsStruc{ss}.(ptsObs{indUse(jj)}) = sObs.(ptsObs{indUse(jj)});
        end
        
        if isempty(stageObsStruc{ss})
            error('CCHF_calibrate:obsCurrPt',['No observation data '...
                'were found for stage ' num2str(ss) ' of calibration point.']);
        end
    end
end



%%ALLOW OPTION TO RESUME PREVIOUSLY INTERUPPTED CALIBRATION
%Requires loading previously tested parameter sets from files
blFresh = 1; indStgStrt = 1;
if regexpbl(sMeta.runType,{'calib','resume'},'and')
    if isfield(sPath,'resume')
        filesParam = dir(fullfile(sPath.resume, '*.csv'));
        filesParam = extractfield(filesParam, 'name');
        if ~isempty(filesParam)
            indStage = nan(numel(filesParam),1);
            for ii = 1 : numel(filesParam)
                indUnd = regexpi(filesParam{ii}, '_');
                indStage(ii) = str2double(filesParam{ii}(indUnd(end)+1:end-4));
            end
            
            %Determine which stage to start at (based on parameter files in
            %folder)
            indStgStrt = intersect((1:sOpt.nStage),indStage);
            indStgStrt = indStgStrt(end);
            
            %Find which generation within the stage to begin with (based on 
            %number of parameter sets in the file associated with the starting
            %stage):
            if ~isempty(indStgStrt)
                %Load parameter sets from current and previous stage:
                if indStgStrt ~= 1
                    indStgLd = [indStgStrt-1, indStgStrt];
                else
                    indStgLd = indStgStrt;
                end
                
                for ll = 1 : numel(indStgLd)
                    filePrm = fullfile(sPath.resume, filesParam{indStage == indStgLd(ll)});
                    prmArray = csvread(filePrm,1,0);
                    fidPrm = fopen(filePrm,'r');
                    varPrm = textscan(fidPrm,'%s',numel(prmArray(1,:)),'Delimiter',',');
                        varPrm = varPrm{1};
                    fclose(fidPrm);

                    %All entries in varPrm are variables, except the last,
                    %which is the performance (based on the metric used). 
                    %Check that these match the expected values. Then use 
                    %these to populate coefFit and coefFamily
                    if all(strcmpi(sMeta.coef(:,1),varPrm(1:end-1))) && strcmpi(varPrm{end}, sOpt.fitTest)
                        nGenPrev = numel(prmArray(:,end))/nPop;

                        coefFamily  = nan(nGen, nPop, numel(sMeta.coef(:,1)));
                        coefFitness = nan(nGen, nPop);
                        %Reshape loaded parameters into desired order:
                        indRe = 0;
                        for ii = 1 : nGenPrev
                            indRe = (ii-1)*nPop + 1;
                            coefFamily(ii,:,:) = prmArray(indRe:indRe + nPop-1,1:end-1);
                            coefFitness(ii,:) =  prmArray(indRe:indRe + nPop-1,end);
                        end 
    
                        indGenRem = find(isnan(coefFitness(:,1)));
                        coefFitness = coefFitness(1 : indGenRem(1)-1,:);
                        coefFamily = coefFamily(1 : indGenRem(1)-1,:,:);

                        %If stage is not 1st, find the best performing
                        %parameter set from the previous stage (this is used
                        %later)
                        if numel(indStgLd) > 1 && ll < numel(indStgLd)
                            [~, rMin, cMin] = min2d(coefFitness);

                            stageBestPrm{ll} = num2cell(squeeze(coefFamily(rMin, cMin,:)));
                        end
                        
                        indGenStrt = ones(sOpt.nStage,1);
                        indGenStrt(indStgStrt) = nGenPrev+1;
                        blFresh = 0;
                    else
                        warning('CCHF_calibrate:prmMismatch',['Not all '...
                            'parameters appear to match between those in '...
                            'unfinished calibration file and those needed '...
                            'for the current model. Therefore calibration '...
                            'will begin from the beginning.']);
                        blFresh = 1;
                    end
                end
            else
                warning('CCHF_calibrate:noStage',['No intersecting '...
                    'stages were found. Therefore calibration '...
                    'will begin from the beginning.']);
                blFresh = 1;
            end
            

        else
            warning('CCHF_calibrate:noPrmFiles',['No csv files '...
                'containing parameters from the unfinished '...
                'calibration were found. Therefore calibration '...
                'will begin from the beginning.']);
            blFresh = 1;
        end
    else
        warning('CCHF_calibrate:noResumePath',['No path to the '...
            'unfinished directory was found. Therefore calibration '...
            'will begin from the beginning.']);
        blFresh = 1;
    end
end

%If calibration starting from beginning or if resume did not find path to 
%start from
if blFresh
    indStgStrt = 1; %First stage of calibration.
    indGenStrt = ones(sOpt.nStage,1);
end






%Check that all indices have been located
checkInd = 0;
for ss = 1 : sOpt.nStage
   checkInd = checkInd + numel(stagePrmInd{ss});
end
if checkInd < numel(sMeta.coef(:,1))
   error('CCHF_calibrate:nIndWrong',[num2str(numel(sMeta.coef(:,1))) ...
       ' variables are present but only ' num2str(checkInd) ...
       ' have been located in the multi-stage calibration algorithm.']); 
end
    
    

%Initialize arrays if uniform sampling being used:
if regexpbl(sOpt.type, 'uniform') %If uniform parameter space sampling, generate all unique sets of coefficients:
    [coefTemp, nRun] = uniform_sample(prmBnds, nPop*nGen);

    %Reassign members and generation based on runs specified by uniform sampling: 
    f = factor(nRun);
    cntr = 0; nMemTemp = 0;
    while nMemTemp <= 20 
        cntr = cntr+1;
        if cntr > numel(f)
           cntr = numel(f); 
        end
        nMemTemp = prod(f(1:cntr));
    end
    nPop = prod(f(1:cntr));
    nGen = prod(f(cntr+1:end));

    %Initialize arrays:
    coefFamily  = nan(nGen, nPop, numel(sMeta.coef(:,1)));
    coefFitness = nan(nGen, nPop);

    %Reorganize coeffients based on new parameters:
    for ii = 1 : nGen
       coefFamily(ii,:,:) = coefTemp(nPop*(ii-1) + 1 : nPop*ii,:);
    end
end



%Create path and header for writing parameter sets:
%Set output path and header for writing coefficients:
pathGenParam = fullfile(sPath.output,['parameter_each_gen_' sMeta.region '.txt']);
pathBstCoef = fullfile(sPath.output,['coefficients_' sMeta.region '_' 'fittest' '.txt']);

hdrFull = {['Module_set: ' sMeta.strModule]; ...
    ['Method: ' sMeta.runType '_' sOpt.type]; ...
    ['Data_temporal_res: ' sMeta.dataTsRes]; ...
    ['Model_timestep: ' sMeta.dt]; ...
    blanks(0); ...
    'Parameter, value'};

%Define attributes to be used in cutting model runs short:
if isfield(sOpt,'dateEval')
    dateEval = sOpt.dateEval;
else
    dateEval = nan(1,3);
    warning('CCHF_calibrate:dateEval',['dateEval is not a field. '...
        'Therefore no short parameter set analysis will be run.']);
end


shortEvalAtt = {'fitTest', fitTest; 'dateEval', dateEval; 'fitReq', 1000};

%Open "Matlab pool" for parallel calibration:
if isempty(gcp('nocreate'))
    localCluster = parcluster('local'); %Check number of possible "workers"
    if isfield(sOpt,'maxWorker') && localCluster.NumWorkers > sOpt.maxWorker
        workers = sOpt.maxWorker;
    else
        workers = localCluster.NumWorkers;
    end
    parpool(workers); %Dedicate all available cores to parallel calibration
end




tic %Start clock that relays total time spent for calibration.
%Loop over calibration stages
for ss = indStgStrt : sOpt.nStage
    %Display messgage about start of current stage:
    if numel(stagePrmInd{ss}) == 0
        disp(['0 parameters have been identified to calibrate this '...
            'stage. Therefore, this stage will be skipped.']);
        continue
    else
        if blFresh
            disp(['Calibration stage ' num2str(ss) ' of ' ...
                num2str(sOpt.nStage) ' starting.' char(10) ...
                num2str(numel(stagePrmInd{ss})) ' of ' ...
                num2str(numel(sMeta.coef(:,1))) ' parameters are being '...
                'calibrated during this stage.']);
        else
            disp(['Calibration stage ' num2str(ss) ' of ' ...
                num2str(sOpt.nStage) ' resuming (starting with stage ' num2str(ss) ', iteration '...
                num2str(indGenStrt(ss)) ').' char(10) ...
                num2str(numel(stagePrmInd{ss})) ' of ' ...
                num2str(numel(sMeta.coef(:,1))) ' parameters are being '...
                'calibrated during this stage.']);
        end
    end
    %Assign new variable for observations to use in current stage
    %(otherwise MATLAB complains that variable is a broadcast variable.
    obsCurr = stageObsStruc{ss};
    
    %Initialize array to store all previously tried parameter sets from current
    %stage:
    coefRepo = nan(nGen*nPop, numel(sMeta.coef(:,1)));
    fitRepo  = nan(nGen*nPop, 1);
    nsf = 3; %Number of significant figures to use in comparing parameter sets (determine if they're the same)
    indSkip = [];
    indReset = 200;
    
    %%Initialize parameters (for all methods except uniform):
    if ~regexpbl(sOpt.type, 'uniform') && blFresh == 1
        %Initialize arrays:
        coefFamily  = nan(nGen, nPop, numel(sMeta.coef(:,1)));
        coefFitness = nan(nGen, nPop);

        %Generate n=(pop-1) random paremeter sets:
        coefFamily(1,1:end,stagePrmInd{ss}) = Monte_Carlo(prmBnds(stagePrmInd{ss},:), nPop);
        %Include the model's default parameter set as final member of population:
        coefFamily(1,end,stagePrmInd{ss}) = cell2mat(sMeta.coef(stagePrmInd{ss},4));
    end
    
    if ~regexpbl(sOpt.type, 'uniform') && (blFresh == 1 || (blFresh == 0 && ss == indStgStrt) )
        %Set parameters not being calibrated this stage or previous stages to default values:
        if ~isempty(stageDInd{ss})
            for ll = 1 : nGen
                for mm = 1 : nPop
                    coefFamily(ll, mm, stageDInd{ss}) = cell2mat(sMeta.coef(stageDInd{ss},4));
                end
            end
        end
        
        %Set parameters that were previously calibrated to the found values:
        if ~isempty(stagePInd{ss}) && ss ~= 1
            for ll = 1 : nGen
                for mm = 1 : nPop
                    coefFamily(ll, mm, stagePInd{ss}) = cell2mat(stageBestPrm{ss-1}(stagePInd{ss}));
                end
            end
        end
    end
    
    %Determine iterations after which early stagnation will be allowed:
    nParam = numel(stagePrmInd{ss});
    if nParam <= 4
        genMinStag = round(150/nPop + 3);
    elseif nParam <= 7
        genMinStag = round(300/nPop + 3);
    elseif nParam > 7 && nParam < 10
        genMinStag = round(500/nPop + 3);
    else
        genMinStag = round(800/nPop + 3);
    end
    
    
    
    %Loop over generations
    valMostFit = nan(nGen,1);
    for ii = indGenStrt(ss) : nGen
        %Calculate coefficients for current generation:
        if ii ~= 1 && ii ~= nGen && ~regexpbl(sOpt.type, 'uniform') %Don't optimize if uniform sampling used because all parameter sets pre-determined.
            coefFamily(ii,:,stagePrmInd{ss}) ...
                = CCHF_next_gen(coefFamily(:,:,stagePrmInd{ss}), coefFitness(:,:), ii-1, prmBnds(stagePrmInd{ss},:), sOpt);
        end
        
        %Array of current parameters:
        prmCurr = squeeze(coefFamily(ii,:,:));
        
        %Identify parameter sets that have already been run (so that they
        %can be skipped)
        if ii ~= 1
            indRepo = find(isnan(coefRepo(:,1)),1,'first');
            coefRepo(indRepo:indRepo + nPop - 1, : ) = roundsd(squeeze(coefFamily(ii-1,:,:)), nsf); 
            fitRepo(indRepo:indRepo + nPop - 1) = coefFitness(ii-1,:)';
            
            [blSkip, indRepo] = ismember(roundsd(prmCurr, nsf), coefRepo,'rows');
            indSkip = find(blSkip == 1);
            indRun = setdiff((1:nPop), indSkip);
            
            %Assign skipped variable scores to current generation:
            coefFitness(ii,indSkip) = fitRepo(indRepo(indSkip));
        else
            indRun = (1:nPop);
        end
        
        %Assign parameter sets to run this iteration:
        prmCurr = prmCurr(indRun,:);
        fitCurr = nan(numel(indRun),1);
        %Loop over each parameter set combination (in parallel if possible)
        shortCnt = 0;
            if ~isempty(gcp('nocreate'))
                pctRunOnAll warning('off','all') %Turn off warning during parloop
            end

        parfor jj = 1 : numel(indRun)
            if all(~isnan(dateEval))
                [sModTemp, fitTemp] = CCHF_engine_v3(sPath, sHydro, sMeta, prmCurr(jj,:)', obsCurr, shortEvalAtt);
                if isnan(fitTemp)
                    fitCurr(jj) = nanmean(mod_v_obs_v2(obsCurr, sModTemp, fitTest, 'combineType'));
                else
                    shortCnt = shortCnt + 1;
                    fitCurr(jj) = fitTemp;
                end
            else
                sModTemp = CCHF_engine_v3(sPath, sHydro, sMeta, prmCurr(jj,:)');
                fitCurr(jj) = nanmean(mod_v_obs_v2(obsCurr, sModTemp, fitTest, 'combineType'));
            end
        end
            if ~isempty(gcp('nocreate'))
                pctRunOnAll warning('on','all') %Turn warnings back on
            end
            
        %Assign fitness scores from current generation:
        coefFitness(ii,indRun) = fitCurr;
            
        if all(~isnan(dateEval)) % & && ii >= genStagTest - 1
           if sum2d(coefFitness < 50) >= nPop && shortEvalAtt{3,2} > 200
               shortEvalAtt{3,2} = 200;
               indReset = ii + 2;
           elseif sum2d(coefFitness < 1) >= 2*nPop && shortEvalAtt{3,2} > 100 %&& ii >= indReset
               shortEvalAtt{3,2} = 100;
               indReset = ii + 2;
           elseif sum2d(coefFitness < 1) >= 3*nPop && shortEvalAtt{3,2} > 50 %&& ii >= indReset
               shortEvalAtt{3,2} = 50;
               indReset = ii + 2;
           elseif sum2d(coefFitness < 1) >= 4*nPop && shortEvalAtt{3,2} > 10 %&& ii >= indReset
               shortEvalAtt{3,2} = 10;
               indReset = ii + 2;
           elseif sum2d(coefFitness < 1) >= 5*nPop && shortEvalAtt{3,2} > 5 %&& ii >= indReset
               shortEvalAtt{3,2} = 5;
               indReset = ii + 2;
           end
        end
        



        %Calculate best fitness from current generation:
        if ii == indGenStrt(ss) && ii ~= 1
            for zz = 1 : ii
                [valMostFit(zz), memFit] = min(coefFitness(zz,:));
                %If NSE used, convert to regular expression (1 - error); in 
                %'fitness' NSE is calculate only as error)
                if regexpbl(sOpt.fitTest, {'NSE','KGE'}) 
                    valMostFit(zz) = 1 - valMostFit(zz);
                end
            end
        else
            [valMostFit(ii), memFit] = min(coefFitness(ii,:));
            %If NSE used, convert to regular expression (1 - error); in 
            %'fitness' NSE is calculate only as error)
            if regexpbl(sOpt.fitTest, {'NSE','KGE'}) 
                valMostFit(ii) = 1 - valMostFit(ii);
            end
        end
        stdError = std(coefFitness(ii,:));


        
        %Write best parameter set from current generation to file:
        hdrShortCoef = {blanks(0); ...
            ['Stage: ' num2str(ss) ' of ' num2str(sOpt.nStage)];...
            ['Generation: ' num2str(ii)];...
            [sOpt.fitTest ': ' num2str(valMostFit(ii))]};
        if ss == 1 && ii == 1
            write_CCHF_coef(pathGenParam, [hdrFull(:); blanks(0); hdrShortCoef(:)], [sMeta.coef(:,1), num2cell(squeeze(coefFamily(ii,memFit,:)))]); 
        else
            write_CCHF_coef(pathGenParam, hdrShortCoef, [sMeta.coef(:,1), num2cell(squeeze(coefFamily(ii,memFit,:)))],'permissions','a'); 
        end
        
        %Record best set of parameters from all generations in current
        %stage:
        [~, indGenFit, indPopFit] = min2d(coefFitness(1:ii,:));
        stageBestPrm{ss} = num2cell(squeeze(coefFamily(indGenFit, indPopFit,:)));

        %Display running time for calibration at end of each generation:
        c = clock;
        disp([char(10) num2str(ii) ' of ' num2str(nGen) ' iterations have completed.'...
            char(10) 'The best fit ' sOpt.fitTest ' value for the current generation is ' ...
            num2str(round2(valMostFit(ii),2)) '.' char(10) ...
            'The standard deviation of fitnesses is ' num2str(round2(stdError,3)) '.' char(10) ...
            'Total elapsed time is ' num2str(round2(toc/60,1)) ' minutes. ' ...
            'The current time is ' num2str(round2(c(4),2)) ':' num2str(round2(c(5),2)) '.']);
        if shortCnt > 0
           disp(['GOOD NEWS! ' num2str(shortCnt) ' of ' num2str(nPop) ...
               ' model parameter runs were shortened. The error threshold is ' ...
               num2str(shortEvalAtt{3,2}) ' More efficiency!']); 
        end
        if numel(indSkip) ~= 0
            disp(['GOOD NEWS! ' num2str(numel(indSkip)) ' of ' num2str(nPop) ...
                ' model parameter runs were skipped because they have already been run. More efficiency!']);
        end
        

        %%WRITE all parameter sets to file:
        %Set output path and header for writing coefficients:
        if ii == 1
            pathCalib = fullfile(sPath.output,['all_parameter_sets-stage_' num2str(ss) '.csv']);
            fParam = fopen(pathCalib,'w+');
            
            numColumns = size(coefFamily,3)+1;
            headerFmt = repmat('%s,',1,numColumns-1);
            numFmt = repmat('%f,',1,numColumns-1);
            hdr = [sMeta.coef(:,1)', sOpt.fitTest];

            fprintf(fParam,[headerFmt,'%s\n'], hdr{:});
            fclose(fParam);
        elseif blFresh == 0 && ss == indStgStrt && ii == indGenStrt(ss)
            pathCalib = fullfile(sPath.output,['all_parameter_sets-stage_' num2str(ss) '.csv']);
            numColumns = size(coefFamily,3)+1;
            numFmt = repmat('%f,',1,numColumns-1);
        end
        fParam = fopen(pathCalib,'a+');
        %coefFamily  = nan(nGen, nPop, numel(sMeta.coef(:,1)));
        paramWrite = [squeeze(reshape(coefFamily(ii,:,:), 1, [], numel(sMeta.coef(:,1)))), reshape(coefFitness(ii,:),[],1)]; 
        fprintf(fParam,[numFmt,'%f\n'], paramWrite');
        fclose(fParam);

        %Record best parameter set so far to file (this way if calibration
        %quits unexpectedly, the parameters will still be there.
        currTime = clock;
        strNow = [num2str(currTime(1)) '-' num2str(currTime(2)) '-' ...
            num2str(currTime(3)) '-' num2str(currTime(4))];
        hdrBst = [hdrFull(1:4,:); ...
            [sOpt.fitTest ' fitness score: ' num2str(coefFitness(indGenFit, indPopFit))]; ...
            ['Date: ' strNow]; ...
            hdrFull(5:end,:)];
        write_CCHF_coef(pathBstCoef, hdrBst, [sMeta.coef(:,1), stageBestPrm{ss}],'permissions','w');
   
        
        %Allow early convergence if improvement has stagnated:
        if isfield(sOpt,'stagnate') && ii > genMinStag
            if regexpbl(sOpt.fitTest,{'NSE','KGE','Pearson'})
                [~, indGenFit] = max(valMostFit(1:ii));
            else
                [~, indGenFit] = min(valMostFit(1:ii));
            end
            
            stagGen = ii - sOpt.stagnate;
            if stagGen > 0 && ii - genMinStag > 0
                %Clause checks that best fit values are same for x generations
                if all(round2(valMostFit(stagGen), 2) == round2(valMostFit(stagGen : ii), 2))
%                 if all(round2(valMostFit(stagGen),2) >= round2(valMostFit(stagGen + 1 : ii),2)) && sum(round2(valMostFit(ii-genStagTest:ii),2) >= round2(valMostFit(stagGen),2)) > 1
                    disp([char(10) 'This stage of optimization is ending early.' char(10) 'The best fit ' sOpt.fitTest ' value has stagnated to ' ...
                        num2str(round2(valMostFit(indGenFit),2)) '.']);
                    blFresh = 1; %Reset blFresh if set to 0 using resume section of code
                    break
                end
            elseif std(reshape(coefFitness(1:ii,:),1,[])) < 0.1
                    disp([char(10) 'This stage of optimizating ending early because results are insenitive to parameters.' ...
                        char(10) 'The best fit ' sOpt.fitTest ' value has stagnated to ' ...
                        num2str(round2(valMostFit(indGenFit),2)) '.']);
                    blFresh = 1; %Reset blFresh if set to 0 using resume section of code
                    break
            end
        end
    end %End of generation loop
    
    
    %Record best value to display
    if regexpbl(sOpt.fitTest,{'NSE','KGE','Pearson'})
        [~, indGenFit] = max(valMostFit(1:ii));
    else
        [~, indGenFit] = min(valMostFit(1:ii));
    end
    
    %Stage completion message:
    disp(['CCHF model calibration stage ' num2str(ss) ' of ' num2str(sOpt.nStage) ' has completed.' ...
        char(10) 'The fittest set of parameter coefficients tested yielded an ' sOpt.fitTest ...
        ' value of ' num2str(round2(valMostFit(indGenFit),2)) ...
        ', ' char(10) 'which occured during iteration ' num2str(indGenFit) ...
        ' of ' num2str(nGen) ' total iterations.' char(10)]);
    
    
    %RECORD N BEST PERFORMING PARAMETER SETS:
    %This is potentially useful for assessing equifinality
    bestPrmSets = reshape(coefFamily,[], numel(sMeta.coef(:,1))); %indices are: (generation, pop member, parameter)
    [bestPrfm, indRank] = sort(reshape(coefFitness,[],1));
    bestPrmSets = bestPrmSets(indRank,:);
    
    %If number of best parameter options not specified, record 0.5% of all
    %test or 50, whichever is greater.
    if ~isfield(sOpt, 'nbest')
       sOpt.nbest = max(0.005*nGen*nPop, 50); 
    end
    bestPrmSets = bestPrmSets(1:sOpt.nbest,:);
    bestPrfm = bestPrfm(1:sOpt.nbest);
    
    pathBest = fullfile(sPath.output,[num2str(sOpt.nbest) '_parameter_sets-stage_' num2str(ss) '.csv']);
    fParamBest = fopen(pathBest,'w+');

    numColumnsBest = size(coefFamily,3)+1;
    headerFmtBest = repmat('%s,',1,numColumnsBest-1);
    numFmtBest = repmat('%f,',1,numColumnsBest-1);
    hdrBest = [sMeta.coef(:,1)', sOpt.fitTest];

    fprintf(fParamBest,[headerFmtBest,'%s\n'], hdrBest{:});
    paramWriteBest = [bestPrmSets, bestPrfm]; 
    fprintf(fParamBest,[numFmtBest,'%f\n'], paramWriteBest');
    fclose(fParamBest);
    
    %Calculate and record correlation between each parameters in N best sets:
    threshSig = 0.01; %correlation must be significant to less than this value
    prmCorr = nan(numel(bestPrmSets(1,:)));
    prmSig = nan(numel(bestPrmSets(1,:)));
    for mm = 1 : numel(bestPrmSets(1,:))
        for nn = 1 : numel(bestPrmSets(1,:))
            [r,p] = corrcoef(bestPrmSets(:,mm), bestPrmSets(:,nn));
            prmCorr(mm,nn) = r(2);
            prmSig(mm,nn) = p(2);
        end
    end
    prmCorr(prmSig > threshSig) = nan;

    %Write to file:
    [foldCf, ~, ~] = fileparts(sPath.coef);
    pathCorr = fullfile(foldCf, ['Pearson_corr-' num2str(sOpt.nbest) '_parameter_sets-stage_' num2str(ss) '.csv']);
    fCorr = fopen(pathCorr,'w+');

    nCol = numel(bestPrmSets(1,:)) + 1;
    headerFmt = repmat('%s,',1,nCol-1);
    numFmt = repmat('%f,',1,nCol-2);

    hdrCorr = [blanks(1), sMeta.coefAtt(:,1)'];

    fprintf(fCorr, [headerFmt,'%s\n'], hdrCorr{:});
    for mm = 1 : numel(bestPrmSets(1,:))
        fprintf(fCorr, '%s,', sMeta.coefAtt{mm,1});
        fprintf(fCorr, [numFmt,'%f\n'], prmCorr(mm,:));
    end
    fclose(fCorr);
    
    
    blFresh = 1; %Reset blFresh if set to 0 using resume section of code
end


%Close dedicated workers used for parallel processing:
poolobj = gcp('nocreate');
delete(poolobj);


%Rerun model for best parameter set:
sModOut = CCHF_engine_v3(sPath, sHydro, sMeta, stageBestPrm{end});


%Write and Plot parameters:
plot_parameters( stageBestPrm{end}, prmBnds, sMeta.coef, pathBstCoef);
hdrBst = [hdrFull(1:4,:); ...
    [sOpt.fitTest ' fitness score: ' num2str(valMostFit(indGenFit))]; ...
    ['Date: ' strNow]; ...
    hdrFull(5:end,:)];
write_CCHF_coef(pathBstCoef, hdrBst, [sMeta.coef(:,1), stageBestPrm{end}],'permissions','w');  

outParam = [sMeta.coef(:,1), stageBestPrm{end}, sMeta.coef(:,5:end)];