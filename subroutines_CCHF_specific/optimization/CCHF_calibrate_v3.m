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

function [sModOut, outParam] = CCHF_calibrate_v3(sObs, sOpt, sPath, sHydro, sMeta)


%Number of sites used in calibration
nSites = numel(sPath(:));

%Load option values:
sMeta.mode = 'calibrate';
nPop = find_att(sOpt.att,'population');
nGen = find_att(sOpt.att,'generations');
nStage = sOpt.nStage;
bestPrmSets = cell(nStage, 1);

if isfield(sOpt, 'methodSiteCombine')
    siteCombine = sOpt.methodSiteCombine;
else
    siteCombine = 'linear';
    if nSites ~= 1
        warning('CCHF_calibrate:unknownSiteCombine',['There is no '...
            'methodSiteCombine field, so the performance between model '...
            'sites is being set to mean']);
    end
end

prmBnds = cell2mat(sMeta.coefAtt(:,2:3));
fitTest = sOpt.fitTest; %Needed for parfor.


%If multiple stages of calibration, (1) determine which parameters fit each
%category and (2) which observation data to use
%Initialize:
stagePrmNm = cell(nStage,1);
stageBestPrm = cell(nStage,1); %Values of calibrated parameters from each stage
stagePrmInd = cell(nStage,1);
stageDInd = cell(nStage,1); %Indices of all parameters to be calibrated in future stages
stageNInd = cell(nStage,1); %Indices of all parameters not used during current calibration stage
stagePInd = cell(nStage,1); %Indices of all parameters calibrated in previous stage (cumulative)
stageObsStruc = cell(nStage,1);

disp(['Calibration is initiating using the ' char(39) sOpt.type ...
    char(39) ' method.']);


%Populate a series of cell arrays that contain information such as which 
%parameters to calibrate during each stage:
if nStage == 1 %All calibration being done in single stage
    stagePrmNm{1} = sMeta.coefAtt(:,1); %Contains names of parameters calibrated in each stage.
    stageObsStruc{1} = sObs; %Contains observations used in each calibration stage.
    stagePrmInd{1} = (1:numel(stagePrmNm{1}))'; %Indices (from sMeta.coef) of parameters calibrated in each stage.
    stageDInd{1} = []; %Indices of all parameters to be calibrated in future stages
    stageNInd{1} = []; %Indices of all parameters not used during current calibration stage
    stagePInd{1} = []; %Indices of all parameters calibrated in previous stage (cumulative)
else %Calibration conducted in multiple stages
    blObs = cell(nStage, 1);
    
    for ss = 1 : nStage
        %Find parameters for each stage:
        stageTemp = cell(0,1);
        for jj = 1 : numel(sOpt.stageVar(ss,:))
            indVarUse = cellfun(@(x) ~isempty(regexpi(x,sOpt.stageVar{ss,jj})), sMeta.coef(:,6));
            stageTemp = [stageTemp; sMeta.coef(indVarUse,1)];
        end
        clear jj
        stagePrmNm{ss} = stageTemp;
        
        %Find all parameter indices not in current stage:
        for jj = 1 : numel(sMeta.coef(:,1))
            if any(strcmpi(sMeta.coef(jj,1),stagePrmNm{ss}(:)))
                stagePrmInd{ss} = [stagePrmInd{ss}; jj];
            else
                stageNInd{ss} = [stageNInd{ss}; jj];
            end
        end
        clear jj
        
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
            clear ll
        end
        
        stageDInd{ss} = setdiff(stageNInd{ss}, stagePInd{ss});
        
        if ~iscell(sObs)
            error('CCHFcalibrate:obsNotCell','The observation input is expected to be of type cell.');
        end
        nReg = numel(sObs(:));
        
        ptsObs = cell(nReg, 1);
        blObsUse = ptsObs;
        stageObsStruc{ss} = cell(nReg, 1);
        blObs{ss} = ones(nReg, 1);
        for kk = 1 : nReg
            %Find observation data for each stage:
            ptsObs{kk} = fieldnames(sObs{kk});
            nPts = numel(ptsObs{kk}(:));
            
%             blNoObs{ss}{kk} = zeros(nPts, 1);
            blObsUse{kk} = zeros(numel(ptsObs),1);
            for jj = 1 : nPts
                typKeepCurr = cell(0,1);
                
                fieldsCurrKp = fieldnames(sObs{kk}.(ptsObs{kk}{jj}));
                %Remove metadata fields:
                for ll = numel(fieldsCurrKp(:)) : -1 : 1
                   if regexpbl(fieldsCurrKp{ll}, {'date', 'path', 'type','flag','att','lat','lon'})
                       fieldsCurrKp(ll) = [];
                   end
                end
                clear ll

                %Loop over each observation type:
                for ll = 1 : numel(fieldsCurrKp(:))
                    %Four categories are 'input', 'cryo', 'land', 'routing'
                    if any(strcmpi(fieldsCurrKp{ll},{'flow','flowrate','streamflow'}))
                        typeCurr = 'routing';
                    elseif any(strcmpi(fieldsCurrKp{ll},{'stake','swe','snowradar','sca','casi','snowcover','geodetic','mb','massbalance'}))
                        typeCurr = 'cryo';
                    elseif any(strcmpi(fieldsCurrKp{ll},{'tmp','tas','tmx','tmn','pre','pr'}))
                        typeCurr = 'input';
                    else
                        warning('CCHF_calibrate:unkownObs',[char(39) ...
                            fieldsCurrKp{1} char(39) ' observation data has not '...
                            'been programmed for in the calibration routine. '...
                            'Please add.']);
                    end

                    if any(strcmpi(typeCurr, sOpt.stageVar(ss,:)))   
                        blObsUse{kk}(ll) = 1;
                        typKeepCurr{end+1} = fieldsCurrKp{ll};
                    end  
                end
                clear ll
                
                %Assign fields to observations for current stage:
                fieldsCurrAll = fieldnames(sObs{kk}.(ptsObs{kk}{jj}));
                if numel(typKeepCurr(:)) > 0
                    indUse = nan(0,1);
                    for ll = 1 : numel(typKeepCurr(:))
                        for zz = 1 : numel(fieldsCurrAll(:))
                            if regexpbl(fieldsCurrAll{zz}, typKeepCurr{ll})
                                indUse(end+1) = zz;
                            end
                        end
                        clear zz
                    end
                    clear ll
                    
                    if numel(indUse) > 0
                        indLat = find(strcmpi(fieldsCurrAll, 'lat') == 1);
                        if ~isempty(indLat)
                            indUse(end+1:end+numel(indLat)) = indLat;
                        end
                        indLat = find(strcmpi(fieldsCurrAll, 'latitude') == 1);
                        if ~isempty(indLat)
                            indUse(end+1:end+numel(indLat)) = indLat;
                        end
                        indLon = find(strcmpi(fieldsCurrAll, 'lon') == 1);
                        if ~isempty(indLon)
                            indUse(end+1:end+numel(indLon)) = indLon;
                        end
                        indLon = find(strcmpi(fieldsCurrAll, 'longitude') == 1);
                        if ~isempty(indLon)
                            indUse(end+1:end+numel(indLon)) = indLon;
                        end
                    end
                    
                    %Assign to array:
                    for ll = 1 : numel(indUse)
                        stageObsStruc{ss}{kk}.(ptsObs{kk}{jj}).(fieldsCurrAll{indUse(ll)}) ...
                            = sObs{kk}.(ptsObs{kk}{jj}).(fieldsCurrAll{indUse(ll)});
                    end
                    clear ll
%                 else
%                     disp(['no observations to use in stage ' num2str(ss) ', region ' num2str(kk) ', point ' num2str(jj) '.']);
                end
                 
            end
            clear jj

            if isempty(stageObsStruc{ss}{kk})
                blObs{ss}(kk) = 0;
            end
                
%             indUse = find(blUse == 1);
%             for jj = 1 : numel(indUse)
%                 stageObsStruc{ss}.(ptsObs{indUse(jj)}) = sObs.(ptsObs{indUse(jj)});
%             end
%             clear jj
        end
        clear kk
        
        if all(blObs{ss} == 0)
            error('CCHF_calibrate:obsCurrPt',['No observation data '...
                'were found for stage ' num2str(ss) ' of calibration point.']);
        end
    end
    clear ss
end


%Edit parameter options if any parameters are constant
indCfConst = [];
valCfConst = [];
for ii = 1 : numel(sMeta.coef(:,1))
    if sMeta.coef{ii,2} == sMeta.coef{ii,3} && sMeta.coef{ii,2} == sMeta.coef{ii,4}
        indCfConst(end+1) = ii;
        valCfConst(end+1) = sMeta.coef{ii,2};
    end
end

%Format constant parameters for use in stage-indice format (and remove from
%arrays of indices to be calibrated)
stageCInd = cell(nStage, 1);
valCInd   = cell(nStage, 1);
if ~isempty(indCfConst)
    for ii = 1 : nStage
        stageCInd{ii} = indCfConst;
        valCInd{ii}   = valCfConst;
        
        stagePrmInd{ii} = setdiff(stagePrmInd{ii}, stageCInd{ii}); 
        stageDInd{ii}   = setdiff(  stageDInd{ii}, stageCInd{ii}); %Indices of all parameters to be calibrated in future stages
        stageNInd{ii}   = union(    stageNInd{ii}, stageCInd{ii}); %Indices of all parameters not used during current calibration stage
        stagePInd{ii}   = setdiff(  stagePInd{ii}, stageCInd{ii}); %Indices of all parameters calibrated in previous stage (cumulative)
    end
end



%%ALLOW OPTION TO RESUME PREVIOUSLY INTERUPPTED CALIBRATION
%Requires loading previously tested parameter sets from files
blFresh = 1; indStgStrt = 1;
if regexpbl(sMeta.runType,{'calib','resume'},'and')
    if isfield(sPath{1},'resume')
        filesParam = dir(fullfile(sPath{1}.resume, '*.csv'));
        filesParam = extractfield(filesParam, 'name');

        if ~isempty(filesParam)
            if numel(filesParam) > 1
                %Remove parameter files that don't have 'all' in the name
                for ii = numel(filesParam) : -1 : 1
                    if ~regexpbl(filesParam{ii}, 'all')
                        filesParam(ii) = [];
                    end
                end
            end
            
            indStage = nan(numel(filesParam),1);
            for ii = 1 : numel(filesParam)
                indUnd = regexpi(filesParam{ii}, '_');
                indStage(ii) = str2double(filesParam{ii}(indUnd(end)+1:end-4));
            end
            clear ii
            
            %Determine which stage to start at (based on parameter files in
            %folder)
            indStgStrt = intersect((1:nStage),indStage);
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
                    filePrm = fullfile(sPath{1}.resume, filesParam{indStage == indStgLd(ll)});
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
                        clear ii
    
                        for ii = nGen : -1 : 1
                            if all(isnan(coefFitness(ii,:)))
                                coefFitness = coefFitness(1 : ii-1,:);
                                coefFamily = coefFamily(1 : ii-1,:,:);
                            end
                        end
                        clear ii


                        %If stage is not 1st, find the best performing
                        %parameter set from the previous stage (this is used
                        %later)
                        if numel(indStgLd) > 1 && ll < numel(indStgLd)
                            [~, rMin, cMin] = min2d(coefFitness);

                            stageBestPrm{ll} = num2cell(squeeze(coefFamily(rMin, cMin,:)));
                        end
                        
                        indGenStrt = ones(nStage,1);
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
                clear ll
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
    indGenStrt = ones(nStage,1);
end






%Check that all indices have been located
checkInd = 0;
for ss = 1 : nStage
   checkInd = checkInd + numel(stagePrmInd{ss});
end
clear ss
if checkInd + numel(indCfConst) ~= numel(sMeta.coef(:,1))
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
    clear ii
end


%Create path and header for writing parameter sets:
%Set output path and header for writing coefficients:
pathGenParam = fullfile(sPath{1}.outputMain,['parameter_each_gen' '.txt']);
pathBstCoef = fullfile(sPath{1}.outputMain,['coefficients_' 'fittest' '.txt']);

disp(['All output from model calibration will be written to:' char(39) sPath{1}.outputMain]);

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
for ss = indStgStrt : nStage
    %Display messgage about start of current stage:
    if numel(stagePrmInd{ss}) == 0
        disp(['0 parameters have been identified to calibrate this '...
            'stage. Therefore, this stage will be skipped.']);
        continue
    else
        if blFresh
            disp(['Calibration stage ' num2str(ss) ' of ' ...
                num2str(nStage) ' starting.' char(10) ...
                num2str(numel(stagePrmInd{ss})) ' of ' ...
                num2str(numel(sMeta.coef(:,1))) ' parameters are being '...
                'calibrated during this stage.']);
        else
            disp(['Calibration stage ' num2str(ss) ' of ' ...
                num2str(nStage) ' resuming (starting with stage ' num2str(ss) ', iteration '...
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
        %Set parameters being calibrated in future stages to default values:
        if ~isempty(stageDInd{ss})
            for ll = 1 : nGen
                for mm = 1 : nPop
                    coefFamily(ll, mm, stageDInd{ss}) = cell2mat(sMeta.coef(stageDInd{ss},4));
                end
                clear mm
            end
            clear ll
        end
        
        %Set parameters that were previously calibrated to the found values:
        if ~isempty(stagePInd{ss}) && ss ~= 1
            for ll = 1 : nGen
                for mm = 1 : nPop
                    coefFamily(ll, mm, stagePInd{ss}) = cell2mat(stageBestPrm{ss-1}(stagePInd{ss}));
                end
                clear mm
            end
            clear ll
        end
        
        %Set parameters that are assigned constant values:
        if ~isempty(stageCInd{ss})
            for ll = 1 : nGen
                for mm = 1 : nPop
                    coefFamily(ll, mm, stageCInd{ss}) = valCInd{ss};
                end
                clear mm
            end
            clear ll
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
    

    %LOOP OVER GENERATIONS OF CALIBRATION
    valMostFit = nan(nGen,1);
    maxGen = 0;
    
    for ii = indGenStrt(ss) : nGen
        %Calculate coefficients for current generation:
        if ii ~= 1 && ~regexpbl(sOpt.type, 'uniform') %Don't optimize if uniform sampling used because all parameter sets pre-determined.
            coefFamily(ii,:,stagePrmInd{ss}) ...
                = coef_next_gen(coefFamily(:,:,stagePrmInd{ss}), coefFitness(:,:), ii-1, prmBnds(stagePrmInd{ss},:), sOpt);
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
        fitCurr = nan(numel(indRun), 1);
        %Loop over each parameter set combination (in parallel if possible)
        shortCnt = 0;

        %Turn off warning during parloop
%         warning('off','all');
%         if ~isempty(gcp('nocreate'))
%             pctRunOnAll warning('off','all')
% %                 spmd
% %                   warning('off','all')
% %                 end
%         end

        parfor jj = 1 : numel(indRun)
            %Turn off warning during parloop
            warning('off','all')
            if ~isempty(gcp('nocreate'))
                pctRunOnAll warning('off','all')
    %                 spmd
    %                   warning('off','all')
    %                 end
            end
            warning('off', 'mod_v_obs:unknownEvalType');
            fitCurrSites = nan(nSites, 1);
            
            if all(~isnan(dateEval))
                [sModTemp, fitTemp] = CCHF_engine_v4(sPath, sHydro, sMeta, 'cf', prmCurr(jj,:)', 'obs', obsCurr, 'evalatt', shortEvalAtt);
                if isnan(fitTemp)
                    for mm = 1 : nSites
                        fitCurrSites(mm) = nanmean(cell2mat(mod_v_obs_v2(obsCurr{mm}, sModTemp{mm}, fitTest, 'combineType')));
                    end
                else
                    shortCnt = shortCnt + 1;
                    fitCurrSites = fitTemp(:)';
                end
            else
                sModTemp = CCHF_engine_v4(sPath, sHydro, sMeta, 'cf', prmCurr(jj,:)');
                for mm = 1 : nSites
                    fitCurrSites(mm) = nanmean(cell2mat(mod_v_obs_v2(obsCurr{mm}, sModTemp{mm}, fitTest, 'combineType')));
                end
            end
            
            switch siteCombine
                case 'linear'
                    fitCurr(jj) = mean(fitCurrSites);
                case 'square'
                    fitCurr(jj) = sqrt(sum(fitCurrSites.^2));
            end
        end
        clear jj
        
        %Turn warnings back on
        warning('on', 'mod_v_obs:unknownEvalType');
        warning('on','all')
        if ~isempty(gcp('nocreate'))
            pctRunOnAll warning('on','all') 
            spmd
              warning('on','all')
            end
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
            clear zz
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
            ['Stage: ' num2str(ss) ' of ' num2str(nStage)];...
            ['Generation: ' num2str(ii)];...
            [sOpt.fitTest ': ' num2str(valMostFit(ii))]};
        if ss == 1 && ii == 1
            write_CCHF_coef(pathGenParam, [hdrFull(:); blanks(0); hdrShortCoef(:)], [sMeta.coef(:,1), num2cell(squeeze(coefFamily(ii,memFit,:)))]); 
        else
            write_CCHF_coef(pathGenParam, hdrShortCoef, [sMeta.coef(:,1), num2cell(squeeze(coefFamily(ii,memFit,:)))],'permissions','a'); 
        end
        
        %Record best set of parameters from all generations in current
        %stage:
        if ii == 1
            indGenFit = 1;
            [~, indPopFit] = nanmin(coefFitness(1,:));
        else
            [~, indGenFit, indPopFit] = min2d(coefFitness(1:ii,:));
        end
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
            pathCalib = fullfile(sPath{1}.outputMain,['all_parameter_sets-stage_' num2str(ss) '.csv']);
            fParam = fopen(pathCalib,'w+');
            
            numColumns = size(coefFamily,3)+1;
            headerFmt = repmat('%s,',1,numColumns-1);
            numFmt = repmat('%f,',1,numColumns-1);
            hdr = [sMeta.coef(:,1)', sOpt.fitTest];

            fprintf(fParam,[headerFmt,'%s\n'], hdr{:});
            fclose(fParam);
        elseif blFresh == 0 && ss == indStgStrt && ii == indGenStrt(ss)
            pathCalib = fullfile(sPath{1}.outputMain,['all_parameter_sets-stage_' num2str(ss) '.csv']);
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
        %Record highest generation reached
        maxGen = ii;
    end %End of generation loop
    clear ii
    
    
    %Record best value to display
    if regexpbl(sOpt.fitTest,{'NSE','KGE','Pearson'})
        [~, indGenFit] = max(valMostFit(1:maxGen));
    else
        [~, indGenFit] = min(valMostFit(1:maxGen));
    end
    
    %Stage completion message:
    disp(['CCHF model calibration stage ' num2str(ss) ' of ' num2str(nStage) ' has completed.' ...
        char(10) 'The fittest set of parameter coefficients tested yielded an ' sOpt.fitTest ...
        ' value of ' num2str(round2(valMostFit(indGenFit),2)) ...
        ', ' char(10) 'which occured during iteration ' num2str(indGenFit) ...
        ' of ' num2str(nGen) ' total iterations.' char(10)]);
    
    
    %RECORD N BEST PERFORMING PARAMETER SETS:
    %This is potentially useful for assessing equifinality
    bestPrmSets{ss} = reshape(coefFamily,[], numel(sMeta.coef(:,1))); %indices are: (generation, pop member, parameter)
    [bestPrfm, indRank] = sort(reshape(coefFitness,[],1));
    bestPrmSets{ss} = bestPrmSets{ss}(indRank,:);
    
    %If number of best parameter options not specified, record 0.5% of all
    %test or 50, whichever is greater.
    if ~isfield(sOpt, 'nbest')
       sOpt.nbest = max(0.005*nGen*nPop, 50); 
    end
    bestPrmSets{ss} = bestPrmSets{ss}(1:sOpt.nbest,:);
    bestPrfm = bestPrfm(1:sOpt.nbest);
    
    pathBest = fullfile(sPath{1}.outputMain,[num2str(sOpt.nbest) '_parameter_sets-stage_' num2str(ss) '.csv']);
    fParamBest = fopen(pathBest,'w+');

    numColumnsBest = size(coefFamily,3)+1;
    headerFmtBest = repmat('%s,',1,numColumnsBest-1);
    numFmtBest = repmat('%f,',1,numColumnsBest-1);
    hdrBest = [sMeta.coef(:,1)', sOpt.fitTest];

    fprintf(fParamBest,[headerFmtBest,'%s\n'], hdrBest{:});
    paramWriteBest = [bestPrmSets{ss}, bestPrfm]; 
    fprintf(fParamBest,[numFmtBest,'%f\n'], paramWriteBest');
    fclose(fParamBest);
    
    %Calculate and record correlation between each parameters in N best sets:
    threshSig = 0.01; %correlation must be significant to less than this value
    prmCorr = nan(numel(bestPrmSets{ss}(1,:)));
    prmSig = nan(numel(bestPrmSets{ss}(1,:)));
    for mm = 1 : numel(bestPrmSets{ss}(1,:))
        for nn = 1 : numel(bestPrmSets{ss}(1,:))
            [r,p] = corrcoef(bestPrmSets{ss}(:,mm), bestPrmSets{ss}(:,nn));
            prmCorr(mm,nn) = r(2);
            prmSig(mm,nn) = p(2);
        end
        clear nn
    end
    clear mm
    prmCorr(prmSig > threshSig) = nan;

    %Write to file:
    pathCorr = fullfile(sPath{1}.outputMain, ['Pearson_corr-' num2str(sOpt.nbest) '_parameter_sets-stage_' num2str(ss) '.csv']);
    fCorr = fopen(pathCorr,'w+');

    nCol = numel(bestPrmSets{ss}(1,:)) + 1;
    headerFmt = repmat('%s,',1,nCol-1);
    numFmt = repmat('%f,',1,nCol-2);

    hdrCorr = [blanks(1), sMeta.coefAtt(:,1)'];

    fprintf(fCorr, [headerFmt,'%s\n'], hdrCorr{:});
    for mm = 1 : numel(bestPrmSets{ss}(1,:))
        fprintf(fCorr, '%s,', sMeta.coefAtt{mm,1});
        fprintf(fCorr, [numFmt,'%f\n'], prmCorr(mm,:));
    end
    clear mm
    
    fclose(fCorr);
    
    blFresh = 1; %Reset blFresh if set to 0 using resume section of code
end
clear ss


%Close dedicated workers used for parallel processing:
poolobj = gcp('nocreate');
delete(poolobj);


%Rerun model for best parameter set:
sModOut = CCHF_engine_v4(sPath, sHydro, sMeta, 'cf', stageBestPrm{end});


%Write and Plot parameters:
plot_cal_prm(bestPrmSets, stagePrmInd, prmBnds, sMeta.coef, pathBstCoef);
hdrBst = [hdrFull(1:4,:); ...
    [sOpt.fitTest ' fitness score: ' num2str(valMostFit(indGenFit))]; ...
    ['Date: ' strNow]; ...
    hdrFull(5:end,:)];
write_CCHF_coef(pathBstCoef, hdrBst, [sMeta.coef(:,1), stageBestPrm{end}],'permissions','w');  

outParam = [sMeta.coef(:,1), stageBestPrm{end}, sMeta.coef(:,5:end)];