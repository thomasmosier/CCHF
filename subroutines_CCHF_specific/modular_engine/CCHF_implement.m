function varargout = CCHF_implement(sMeta, sOpt)

if iscell(sMeta.region)
    nSites = numel(sMeta.region(:));
elseif ischar(sMeta.region)
    nSites = 1;
    sMeta.region = {sMeta.region};
else
    error('CCHF_main:nmRegions',['The model dies not recognize the '...
        'format of the region input variable.']);
end

sMeta.('rtDir')    = cell(nSites, 1);
sMeta.('output')   = cell(nSites, 1);

sObs = cell(nSites, 1);

sHydro = cell(nSites, 1);

% if nSites ~= numel(regions(:))
%    error('CCHF_main:diffNumRegions',['The model is set to evaluate for ']) 
% end


varLat = 'latitude';
varLon = 'longitude';
%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
cd(pathScript);
addpath(genpath(pathScript));

%Initialize structure arrays for passing variables to downscaling script:
clear global sHydro sAtm sCryo sLand %Possibly simply initialize instead of clearing

%Create local variables used only within this script:
%Create output directory for model run output because
%multi-validation processing can eat up a lot of memory
%otherwise
varValOut = 'multival';
valRt = 'output_';
            

%'sMeta.useprevrun' provides option to load output from previous run (not normal
%condition)
if sMeta.useprevrun == 0
    
    %Either load path inputs or initialize arrays for assignment
    if ~isempty(sMeta.pathinputs)
        load(sMeta.pathinputs)
    else
        sPath = cell(nSites, 1);
        pathGage = cell(nSites, 1);
    end

    %Create time vector to loop over based on start date and time-step
    %If this field populated here, it wont be populated inside each model call,
    %which saves time
    sMeta = dates_run(sMeta, 'spin');
    sMeta.progress = 'year';

    %Load global constants from funtion:
    sMeta.global = global_params; %contains global parameter values (albedo of ice, etc.)


    %LOCATE INPUTS FOR EACH SITE:
    startPath = pwd;
    %Loop over sites:
    for ii = 1 : nSites
        %Only open UI if inputs paths not already known
        if isempty(sMeta.pathinputs)
            %Initialize structure for current site
            sPath{ii} = struct;
                sPath{ii}.('dem')    = blanks(0);
                sPath{ii}.('pr')     = blanks(0);
                sPath{ii}.('tas')    = blanks(0);
                sPath{ii}.('tasmin') = blanks(0);
                sPath{ii}.('tasmax') = blanks(0);
                sPath{ii}.('regclip')= blanks(0);
                sPath{ii}.('output') = blanks(0);
                sPath{ii}.('coef')   = blanks(0);
            sHydro{ii} = struct;

            %Ask to load information from directory of calibration run that was
            %interrupted
            if ii == 1 && regexpbl(sMeta.runType,{'calib','resume'},'and')
                uiwait(msgbox(sprintf(['Select the folder containing containing '...
                    'the unfinished calibration for ' sMeta.region{ii} '.\n']), ...
                    '(Click OK to Proceed)','modal'));
                sPath{ii}.resume = uigetdir(startPath,['Select the folder containing '...
                    'the unfinished calibration for ' sMeta.region{ii}]);
            end

            %Digital Elevation Model (DEM) selection and display:
            uiwait(msgbox(sprintf(['Select the Digital Elevation Model (DEM) for ' ...
                sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
            [fileDem, foldDem] = uigetfile({'*.asc';'*.txt'},['Select the Digital Elevation Model '...
                '(DEM) for ' sMeta.region{ii}], startPath);
            sPath{ii}.dem = fullfile(foldDem, fileDem);
            disp([char(39) sPath{ii}.dem char(39) ' has been chosen as the DEM.']);

            %Update search path:
            [startPath, ~, ~] = fileparts(sPath{ii}.dem);

            if isempty(fileDem) || isequal(fileDem, 0)
               error('CCHF_main:noDEM','A DEM has not been selected. Therefore, the program is aborting.'); 
            end

            %Load manual flow direction grid (ESRI ASCII Format):
            if sMeta.ldFdr == 1
                %Flow direction selection and display:
                uiwait(msgbox(sprintf(['Select the flow direction grid (typically created using ArcGIS) for ' ...
                    sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
                [fileFdr, foldFdr] = uigetfile({'*.asc';'*.txt'},['Select the flow direction grid for ' ...
                    sMeta.region{ii}], startPath);
                sPath{ii}.fdr = fullfile(foldFdr, fileFdr);
                disp([char(39) sPath{ii}.fdr char(39) ' has been chosen as the flow direction grid.']);

                if isempty(fileFdr) || isequal(fileFdr, 0)
                   error('CCHF_main:noFDR',['No flow direction grid '...
                       'has been selected, even though the option was chosen.' ...
                       ' Therefore, the program is aborting.']); 
                end
            end


            %FIND GLACIER DATA PATH:
            if regexpbl(sMeta.iceGrid, 'fine') && isempty(find_att(sMeta.module, 'glacier','no_warning'))
                disp('The glacier grid resolution is fine but glaciers are being treated as static.');
            end

            if regexpbl(sMeta.iceGrid, 'same')
                uiwait(msgbox(sprintf(['Select the glacier presence '...
                    'geo-referenced grid or shapefile for ' ...
                    sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
                [fileGlac, foldGlac] = uigetfile({'*.*'},['Select the binary glacier presence grid for ' ...
                    sMeta.region{ii}], startPath);
                sPath{ii}.ice = fullfile(foldGlac, fileGlac);
                disp([char(39) sPath{ii}.ice char(39) ' has been chosen as the binary glacier presence grid.']);

                if isempty(fileGlac) || isequal(fileGlac, 0)
                   error('CCHF_main:noGlacier',['No glacier presence grid '...
                       'has been selected, even though the option was chosen.' ...
                       ' Therefore, the program is aborting.']); 
                end
            elseif regexpbl(sMeta.iceGrid, 'fine')
                iceAuxDisp = ['on a finer grid that will determine the ' ...
                    'spatial scale of glacier process evaluation'];

                uiwait(msgbox(sprintf(['Select the glacier presence geo-referenced grid or shapefile for ' ...
                    sMeta.region{ii} ' ' iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
                [fileGlac, foldGlac] = uigetfile({'*.*'},['Select the glacier presence geo-referenced grid or shapefile for ' ...
                    sMeta.region{ii} ' ' iceAuxDisp '.\n'], startPath);
                sPath{ii}.ice = fullfile(foldGlac, fileGlac);
                disp([char(39) sPath{ii}.ice char(39) ' has been chosen as the glacier presence grid.']);

                if isempty(fileGlac) || isequal(fileGlac, 0)
                   error('CCHF_main:noGlacier',['No glacier presence grid '...
                       'has been selected, even though the option was chosen.' ...
                       ' Therefore, the program is aborting.']); 
                end

                uiwait(msgbox(sprintf(['Select the glacier surface DEM for ' ...
                    sMeta.region{ii} ' ' iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
                [fileGlacDem, foldGlacDem] = uigetfile({'*.*'},['Select the glacier surface DEM for ' ...
                    sMeta.region{ii} ' ' iceAuxDisp '.\n'], startPath);
                sPath{ii}.iceDem = fullfile(foldGlacDem, fileGlacDem);
                disp([char(39) sPath{ii}.ice char(39) ' has been chosen as the glacier DEM grid.']);

                if isempty(fileGlacDem) || isequal(fileGlacDem, 0)
                   error('CCHF_main:noGlacierDem',['No glacier DEM grid '...
                       'has been selected, even though the option was chosen.' ...
                       ' Therefore, the program is aborting.']); 
                end
            elseif ~regexpbl(sMeta.iceGrid, 'none')
                error('CCHF_main:unknownIceGridType',['The ice grid method is ' sMeta.iceGrid ', which is not known.']);
            end


            %Load debris cover grid (if single value, assumed uniform thickness):
            if sMeta.blDebris && ~regexpbl(sMeta.iceGrid, 'none')
                uiwait(msgbox(sprintf(['Select the debris cover thickness grid (meters depth) for ' ...
                    sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
                [fileDeb, foldDeb] = uigetfile({'*.asc';'*.txt'},['Select the debris cover thickness grid for ' ...
                    sMeta.region{ii}], startPath);
                sPath{ii}.icdbr = fullfile(foldDeb, fileDeb);
                disp([char(39) sPath{ii}.icdbr char(39) ' has been chosen as the debris cover thickness grid.']);

                if isempty(fileDeb) || isequal(fileDeb, 0)
                   error('CCHF_main:noDebris',['No glacier debris cover thickness grid '...
                       'has been selected, even though the option was chosen.' ...
                       ' Therefore, the program is aborting.']); 
                end
            end

            %UI to locate ice thickness grid (if specified)
            iceWEMod = find_att(sMeta.module, 'glacier0');
            if strcmpi(iceWEMod, 'external') && ~regexpbl(sMeta.iceGrid, 'none')
                uiwait(msgbox(sprintf(['Select the ice thickness grid (expressed as water equivalent, meters depth) for ' ...
                    sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
                [fileIcWe, foldIcWe] = uigetfile({'*.asc';'*.txt'},['Select the ice thickness grid (water equivalent) for ' ...
                    sMeta.region{ii}], startPath);
                sPath{ii}.icwe = fullfile(foldIcWe, fileIcWe);
                disp([char(39) sPath{ii}.icwe char(39) ' has been chosen as the ice thickness grid (expressed as water equivalent).']);

                if isempty(fileIcWe) || isequal(fileIcWe, 0)
                   error('CCHF_main:noIceWaterEq',['No glacier thickness grid '...
                       'has been selected, even though the option was chosen.' ...
                       ' Therefore, the program is aborting.']); 
                end
            end

            %Load ice lake fraction
            if sMeta.blIcePond == 1  && ~regexpbl(sMeta.iceGrid, 'none')
                uiwait(msgbox(sprintf(['Select the glacier pond fraction grid (range = 0 to 1) for ' ...
                    sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
                [fileIcePnd, foldIcePnd] = uigetfile({'*.asc';'*.txt'},['Select the glacier pond fraction grid for ' ...
                    sMeta.region{ii}], startPath);
                sPath{ii}.icpndx = fullfile(foldIcePnd, fileIcePnd);
                disp([char(39) sPath{ii}.icpndx char(39) ' has been chosen as the glacier pond fraction grid.']);

                if isempty(fileIcePnd) || isequal(fileIcePnd, 0)
                   error('CCHF_main:noIceLake',['No glacier pond fraction grid '...
                       'has been selected, even though the option was chosen.' ...
                       ' Therefore, the program is aborting.']); 
                end
            end


    %         %If a binary grid is being used to clip the region:
    %         if blClip == 1
    %             uiwait(msgbox(sprintf(['Select the binary region grid that will be used to clip the ' ...
    %                 sMeta.region{ii} ' region.\n']), '(Click OK to Proceed)','modal'));
    %             [fileRegClip, foldRegClip] = uigetfile({'*.asc';'*.txt'},['Select the binary region clip grid for ' ...
    %                 sMeta.region{ii}], startPath);
    %             sPath{ii}.regionClip = fullfile(foldRegClip, fileRegClip);
    %             disp([char(39) sPath{ii}.regionClip char(39) ' has been chosen as the binary region clip grid.']);
    % 
    %             if isempty(fileRegClip) || isequal(fileRegClip, 0)
    %                error('CCHF_main:noClip',['No region clipping file '...
    %                    'has been selected, even though the option was chosen.' ...
    %                    ' Therefore, the program is aborting.']); 
    %             end
    %         end


            %Identify climate variables to load:
            [sMeta.varLd, sMeta.varLdDisp] = clm_var_load(sMeta.module);
            %Load climate variables:
            sPath{ii} = clm_path_ui(sPath{ii}, sMeta, sMeta.region{ii});
        end %End inputs UIs
        

        %Create unique 'main' output directory based upon inputs:
        if ii == 1   
            if isfield(sPath{ii},'resume')
                foldOutputMain = sPath{ii}.resume;
                [~, sMeta.strModule] = CCHF_out_dir(sPath{ii}, sMeta);
            else
                [foldOutputMain, sMeta.strModule] = CCHF_out_dir(sPath{ii}, sMeta);
            end
            
            if~exist(foldOutputMain, 'dir')
               mkdir(foldOutputMain); 
            end
        end
        sPath{ii}.outputMain = foldOutputMain;
        sPath{ii}.output = fullfile(sPath{ii}.outputMain, sMeta.region{ii});


        %Load DEM and calculate watershed geometry fields:
        sDem = read_geodata_v2(sPath{ii}.dem, 'data', nan(1,2), nan(1,2), nan(1,2), '0', 'out', 'none', 'onefile');
        sHydro{ii}.dem = sDem.data;
        sHydro{ii}.(varLat) = sDem.(varLat);
        sHydro{ii}.(varLon) = sDem.(varLon);
        if isfield(sPath{ii}, 'fdr')
            sFdr = read_geodata_v2(sPath{ii}.fdr, 'data', nan(1,2), nan(1,2), nan(1,2), '0', 'out', 'none', 'onefile');
            sHydro{ii}.fdrESRI = sFdr.data;
            clear('sFdr');
        else
            warning('cchfMain:noFDR',['CCHF currently does not have the '...
                'capacity to calculate a flow direction grid. Instead, '...
                'this must be generated in a standalone program such as ArcMap.']);
        end
        clear('sDem');

        %Load observation data to use for calibration or validation:
        if regexpbl(sMeta.runType,{'calibrat','validat','default'}) && isempty(sMeta.pathinputs)
            %Calibration/validation data required:
            uiwait(msgbox(sprintf(['Select the file(s) with observation data for ' ...
                sMeta.region{ii} '.  If data does not include information such as location, you will '...
                'be required to enter that in the main panel.']), ...
                '(Click OK to Proceed)','modal'));

            pathGage{ii} = cell(sMeta.nGage,1);
            for jj = 1 : sMeta.nGage
                [fileGage, foldGage] = uigetfile({'*.*'}, ...
                    ['Select observation data ' num2str(jj) ' of ' ...
                    num2str(sMeta.nGage) ' for ' sMeta.region{ii}], startPath);
                pathGage{ii}{jj} = fullfile(foldGage, fileGage);

                disp([pathGage{ii}{jj} ' has been selected as observation data file ' num2str(jj) ' of ' num2str(sMeta.nGage) ' for ' sMeta.region{ii} '.']);
            end
        elseif ~regexpbl(sMeta.runType,{'calibrat','validat','default'}) && ~isempty(sMeta.pathinputs)
            pathGage{ii} = cell(0,1);
        end
    end %End of loop to select sites 
    clear ii


    %START PROCESSING LOG to record all messages displayed on screen to file in output directory
    [fileDiary] = command_log(sPath{1}.outputMain);
    disp(['All displayed text is being written to ' char(39) fileDiary ...
        char(39) '.']);


    %display message for location of output folder:
    if numel(nSites) > 1
       uiwait(msgbox(sprintf(['Model outputs will be written to a subfolder of ' ...
           sMeta.regions{1} ' because this region is listed first.\n']), ...
           '(Click OK to Proceed)','modal'));
    end


    %GET PARAMETER FILE (IF VALIDATION OR SIMULATION RUN)
    if regexpbl(sMeta.runType,{'valid','sim'})
        %Load parameters from file:
        uiwait(msgbox(sprintf(['Select the set of parameter coefficients ' ...
            'written during a calibration run of the CCHF model for ' ...
            sMeta.region{1} '.\n']), '(Click OK to Proceed)','modal'));
        [fileCoef, foldCoef] = uigetfile({'*.txt'; '*.csv'}, ['Select the set of parameter coefficients for ' ...
            sMeta.region{1}], startPath);
        if isempty(fileCoef) || isempty(foldCoef)
            error('CCHF_main:noCoefFile',['No file containing a set of '...
                'model coefficients has been selected. This is a requirement.' ...
               ' Therefore the program is aborting.']); 
        end
        
        for ii = 1 : nSites
            sPath{ii}.coef = fullfile(foldCoef, fileCoef);
        end
        clear ii
        disp([char(39) sPath{1}.coef char(39) ' has been chosen as the set of parameter coefficients.']);
    end
    
    
    %LOAD ALL OBSERVATION DATA
    for ii = 1 : nSites
        sMeta.siteCurr = ii;

        if regexpbl(sMeta.runType,{'calibrat','validat','default'})
            if ischar(pathGage{ii})
                pathGage{ii} = {pathGage{ii}};
            end  
            
            %Remove any duplicate or empty elements:
            pathGage{ii} = unique(pathGage{ii}(~cellfun('isempty',pathGage{ii})));  
            
            for zz = numel(pathGage{ii}(:)) : -1 : 1
                if isnumeric(pathGage{ii}{zz}) || isempty(regexpi(pathGage{ii}{zz}, '[a-z]'))
                    pathGage{ii}(zz) = [];  
                end
            end
            
            %Load observation data:
            [sObs{ii}, pathCchfGage] = read_gagedata(pathGage{ii}, sHydro{ii}.(varLon), sHydro{ii}.(varLat), ...
                'time',[sMeta.dateStart;sMeta.dateEnd], ...
                'mask',sHydro{ii}.dem);
            
            %Save path to CCHF formatted obersvation data (to save and use
            %in future runs)
            if ~isempty(pathCchfGage)
                pathGage{ii} = pathCchfGage;
            end
        else
            sObs{ii} = struct;
        end

        if isempty(sMeta.output{ii})
           sMeta.output{ii} = sMeta.addoutput;
        end
        [sMeta.output{ii}, sObs{ii}] = output_all_gage(sMeta.output{ii}, sObs{ii}, sHydro{ii}.(varLon), sHydro{ii}.(varLat));
    end
    clear ii


    %Calculate watershed properties and load/initialize glacier data:
    for ii = 1 : nSites
        sMeta.siteCurr = ii;

%         if blClip == 1 && isfield(sPath{ii}, 'regionClip') 
%             sRegClip = read_geodata_v2(sPath{ii}.regionClip, 'data', nan(1,2), nan(1,2), nan(1,2), '0', 'out', 'none');
%             sHydro{ii}.regionClip = sRegClip.data;
%             clear('sRegClip');
%         end

        %Calculate fdr, fac, flow calculation order, slope, aspect, distance between neighboring centroids
        sHydro{ii} = watershed(sHydro{ii});

        %Initialize ice grid:
        sHydro{ii}.sIceInit = ice_grid_load_v2(sHydro{ii}, sPath{ii}, sMeta);
    end
    clear ii


    %%Find paremeters needed for current version of model:
    %Do this by running CCHF model in 'parameter' mode
    sMeta.mode = 'parameter';
    sMeta.coefAtt = CCHF_engine_v4(sPath, sHydro, sMeta);

    
    %IF VALIDATION, LOAD COEFFICIENT SET AND COMPARE TO PARAMETERS NEEDED FOR CURRENT MODEL:
    if regexpbl(sMeta.runType,{'valid','sim'})
        sMeta.coef = CCHF_ld_prms(sPath{1}.coef, sMeta);
    elseif ~regexpbl(sMeta.runType,'calibrate') && numel(sMeta.coefAtt(1,:)) >= 6
        %Reduce dimensions of sMeta.coef (loaded during 'parameter' run)
        sMeta.coef = sMeta.coefAtt(:,[1,4,5]); 
    end

    %Set paramters specified in main function header:
    if ~regexpbl(sMeta.runType,'calibrate') %Any run type besides calibration
        sMeta.coef = CCHF_set_prm_val(sMeta.coef, sMeta.fixedPrm);
    else %Calibration
        sMeta.coefAtt = CCHF_set_prm_val(sMeta.coefAtt, sMeta.fixedPrm);
    end



    %%DISPLAY RELEVENT CONTENT TO CONSOLE
    %Display modules chosen:
    disp('The chosen set of module representations is: ');
    for ii = 1 : numel(sMeta.module(:,1))
        disp([sMeta.module{ii,1} ' = ' sMeta.module{ii,2}]);
    end
    disp(blanks(1));

    %Display message with citation information:
    [~] = MHS_cite_v2('CCHF');

    %Display modeling package being used and user choices:
    [~,dFold1,dFold2] = fileparts(pwd);
    disp([char(39) dFold1 dFold2 char(39) ' is being used.' char(10)]);
    disp_CCHF_meta_v2(sPath, sMeta)  


    %Find individual files to load each timestep (saves time later)
    for ii = 1 : nSites
        sPath{ii} = path_find_files(sPath{ii}, sMeta);    
    %     %Record root directory in sMeta (for use during simulations)
    %     indOutRt = regexpi(sPath{ii}.output,filesep);
    %     sMeta.rtDir{ii} = sPath{ii}.output(1:indOutRt(end)-1);

        [sMeta.rtDir{ii}, ~, ~] = fileparts(sPath{ii}.dem);
    end
    clear ii

    %Save input paths in file for possible future use
    if isempty(sMeta.pathinputs)
        [dirInputs, ~, ~] = fileparts(sMeta.rtDir{1});
        fileInputs = 'CCHF_saved_input_paths_4';
        for ii = 1 : nSites
            fileInputs = [fileInputs '_' sMeta.region{ii}];
        end
        fileInputs = [fileInputs '.mat'];
        
        pathInputs = fullfile(dirInputs,fileInputs);
        save(pathInputs, 'sPath', 'pathGage');
        disp(['Paths of all input data used saved to ' pathInputs char(10) 'This path can be set in main script to supress UI for loading inputs. You can customize the file name and location.'])
    end
    
    
    %%RUN THE HYDROLOGIC MODEL:
    tStrt = now;
    % sOutput = cell(1,1);
    if regexpbl(sMeta.runType,'calib')
        sMeta.mode = 'calibrate';
        sMeta.coef = sMeta.coefAtt;
            sMeta = rmfield_x(sMeta,'progress');
        [sOutput, sMeta.coef] = CCHF_calibrate_v3(sObs, sOpt, sPath, sHydro, sMeta);
    elseif regexpbl(sMeta.runType,{'valid','sim'})
        sMeta.mode = 'validate';

        %Validate model using one of two different settings: 
        %(1) Multiple parameter sets or (2) single parameter set
        if iscell(sMeta.coef) && ~ischar(sMeta.coef{1}) && ~all2d(size(sMeta.coef{1}) == 1)
            nReg = numel(sMeta.region(:));
            pathValOut = cell(nReg, 1);
            for ii = 1 : nReg
                pathValOut{ii} = fullfile(sPath{ii}.output, 'model_output');
                if ~exist(pathValOut{ii}, 'dir')
                    mkdir(pathValOut{ii});
                end
            end
            
            sMeta = rmfield_x(sMeta,'progress');

            %Ensure no output data written to file:
            sMeta.wrtGridOut = 0;
            sMeta.blDispOut = 0;
            
            %Reformat all parameter sets:
            coefAll = nan([numel(sMeta.coef{1}(:,1)), numel(sMeta.coef(:,1))]);
            for ii = 1 : numel(coefAll(1,:))
                coefAll(:,ii) = cell2mat(sMeta.coef{ii,:}(:,2));
            end
            clear ii

            %Open "Matlab pool" for parallel validation:
            if isempty(gcp('nocreate'))
                localCluster = parcluster('local'); %Check number of possible "workers"
                if exist('maxWorkersVal', 'var') && localCluster.NumWorkers > sMeta.maxworker
                    maxWorkersCurr = sMeta.maxworker;
                else
                    maxWorkersCurr = localCluster.NumWorkers;
                end
                parpool(maxWorkersCurr); %Dedicate all available cores to parallel calibration
            end

            if ~isempty(gcp('nocreate'))
               pctRunOnAll warning('off', 'all') %Turn off warnings during parfor
            end

            sOutput = cell(numel(coefAll(1,:)),1);
            nRuns = numel(coefAll(:,1));
            parfor ii = 1 : numel(coefAll(1,:))
                if mod(ii,10) == 0
                   display(['Currently on validation run ' num2str(ii) ' of ' num2str(nRuns)]); 
                end

                %Write output to matlab file
                sOutput{ii} = fullfile(pathValOut, [valRt num2str(ii) '.mat']);
                
                [~] = CCHF_engine_v4(sPath, sHydro, sMeta, 'cf', coefAll(:,ii), 'path', sOutput{ii});
            end
            clear ii

            if ~isempty(gcp('nocreate'))
               pctRunOnAll warning('on', 'all') %Turn off warnings during parfor
            end

    %         %Close dedicated workers used for parallel processing:
    %         poolobj = gcp('nocreate');
    %         delete(poolobj);
        else
            sOutput = CCHF_engine_v4(sPath, sHydro, sMeta);
        end
    elseif regexpbl(sMeta.runType,'default')
        sMeta.mode = 'default';
        sMeta.coef = sMeta.coefAtt(:,[1,4]);
        sOutput = CCHF_engine_v4(sPath, sHydro, sMeta);
    else
        error('CCHF_main:runType', [char(39) sMeta.runType char(39) ...
            ' was selected as the type of model run, but this option is not recognized.']);
    end
    tEnd = now;
    disp(['It took ' num2str(round2((tEnd-tStrt)*24*60,1)) ' minutes for the current model run.']);


    %Save output structures to file:
    if ~exist(sPath{1}.outputMain, 'dir')
       mkdir(sPath{1}.outputMain) 
    end
    save(fullfile(sPath{1}.outputMain, 'model_run_structs.mat'), 'sPath', 'sMeta', 'sHydro', 'sOpt', 'sObs', 'sOutput', '-v7.3');
else %Load output structures saved during previous run
    uiwait(msgbox(sprintf(['You have selected to load output from previous model run(s).' ...
        '\n You will be prompted to select 1 Matlab arrays (*.mat) ' ...
        'corresponding to these run(s).\n']), '(Click OK to Proceed)','modal'));
    
    [fileRun, foldRun] = uigetfile({'*.mat'}, 'Select the Matlab array containing output from the previous run', pwd);
    load(fullfile(foldRun, fileRun));
    nSites = numel(sMeta.region(:));
end


%Create function outputs
if nargout > 0
   varargout{1} = sOutput; 
   if nargout > 1
       varargout{2} = sObs; 
   end
end


%Plot model performance
if ~regexpbl(sMeta.runType,'sim')
    if isfield(sMeta, 'wrtGridEval') && sMeta.wrtGridEval == 1
        wrtGrid = sMeta.wrtTyp;
    else
        wrtGrid = blanks(0);
    end
    warning('off','all'); %Turn off warning that some time-series not being used.
    %Define additional stats to calculate on data
    statsExtra = {'KGEr', 'KGEs', 'KGEb', 'NSE', 'MAPE', 'MAE'};
    if strcmpi(sOpt.fitTest,'kge')
        cellStats = [sOpt.fitTest, statsExtra];
    else
        cellStats = [sOpt.fitTest, 'kge', statsExtra];
    end
    
    if numel(sOutput) == nSites %&& isstruct(sOutput{1})
        for mm = 1 : nSites
            %Create output directory for assessment plots:
            dirModObs = fullfile(sPath{mm}.output, 'mod_v_obs_compare');
            if ~exist(dirModObs, 'dir')
                mkdir(dirModObs);
            end
            
            %Display all stats and write to file
            report_stats_v2(sObs{mm}, sOutput{mm}, cellStats, sPath{mm}.output, sMeta);

            %Create plots: modeled versus observed
            mod_v_obs_v2(sObs{mm}, sOutput{mm}, sOpt.fitTest, 'plot', dirModObs, 'scatter','grid','combineType', 'lon', sHydro{mm}.(varLon), 'lat', sHydro{mm}.(varLat), 'wrtGridTyp', wrtGrid);
            mod_v_obs_v2(sObs{mm}, sOutput{mm}, sOpt.fitTest, 'plot', dirModObs, 'combineType', 'wrtTyp', sMeta.wrtTyp);
        end
        clear mm
    else %This is used when validation is run to assess equifinality
        %This whole conditonal section needs to be updated...
        %Check if output is in memory or saved as path:
        if isstruct(sOutput{1})
            typOut = 'struct';
        elseif ischar(sOutput{1}) || (numel(sOutput{1}(:)) == numel(sMeta.region(:)) && ischar(sOutput{1}{1}))
            typOut = 'path';
        else
            error('CCHF_main:outputTypeUnknown',['The output array has type ' ...
                class(sOutput{1}) ', which has not been programmed for.']);
        end
        
        %Open "Matlab pool" for parallel validation:
        if isempty(gcp('nocreate'))
            localCluster = parcluster('local'); %Check number of possible "workers"
            if isfield(sOpt,'maxWorker') && localCluster.NumWorkers > sOpt.maxWorker
                maxWorkersUse = sOpt.maxWorker;
            else
                maxWorkersUse = localCluster.NumWorkers;
            end
            parpool(maxWorkersUse); %Dedicate all available cores to parallel calibration
        end
        
        %Initialize:
        scoreTemp = cell(numel(sOutput(:)), 1);
        obsTypTemp = scoreTemp;
        
        %Loop over outer-set (either number of sites or number of parameter
        %sets)
        warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
        
        if strcmpi(typOut, 'struct')
            nIter = numel(sOutput{1}(:));
            parfor mm = 1 : numel(sOutput(:)) 
                scoreSameTemp = cell(nIter,1);
                typSameTemp = cell(nIter,1);
                for nn = 1 : nIter
                    [scoreSameTemp{nn}, typSameTemp{nn}] = report_stats_v2(sObs{nn}, sOutput{mm}{nn}, cellStats, '', sMeta, 'no_disp', 'no_write');
                end
                            
                scoreTemp{mm} = scoreSameTemp;
                obsTypTemp{mm} = typSameTemp;
            end
        elseif strcmpi(typOut, 'path')
            parfor mm = 1 : numel(sOutput(:)) 
                scoreSameTemp = cell(nSites,1);
                typSameTemp = cell(nSites,1);
                for nn = 1 : nSites
                    [foldLd, fileLd, ~]= fileparts(sOutput{mm}{nn});
                    temp = load(fullfile(foldLd, [fileLd '_' sMeta.region{nn} '.mat']));
                    nmsTemp = fieldnames(temp);
                    if numel(nmsTemp) == 1
                        [scoreSameTemp{nn}, typSameTemp{nn}] = report_stats_v2(sObs{nn}, temp.(nmsTemp{1}), cellStats, '', sMeta, 'no_disp', 'no_write');
                    else
                        error('CCHFMain:multOutFlds','The output array has an unexpected number of fields.')
                    end
                end
                
                scoreTemp{mm} = scoreSameTemp;
                obsTypTemp{mm} = typSameTemp;
            end
        else
            error('CCHFImplement:unknownOutputType', ['The output type ' class(typOut) ' has not been programmed for.']);
        end
        clear mm
        warning('on', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
          
        if ~isempty(gcp('nocreate'))
           pctRunOnAll warning('on', 'all') %Turn on warnings
        end

        %Close dedicated workers used for parallel processing:
        poolobj = gcp('nocreate');
        delete(poolobj);
        warning('on','all'); %Turn off warning that some time-series not being used.
        
        
        %Initialize output
        scoreOut = cell(numel(sMeta.region(:)), 1);
            
        %Reformat statistics so that each region/obs type is in seperate cell
        %array
        nReg = numel(sMeta.region(:));
        nRuns = numel(scoreTemp(:));
        nStats = numel(cellStats(:));
        %rows are model runs and columns are metrics
        for nn = 1 : nReg %loop over region
            nObsTyp = numel(obsTypTemp{1}{nn}{1}(:));
                        
            scoreOut{nn} = cell(nObsTyp,1);
            [scoreOut{nn}{:}] = deal(nan(nRuns, nStats));
            
            for ll = 1 : nObsTyp %loop over observation type
                for mm = 1 : nRuns %loop over model run
                    for yy = 1 : nStats %loop over metric 
                        scoreOut{nn}{ll}(mm,yy) = scoreTemp{mm}{nn}{yy}(ll);
                    end
                    clear yy
                end
                clear mm
            end
            clear ll
        end
        clear nn

        %Write all validation results to csv files (one for each region - obs type):
        headerFmt = repmat('%s,',1,nStats-1);
        numFmt = repmat('%f,',1,nStats-1);

        hdrStats = cellStats(:)';

        pathStats = cell(nReg, 1);
        statFileRt = cell(nReg, 1);
        for nn = 1 : nReg
            %Create output directory for assessment plots:
            dirMult = fullfile(sPath{nn}.output, ['mod_stats_' num2str(nRuns) '_best_cal_runs']);
            if ~exist(dirMult, 'dir')
                mkdir(dirMult);
            end
            
            nObsTyp = numel(obsTypTemp{1}{nn}{1}(:));
        
            pathStats{nn} = cell(nObsTyp, 1);
            statFileRt{nn} = cell(nObsTyp, 1);
            
            for yy = 1 : nObsTyp
                statFileRt{nn}{yy} = [sMeta.region{nn} '_' strrep(obsTypTemp{1}{nn}{1}{yy},' ','_')];
                pathStats{nn}{yy} = fullfile(dirMult, [statFileRt{nn}{yy} '_performance.csv']);
                if exist(pathStats{nn}{yy}, 'file')
                    delete(pathStats{nn}{yy});
                end
                
                fStats = fopen(pathStats{nn}{yy},'w+');
                
                fprintf(fStats, [headerFmt,'%s\n'], hdrStats{:});
                
                for mm = 1 : nRuns
                    fprintf(fStats, [numFmt,'%f\n'], scoreOut{nn}{yy}(mm,:));
                end
                clear mm
                
                fclose(fStats);
            end
            clear yy
        end
        clear nn


        %Calculate validation statistics and write to file:
        aggStats = {'best_cal_run','mean_val', 'high_val', 'low_val', 'median_val', 'mode_val', 'SD_val'};
        nStats = numel(cellStats);
        headerFmt = repmat('%s,',1,nStats);
        numFmt = repmat('%f,',1,nStats-1);

        hdrStats = [blanks(1), cellStats(:)'];

        for nn = 1 : numel(pathStats(:)) %loop over region
            for yy = 1 : numel(pathStats{nn}(:)) %Loop over observation type
                [foldCurr, ~, ~] = fileparts(pathStats{nn}{yy});
                pathAggCurr = fullfile(foldCurr, [statFileRt{nn}{yy} '_agg_stats.csv']);
                if exist(pathAggCurr, 'file')
                    delete(pathAggCurr);
                end
                fStats = fopen(pathAggCurr,'w+');
                fprintf(fStats, [headerFmt,'%s\n'], hdrStats{:});
                
                for oo = 1 : numel(aggStats)
                    fprintf(fStats, '%s,', aggStats{oo});
                    
                    switch aggStats{oo}
                        case 'best_cal_run'
                            statCurr = scoreOut{nn}{yy}(1,:);
                        case 'mean_val'
                            statCurr = mean(scoreOut{nn}{yy}, 1);
                        case 'high_val'
                            statCurr = max(scoreOut{nn}{yy}, [], 1);
                        case 'low_val'
                            statCurr = min(scoreOut{nn}{yy}, [], 1);
                        case 'median_val'
                            statCurr = median(scoreOut{nn}{yy}, 1);    
                        case 'mode_val'
                            statCurr = mode(scoreOut{nn}{yy}, 1);
                        case 'SD_val'
                            statCurr = std(scoreOut{nn}{yy}, 1);
                        otherwise
                            display(['The current loop is being skipped '...
                                'because ' aggStats{oo} ' is not a recognized case.']);
                            continue
                    end

                    fprintf(fStats, [numFmt,'%f\n'], statCurr);
                end
                clear oo
                
                fclose(fStats); 
            end
            clear yy
        end
        clear nn
    end
end


%MAKE PLOTS OF MODEL OUTPUT
warning('off','MATLAB:MKDIR:DirectoryExists');
if sMeta.blDispOut == 1
    strDispOut = '';
elseif sMeta.blDispOut == 0
    strDispOut = 'no_disp';
end
if numel(sOutput) == nSites
    for mm = 1 : nSites
        dirModPlots = fullfile(sPath{mm}.output, 'model_output_plots');
            mkdir(dirModPlots);
        pathOutRt = fullfile(dirModPlots,'model');

        %Plot output
        plot_CCHF(sOutput{mm}, {'flow'}, sMeta, [pathOutRt '_flow'], strDispOut);
        % plot_CCHF(sOutput{mm}, {'et','pet'}, sMeta, [pathOutRt '_ET'], strDispOut);
        plot_CCHF(sOutput{mm}, {'hfnet','hfrs', 'hfrl', 'hft','hfcp', 'hfgc', 'hfsnc'}, sMeta, [pathOutRt '_heat_components'],'avg', strDispOut);
        % plot_CCHF(sOutput{mm}, {'prsn','rain','et','mrro'}, sMeta, [pathOutRt '_inVsOut'],'avg', strDispOut);  
        % plot_CCHF(sOutput{mm}, {'prsn','swe','sndwe','rain','mrro','icdwe'}, sMeta, [pathOutRt '_inVsOut'],'avg', strDispOut);
%         plot_CCHF(sOutput{mm}, {'rain','snlr','iclr'}, sMeta, [pathOutRt '_inVsOut'],'avg', strDispOut);
        plot_CCHF(sOutput{mm}, {'snlr','rnrf','iclr'}, sMeta, [pathOutRt '_water_release'],'avg', strDispOut);
        % plot_CCHF(sOutput{mm}, {'sndwe','icdwe'}, sMeta, [pathOutRt '_cryoChange'], strDispOut);
        plot_CCHF(sOutput{mm}, {'snw','snlw','sndwe'}, sMeta, [pathOutRt '_snw'], strDispOut);
        plot_CCHF(sOutput{mm}, {'snw','sncc'}, sMeta, [pathOutRt '_snow_coldcontent'], strDispOut);
        % plot_CCHF(sOutput{mm}, {'icdwe'}, sMeta, [pathOutRt '_dIce'], strDispOut);
        % plot_CCHF(sOutput{mm}, {'snw'}, sMeta, [pathOutRt '_snowpack'], strDispOut);
        % plot_CCHF(sOutput{mm}, {'sww','lwsnl','snlh','snlr','iclr','lhpme','rnrf'}, sMeta, [pathOutRt '_snowpack'],'avg', strDispOut);
        plot_CCHF(sOutput{mm}, {'lwsnl','snlh','snlr','iclr','lhpme','sndwe'}, sMeta, [pathOutRt '_snowpack'],'avg', strDispOut);
        plot_CCHF(sOutput{mm}, {'tas','tsis','tsn'}, sMeta, [pathOutRt '_temperature'], strDispOut);
        plot_CCHF(sOutput{mm}, {'prsn','rain'}, sMeta, [pathOutRt '_rain'], strDispOut);
        % plot_CCHF(sOutput{mm}, {'albedoS', 'albedoI'}, sMeta, [pathOutRt '_albedo'], strDispOut);
        % plot_CCHF(sOutput{mm}, 'flow', sMeta, [pathOutRt '_flowrate'], strDispOut);
    end
    clear mm
end
warning('on','MATLAB:MKDIR:DirectoryExists');


%Display message that processing complete and turn off diary:
if regexpbl(sMeta.runType,'calib')
    runTypeDisp = 'calibration';
elseif regexpbl(sMeta.runType,'valid')
    runTypeDisp = 'validation';
elseif regexpbl(sMeta.runType,'default')
    runTypeDisp = 'default parameter set';
elseif regexpbl(sMeta.runType,'sim')
    runTypeDisp = 'simulation';
else
    runTypeDisp = 'unknown';
end
regDisp = '';
for ii = 1 : nSites
   regDisp = [regDisp, sMeta.region{ii}];
   if ii ~= nSites
      regDisp = [regDisp '/']; 
   end
end
clear ii

disp(['The ' runTypeDisp ' run(s) for ' regDisp ' have finished.']);
for ii = 1 : nSites
    disp(['Results have been written to ' char(39) sPath{ii}.output ...
        char(39)]);
end
diary off   %Stop recording command history.
