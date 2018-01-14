% Copyright 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the Conceptual Cryosphere Hydrology Modeling Package 
%(referred to as the 'CCHF Package').
% 
% The CCHF Package is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% 
% The CCHF Package is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the CCHF Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%The first section, 'USER INPUTS' contains all of the options for the user
%to edit.  The second section 'IMPLEMENTATION CODE' should not be altered.
%Note: All strings are case insensitive.

clear all   %Clears variables in MATLAB, ensuring no conflicts.
clc         %Clears all existing text in the command window.
close all   %Closes all figure windows. 





%%USER INPUTS:
%Select the time period to downscale:
runType = 'calibrate';	%Either 'default' (guess a parameter set), 
                            %'calibrate' (for optimizing parameters), 
                            %'calibrate_resume' (use this if a calibration 
                                %routine was interuptted)
                            %'validate' (apply optimized parameter set and 
                                %compare to observation data)
                            %'simulate' (apply optimized parameter set but
                                %do not compare to observations)
                            %FOR 'VALIDATE' OR 'SIMULATE': can append 
                                %'_[Number]' to perform equifinality test

%Select the name for the region (used for naming output files):
region = 'wolverine'; %{'Hunza','Karora'}; %{'Naltar','Karora'}; %{'Donyian','Hunza'}; 

startDate = [2000, 9];  %Vector with format '[year, month]' specifying date 
                        %to begin model run (run will include this date).

endDate = [2003, 8];   %Vector with format '[year, month]' specifying date 
                        %to end model run (run will include this date).                          



     
%Number of gage files to load
%On the first run, must select multiple gage files (if desired) and enter
%metadata. These will be compiled such that during subsquent runs the user
%only locates teh compiled file. The naming scheme for the compiled file is
%'CCHF_gage-*.txt'
nGage = 1;  

                    
%REQUIRED MODULE CHOICES (ONE FROM EACH CATEGORY):
%'heat': 'STI' (simple temperature index), 'Pelli' (Pellicciotti version of ETI), 'Hock', 'LST', or 
    %'SETI'
%'mass': 'step', 'cc', 'enbal', 'Liston'
%'icmlt': 'ratio' (ratio of heat on snow and equivalent ice surface 
    %calculated, then any energy remaining after snowmelt is used for ice 
    %melt), 'time' (determines time left in each time step for ice melt), 
    %'Liston' (glacier melt only occurs during time steps in which there is 
    %no snow and now snow melt has occured), 'Neumann' (implements a
    %constant basal temperature Neumann boundary condition), or 'gradient'
    %(implements ablation gradient method in which ice melt is only 
    %calculated once per year.
%'trans': 'DeWalle' or 'Coops' or 'dem_exp_decay'
%'albedo' (not used if 'cryo-simple'): 'Brock' or 'Pelli'
%'toa': 'DeWalle' (works for daily or monthly) or 'Liston' (works for
    %daily)
%'pet': 'Hammon', 'Hargreaves', 'Makkink'
%'runoff': 'bucket' or 'direct'
%'time' (refers to travel time between cells): 'Johnstone' or 'Liston'
%'flow': 'Muskingum', 'lumped', or 'Liston'
%'snlq': 'percent' (snow liquid water holding capacity equal to fitting
    %parameter representing percentage), 'zero' (all liquid drains
    %immediately), 'density' (there is a snow density threshold, used to
    %determine when drainage occurs
%'density': 'COSIMA' or 'Liston'
%'firn': 'simple' (snow becomes ice once threshold reached) or 'static' (no firn compaction)
%'glaciermove': 'Shea' (model will allow glacier sliding and 
    %calculate changes in ice thickness) or 'static' (glaciers held constant)
moduleSet = { ...
    'heat',  'STI'; ... %surface heat flux module representation
    'mass', 'cc'; ... %snow energy and mass module representation
    'icmlt', 'ratio'; ... %ice melt module
    'runoff', 'bucket'; ... %runoff module representation
    'toa', 'DeWalle'; ... %top-of-atmosphere radiation module representation
    'trans', 'dem_exp_decay'; ... %atmospheric transmissivity module representation
    'albedo', 'Pelli'; ... %albedo module representation
    'glacier0', 'Shea'; ... %glacier water equivalent at time 0
    'pet', 'Hammon'; ... %potential evapotranspiration module representation
    'timelag', 'Johnstone'; ... %timelag (of flow through each grid cell) module representation
    'flow', 'lumped'; ... %flow routing module
    'density', 'Liston'; ... %densification of snow (does not participate in firn to ice processes)
    'snlq', 'percent'; ... %snow holding capacity
    'sca', 'max'; ...
    'sublimate', 'Lutz'; ...
    'firn', 'static'; ... %firn compaction
	'glaciermove', 'static'; ... %glacier movement
    };

%Boolean value to toggle on all glacier processes. If off, glaciers can
%still be included, but they will be static (infinite sources of melt water)
blDynamicGlaciers = 0; %0 = static glaciers, 1 = dynamic (use firn/glacier processes selected in module set) 

%These settings only used in calibration (can to left alone for validation)
optType = 'hybrid';  %optimization method: 
                        %'GA_binary' (Genetic Algorithm using binary chromosomes), 
                        %'GA_real' (genetic algorithm using real numbers), 
                        %'monte_carlo', 
                        %'uniform_sampling', 
                        %'PSO' (Particle Swarm optimization), or 
                        %'hybrid' (combination of PSO, Monte Carlo, and linear sensitivity)
fitType = 'kge_Parajka_mbe'; %fitness score options: 
                    %'KGE' (Kling-Gupta Efficiency), 
                    %'NSE' (Nash-Sutcliffe Efficiency), 
                    %'MAE' (mean absolute error), 
                    %'MAPE' (mean absolute percent error),
                    %'WMAPE' (weighted mean absolute percent error), or
                    %others...
nRuns = 8000;   %Maximum number of runs

%Option to determine how performance at multiple sites is ranked during
%calibration:
methodSiteCombine = 'linear'; %Options: 
                                %'linear' (compute mean of performance at multiple sites)
                                %'square' (compute square root of sum of squares of performance at multiple sites)


%Multiple calibration stages: (if commented out, default is to calibrate
%all paremters simultaneously)
    %Categories that can be calibrated seperately: 
        %'input' (precip and temp), 
        %'cryo' (snow and ice, including heat processes), 
        %'land' (ET, soil), 
        %'routing' (velocity/lag time)
% nStage = 1; %%Number of calibration stages (Use two stages)
% category = {'input','cryo', 'land','routing'};
nStage = 2; %%Number of calibration stages (Use two stages)
category = {'input','cryo'; ...
    'land','routing'}; %Each row corresponds parameter types to optimize in that stage

useFdr = 1; %Set to one to load and use ArcGIS flow direction grid. If 0, 
            %uses algorithm to computer flow direction, which works well on 
            %large watersheds but not well on small watersheds. 
            
%Select the time-series elements to produce:
monthsSpin = '24 months';	%String defining number of months to 
                            %run model for prior to start date
                        %'startDate'.
                        
timeStep = 'daily';     %String specifying time resolution.  Can be 
                        %'daily', 'monthly', 'hourly'.
                        
dataRes = 'daily';    %String specifying resolution of input time-series 
                        %data 
     
maxWorkersCal = 8; %Specifcy maximum number of CPUs to use in parallell optimization.
maxWorkersVal = 4;
                        
%For genetic algorithm only:
perUse = []; %If empty, no roulette selection.
nElite = 2;
rateCross = 0.8;
    nCross = 1;
    siteMutate = 1;
%End of genetic algorithm parameters
         
%POINTS WHERE OUTPUT DESIRED:
%Eklutna:  
%testPts = {[-148.971, 61.296]; [-148.955, 61.147]} ;  
%Wolverine:
% testPts = {[-148.921, 60.4197]; [-148.915, 60.3769]; [-148.907, 60.4042]};
% reportPt = {'avg'; testPts{:}};
% FallsCreek:
% outlet = [-122.349343, 44.379419];

reportPt = {'avg'};

output = {...
    'flow', reportPt; ... %{'all'; reportPt{:}}; ... %surface flowrate through cell
    'mrro', reportPt; ... %surface runoff from cell
    'pr', reportPt; ... %precipitation
    'rain', reportPt; ... %rain
    'prsn', reportPt; ... %snowfall
    'tas', reportPt; ... %2-m air temperature
    'pet', reportPt; ... %Potential evapotranspiration
    'et', reportPt; ... %Actual evapotranspiration
    'snw', reportPt; ... %solid snow (Water equivalent)
    'casi', reportPt; ... %snow/ice covered fraction
    'scx', reportPt; ... %snow covered fraction
    'icx', reportPt; ... %ice covered fraction
    'tsn', reportPt; ... %Snow internal temperature
    'tsis', reportPt; ... %Snow/ice surface temperature
    'sncc', reportPt; ... %Snow cold-content
    'snsb', reportPt; ... %Snow sublimation (we
    'snalb', reportPt; ... %snow albedo
    'sm', reportPt; ... %Soil moisture
    'lhpme', reportPt; ... %Latent heat expressed as potential melt equivalent (meters)
    'hfnet', reportPt; ... %Net heat flux
    'hfrs', reportPt; ... %Shortwave radiation
    'hft', reportPt; ... %degree-index heat flux
    'hfrl', reportPt; ... %Longwave radiation heat flux
    'hfcs', reportPt; ... %sensible convective heat flux
    'hfcp', reportPt; ... %latent precipitation heat flux
    'hfsnc', reportPt; ... %conduction between snow surface and internal snowpack
    'hfgc', reportPt; ... %ground/glacier conductive heat flux
    'icdwe', reportPt; ... %change in ice (water equivalent)
    'sndwe', reportPt; ... %change in snow (water equivalent)
    'rnrf', reportPt; ... %Rain that is added to runoff immediately (i.e. no snow beneath)
    'lwsnl', reportPt; ... %Liquid water content of snow
    'snlh', reportPt; ... %Liquid water holding capacity of snow
    'snlr', reportPt; ... %Liquid water released from snowpack
    'iclr', reportPt; ... %Liquid water released from ice
    'sndt', reportPt; ... %Change in snow temperature
    };
% output = {...
%     'flow', {'all'; reportPt{:}}; ...
%     'mrro', {'all'; reportPt{:}}; ...
%     'pr', reportPt; ... %precipitation
%     'rain', {'all'; reportPt{:}}; ...
%     'prsn', reportPt; ... %snowfall
%     'tas', reportPt; ... %2-m air temperature
%     'pet', reportPt; ... %Potential evapotranspiration
%     'et', reportPt; ... %Actual evapotranspiration
%     'snw', {'all'; reportPt{:}}; ... %Water equivalent of snowpack
%     'casi', reportPt; ... %Water equivalent of snowpack
%     'tsn', reportPt; ... %Snow internal temperature
%     'tss', reportPt; ... %Snow surface temperature
%     'sncc', reportPt; ... %Snow cold-content
%     'albedo', reportPt; ...
%     'sm', reportPt; ... %Soil moisture
%     'lhpme', reportPt; ... %Latent heat melt equivalent of snow.
%     'hfnet', reportPt; ... %Net heat flux
%     'hfrs', reportPt; ... %Shortwave radiation
%     'heatT', reportPt; ... %degree-index heat flux
%     'hfrl', reportPt; ... %Longwave radiation heat flux
%     'hft', reportPt; ... %sensible convective heat flux
%     'hfcp', reportPt; ... %latent precipitation heat flux
%     'hfgc', reportPt; ... %ground/glacier conductive heat flux
%     'icdwe', {'all'; reportPt{:}}; ... %change in ice (water equivalent)
%     'sndwe', {'all'; reportPt{:}}; ... %change in snow (water equivalent)
%     'rnrf', {'all'; reportPt{:}}; ... %Rain that is added to runoff immediately (i.e. no snow beneath)
%     'lwsnl', reportPt; ... %Liquid water content of snow
%     'snlh', reportPt; ... %Liquid water holding capacity of snow
%     'snlr', {'all'; reportPt{:}}; ... %Liquid water released from snowpack
%     'iclr', {'all'; reportPt{:}}; ... %Liquid water released from ice
%     'sndt', reportPt; ... %Change in snow temperature
%     };
            
%Boolean value to display time-series plots of modelled output:
blDispOutput = 1;

%Include glaciers? %Select type of DEM grid for ice. 
iceGrid = 'same'; %Choices are:
                    %'none' = glaciers not used.
                    %'same' = same grid as main DEM.
                    %'fine' = finer grid than main DEM (requires two grids:
                            %one has boolean ice value and corresponding 
                            %high-res DEM)

useDebris = 0;  %Set to 1 to load debris cover array. If 0, all ice presumed 
                %to be clean ice.                         
                            
%Clip to non-rectangular region of interest using reference file:
blClip = 0;     %(1 = yes or 0 = no).  If 1, user prompted to 
                %locate reference file in ESRI ACII format.  All downscaled 
                %output will only be produced for cells in reference file 
                %where value is not NaN.
        

%PRINTING OPTIONS: 
%Determine output file type:
writeMod = 1;       %1 = Gridded files will be written for any outut parameters recorded at grid cells where 'all' observations are recorded (see 'output' and 'reportPts').
                    %0 = No gridded files will be written for modelled output.
writeEval = 1;      %1 = Gridded files of observations and model output will be written during model evaluation.
                    %0 = No gridded files will be written for modelled output.
writeType = 'ascii';    %'ascii' = individual ascii files (ESRI format)
                        %'netCDF' = A single netCDF file (Analogous to CMIP5 format)

%Defines abbreviated calibration period (if ftiness requirement not met,
%cuts model run for parameter set short. 
%Evaluate model two years after start of year to determine worthiness of
%parameter set:
dateEval = nan(1,3);
% dateEval = startDate + [2, zeros(1,numel(startDate)-1)];

blPrevRun = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%IMPLEMENTATION CODE (DO NOT EDIT):    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(region)
    nSites = numel(region(:));
elseif ischar(region)
    nSites = 1;
    region = {region};
else
    error('CCHF_main:nmRegions',['The model dies not recognize the '...
        'format of the region input variable.']);
end

% if nSites ~= numel(regions(:))
%    error('CCHF_main:diffNumRegions',['The model is set to evaluate for ']) 
% end


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
            

%'blPrevRun' provides option to load output from previous run (not normal
%condition)
if blPrevRun == 0
    % global sHydro{ii}
    pathGage = cell(nSites, 1);
    sPath = cell(nSites, 1);
    sHydro = cell(nSites, 1);
    sObs = cell(nSites, 1);
    sMeta = struct;
        sMeta.('runType')  = runType;
        sMeta.('dateStart')= startDate; 
        sMeta.('dateEnd')  = endDate;
        sMeta.('spinup')   = monthsSpin;
        sMeta.('dt')       = timeStep;
        sMeta.('dataTsRes')= dataRes;
        sMeta.('region')   = region;
        sMeta.('module')   = moduleSet;
        sMeta.('mode')     = blanks(0);
        sMeta.('coef')     = cell(1,2);
        sMeta.('wrtTyp')   = writeType;
        sMeta.('wrtGridOut')  = writeMod;
        sMeta.('wrtGridEval') = writeEval;
        sMeta.('output')   = cell(nSites, 1);
        sMeta.('iceGrid')  = iceGrid;
        sMeta.('rtDir')    = cell(nSites, 1);
        sMeta.('glacierDynamics') = blDynamicGlaciers;
        sMeta.varLd        = {'pr','tas','tasmin','tasmax'};


    sOpt = struct;
        sOpt.('fitTest') = fitType;
        sOpt.('methodSiteCombine') = methodSiteCombine;

    if regexpbl(sMeta.runType,'calibrate')
        sOpt.maxWorker = maxWorkersCal;
        sOpt.('type')    = optType; 

        sOpt.nbest = 100; %Record best _n_ performing parameter sets

        %Determine number of generations based on model runs and a
        %fixed population of 30
        nPop = 30;
%         nGen = 5;
        nGen = round(nRuns / nPop);

        if nGen < 1
            nGen = 1;
        end
        if nGen < 20
           warning('CCHF_main:nGenCalibrate',['Currently the number of '...
               'generations is set to ' num2str(nGen) ...
               '.  The recommended minimum is 20.']); 
        end

        sOpt.('att') = {...
            'generations', nGen; ...
            'population', nPop};

        %Set parameters for number of calibration stages:
        sOpt.nStage = nStage;
        if nStage > 1
            sOpt.stageVar = category;
        end

        if regexpbl(sOpt.type, 'uniform') && nStage > 1
            disp(['The number of calibration stages is being set to 1 '...
                'because the otpimization type is ' char(39) sOpt.type ...
                char(39) '.'])
        end

        %Set number of generations until stagnation:
        if regexpbl(sOpt.type,{'genetic','GA'})
            sOpt.stagnate = 20;
        elseif regexpbl(sOpt.type,{'PSO','hybrid'})
            sOpt.stagnate = 10;
        else
            sOpt.stagnate = 15;
        end

        sOpt.dateEval = dateEval;
        %If genetic algorithm used, assign specific fields:
        if regexpbl(sOpt.type,{'genetic','GA'})
            nMemCross = ceil(rateCross*(nPop - nElite));
            if mod(nMemCross, 2) %Ensure even number of parameter sets crossed
                nMemCross = nMemCross + 1;
            end
            nMemMut = nPop - nElite - nMemCross;
            sOpt.att(end+1:end+6,:) = {...
                'fitness requirement', perUse; ...
                'crossover pop', nMemCross; ...
                'crossover sites', nCross; ...
                'mutation pop', nMemMut; ...
                'mutation sites', siteMutate; ...
                'elite pop', nElite};
        end
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
        if useFdr == 1
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
            warning('CCHF_main:fineDemNoGlacier',['The glacier grid resolution '...
                'is being set to be the same as the main grid because no '...
                'glacier dynamics module has been selected.']);
            sMeta.iceGrid = 'same';
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
                sMeta.region{ii} iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
            [fileGlac, foldGlac] = uigetfile({'*.*'},['Select the glacier presence geo-referenced grid or shapefile for ' ...
                sMeta.region{ii} iceAuxDisp '.\n'], startPath);
            sPath{ii}.ice = fullfile(foldGlac, fileGlac);
            disp([char(39) sPath{ii}.ice char(39) ' has been chosen as the glacier presence grid.']);

            if isempty(fileGlac) || isequal(fileGlac, 0)
               error('CCHF_main:noGlacier',['No glacier presence grid '...
                   'has been selected, even though the option was chosen.' ...
                   ' Therefore, the program is aborting.']); 
            end

            uiwait(msgbox(sprintf(['Select the glacier surface DEM for ' ...
                sMeta.region{ii} iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
            [fileGlacDem, foldGlacDem] = uigetfile({'*.*'},['Select the glacier surface DEM for ' ...
                sMeta.region{ii} iceAuxDisp '.\n'], startPath);
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
        if useDebris && ~regexpbl(sMeta.iceGrid, 'none')
            %Flow direction selection and display:
            uiwait(msgbox(sprintf(['Select the debris cover grid for ' ...
                sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
            [fileDeb, foldDeb] = uigetfile({'*.asc';'*.txt'},['Select the debris cover grid for ' ...
                sMeta.region{ii}], startPath);
            sPath{ii}.debris = fullfile(foldDeb, fileDeb);
            disp([char(39) sPath{ii}.debris char(39) ' has been chosen as the debris cover grid.']);

            if isempty(fileDeb) || isequal(fileDeb, 0)
               error('CCHF_main:noDebris',['No glacier surface debris grid '...
                   'has been selected, even though the option was chosen.' ...
                   ' Therefore, the program is aborting.']); 
            end
        end


        %If a binary grid is being used to clip the region:
        if blClip == 1
            uiwait(msgbox(sprintf(['Select the binary region grid that will be used to clip the ' ...
                sMeta.region{ii} ' region.\n']), '(Click OK to Proceed)','modal'));
            [fileRegClip, foldRegClip] = uigetfile({'*.asc';'*.txt'},['Select the binary region clip grid for ' ...
                sMeta.region{ii}], startPath);
            sPath{ii}.regionClip = fullfile(foldRegClip, fileRegClip);
            disp([char(39) sPath{ii}.regionClip char(39) ' has been chosen as the binary region clip grid.']);

            if isempty(fileRegClip) || isequal(fileRegClip, 0)
               error('CCHF_main:noClip',['No region clipping file '...
                   'has been selected, even though the option was chosen.' ...
                   ' Therefore, the program is aborting.']); 
            end
        end


        %Identify climate variables to load:
        [sMeta.varLd, sMeta.varLdDisp] = clm_var_load(sMeta.module);
        %Load climate variables:
        sPath{ii} = clm_path_ui(sPath{ii}, sMeta, sMeta.region{ii});


        %Load DEM and calculate watershed geometry fields:
        sDem = read_geodata(sPath{ii}.dem, sMeta, 'none');
        sHydro{ii}.dem = sDem.data;
        sHydro{ii}.lat = sDem.lat;
        sHydro{ii}.lon = sDem.lon;
        if isfield(sPath{ii}, 'fdr')
            sFdr = read_geodata(sPath{ii}.fdr , sMeta, 'none');
            sHydro{ii}.fdrESRI = sFdr.data;
            clear('sFdr');
        else
            warning('cchfMain:noFDR',['CCHF currently does not have the '...
                'capacity to calculate a flow direction grid. Instead, '...
                'this must be generated in a standalone program such as ArcMap.']);
        end
        clear('sDem');

        %Load observation data to use for calibration or validation:
        if regexpbl(runType,{'calibrat','validat','default'})
            %Calibration/validation data required:
            uiwait(msgbox(sprintf(['Select the file(s) with observation data for ' ...
                sMeta.region{ii} '.  If data does not include information such as location, you will '...
                'be required to enter that in the main panel.']), ...
                '(Click OK to Proceed)','modal'));

            pathGage{ii} = cell(nGage,1);
            for jj = 1 : nGage
                [fileGage, foldGage] = uigetfile({'*.*'}, ...
                    ['Select observation data ' num2str(jj) ' of ' ...
                    num2str(nGage) ' for ' sMeta.region{ii}], startPath);
                pathGage{ii}{jj} = fullfile(foldGage, fileGage);

                disp([pathGage{ii}{jj} ' has been selected as observation data file ' num2str(jj) ' of ' num2str(nGage) ' for ' sMeta.region{ii} '.']);
            end
        end

        %Create unique 'main' output directory based upon inputs:
        if ii == 1 
            if isfield(sPath{ii},'resume')
    %             sPath{ii}.output = sPath{ii}.resume;
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
    end %End of loop to select sites 
    clear ii


    %START PROCESSING LOG to record all messages displayed on screen to file in output directory
    [fileDiary] = downscale_diary(sPath{1}.outputMain);
    disp(['All displayed text is being written to ' char(39) fileDiary ...
        char(39) '.']);


    %display message for location of output folder:
    if numel(nSites) > 1
       uiwait(msgbox(sprintf(['Model outputs will be written to a subfolder of ' ...
           sMeta.regions{1} ' because this region is listed first.\n']), ...
           '(Click OK to Proceed)','modal'));
    end


    %GET AND LOAD PARAMETER FILE (IF VALIDATION OR SIMULATION RUN)
    if regexpbl(runType,{'valid','sim'})
        %Load parameters from file:
        uiwait(msgbox(sprintf(['Select the set of parameter coefficients ' ...
            'written during a calibration run of the CCHF model for ' ...
            sMeta.region{1} '.\n']), '(Click OK to Proceed)','modal'));
        [fileCoef, foldCoef] = uigetfile({'*.txt'; '*.csv'}, ['Select the set of parameter coefficients for ' ...
            sMeta.region{1}], startPath);
        for ii = 1 : nSites
            sPath{ii}.coef = fullfile(foldCoef, fileCoef);
        end
        clear ii
        disp([char(39) sPath{1}.coef char(39) ' has been chosen as the set of parameter coefficients.']);
    end


    %LOAD ALL OBSERVATION DATA
    for ii = 1 : nSites
        sMeta.siteCurr = ii;

        if regexpbl(runType,{'calibrat','validat','default'})
            %Remove any duplicate or empty elements:
            pathGage{ii} = unique(pathGage{ii}(~cellfun('isempty',pathGage{ii})));  
            
            for zz = numel(pathGage{ii}(:)) : -1 : 1
                if isnumeric(pathGage{ii}{zz}) || isempty(regexpi(pathGage{ii}{zz}, '[a-z]'))
                    pathGage{ii}(zz) = [];  
                end
            end
            
            %Load observation data:
            sObs{ii} = read_gagedata(pathGage{ii}, sHydro{ii}.lon, sHydro{ii}.lat, ...
                'time',[sMeta.dateStart;sMeta.dateEnd], ...
                'mask',sHydro{ii}.dem);
        else
            sObs{ii} = struct;
        end

        if isempty(sMeta.output{ii})
           sMeta.output{ii} =  output;
        end
        [sMeta.output{ii}, sObs{ii}] = output_all_gage(sMeta.output{ii}, sObs{ii}, sHydro{ii}.lon, sHydro{ii}.lat);
    end
    clear ii


    for ii = 1 : nSites
        sMeta.siteCurr = ii;

        if blClip == 1 && isfield(sPath{ii}, 'regionClip') 
            sRegClip = read_geodata(sPath{ii}.regionClip, sMeta, 'none');
            sHydro{ii}.regionClip = sRegClip.data;
            clear('sRegClip');
        end

        %Calculate fdr, fac, flow calculation order, slope, aspect, distance between neighboring centroids
        sHydro{ii} = watershed(sHydro{ii});

        %Initialize ice grid:
        [sHydro{ii}.sIceInit, sHydro{ii}.icbl] = ice_grid_load_v2(sHydro{ii}, sPath{ii}, sMeta);
    end
    clear ii


    %%Find paremeters needed for current version of model:
    %Do this by running CCHF model in 'parameter' mode
    sMeta.mode = 'parameter';
    sMeta.coefAtt = CCHF_engine_v4(sPath, sHydro, sMeta);


    %IF VALIDATION, LOAD COEFFICIENT SET AND COMPARE TO PARAMETERS NEEDED FOR CURRENT MODEL:
    if regexpbl(runType,{'valid','sim'})
        if regexpbl(sPath{1}.coef,'.txt')
            [cfTemp, hdrParam] = read_CCHF_coef(sPath{1}.coef);

            %Check that loaded coefficients are same as those needed for current
            %model:
            if ~regexpbl(cfTemp(:,1), sMeta.coefAtt(:,1),'and')
                %Find missing parameter:
                cfBl = zeros(numel(sMeta.coefAtt(:,1)),1);
                for ii = 1 : numel(sMeta.coefAtt(:,1))
                    if regexpbl(cfTemp(:,1),sMeta.coefAtt{ii,1})
                        cfBl(ii) = 1;
                    end
                end
                indCf = find(cfBl == 0);
                if ~isempty(indCf)
                    nmQuit = -9999;
                    uiwait(msgbox(sprintf(['The loaded parameter set does not '...
                        'include ' num2str(numel(indCf)) ' parameters needed '...
                        'for the ' char(39) sMeta.strModule char(39) ' model '...
                        'version.' char(10) 'You will be prompted to enter '...
                        'values for these parameters. If you wish to cancel '...
                        'this model run, enter ' num2str(nmQuit) ' for one of the values.']),...
                        '(Click OK to Proceed)','modal'));

                    for ii = 1 : numel(indCf)
                       cfCurr = input(['Enter a value for ' char(39) ...
                           sMeta.coefAtt{indCf(ii),1} char(39) '. The bounds are ' ...
                           num2str(sMeta.coefAtt{indCf(ii),2}) ' thru ' ...
                           num2str(sMeta.coefAtt{indCf(ii),3}) ' and the parameter is '...
                           'used in ' char(39) sMeta.coefAtt{indCf(ii),5} char(39) ':' char(10)]); 
                       cfTemp(end+1,:) = {sMeta.coefAtt{indCf(ii),1}, cfCurr};
                       if cfCurr == nmQuit
                           error('CCHF_main:missingParam',['The exit command '...
                               'was entered during the process of entinering a'...
                               ' missing parameter value.']);
                       end
                    end
                    clear ii
                end
            end

            sMeta.coef = cfTemp;
        elseif regexpbl(sPath{1}.coef,'.csv')
            prmArray = csvread(sPath{1}.coef,1,0);
            fidPrm = fopen(sPath{1}.coef,'r');
            varPrm = textscan(fidPrm,'%s',numel(prmArray(1,:)),'Delimiter',',');  
                varPrm = varPrm{1};
            fclose(fidPrm);

            if numel(varPrm(:)) - 1 == numel(sMeta.coefAtt(:,1))
                varPrm = varPrm(1:end-1)'; %Last column is evaluation metric
                prmArray = prmArray(:,1:end-1);
            elseif numel(varPrm(:)) ~= numel(sMeta.coefAtt(:,1))
                error('CCHF_main:nParam',['The number of parameters loaded '...
                    'is not equal to the number required by this model formulation.']);
            end

            %Check that all parameters present:
            if ~regexpbl(varPrm, sMeta.coefAtt(:,1))
                error('CCHF_main:missingParam',['There is a missing '...
                    'parameter in the current list of parameter sets.']);
            end

            %Check that order of parameters is correct:
            for ii = 1 : numel(varPrm)
                if ~regexpbl(varPrm{ii}, sMeta.coefAtt{ii,1})
                    error('CCHF_main:outoforderParam',[varPrm{ii} ' in the '...
                        'loaded parameter set is out of order relative to '...
                        'the parameter set identified for the current model '...
                        'formulation.']);
                end
            end
            clear ii

            %If loaded parameters contain multiple sets, put each in seperate
            %cell array and calculate correlation between parameters:
            if all(size(prmArray) ~= 1)
                prmCellTemp = cell(numel(prmArray(:,1)),1);
                for mm = 1: numel(prmCellTemp(:))
                    prmCellTemp{mm} = prmArray(mm, :);
                end

                sMeta.coef = prmCellTemp;
            else
                sMeta.coef = prmArray;
            end
        else
            [~, ~, extCf] = fileparts(sPath{1}.coef);
            error('CCHF_main:unknownCoefFormat',['The coeffienct file is in a ' ...
                char(39) extCf char(39) ' format, which has not been programmed for.']);
        end
    elseif ~regexpbl(sMeta.runType,'calibrate') && numel(sMeta.coefAtt(1,:)) >= 6
        %Reduce dimensions of sMeta.coef (loaded during 'parameter' run)
       sMeta.coef = sMeta.coefAtt(:,[1,4,5]); 
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

        sMeta.rtDir{ii} = sPath{ii}.dem;
    end
    clear ii


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

        if iscell(sMeta.coef) && ~ischar(sMeta.coef{1}) && ~all2d( size(sMeta.coef{1}) == 1 )
            nReg = numel(sMeta.region(:));
            pathValOut = cell(nReg, 1);
            for ii = 1 : nReg
                pathValOut{ii} = fullfile(sPath{ii}.output, 'model_output');
                if ~exist(pathValOut{ii}, 'dir')
                    mkdir(pathValOut{ii});
                end
            end
            
            sMeta = rmfield_x(sMeta,'progress');

            coefAll = nan(numel(sMeta.coef(:,1)), numel(sMeta.coef{1}));
            for ii = 1 : numel(coefAll(:,1))
                coefAll(ii,:) = sMeta.coef{ii,:};
            end
            clear ii

            %Open "Matlab pool" for parallel validation:
            if isempty(gcp('nocreate'))
                localCluster = parcluster('local'); %Check number of possible "workers"
                if exist('maxWorkersVal', 'var') && localCluster.NumWorkers > maxWorkersVal
                    maxWorkersCurr = maxWorkersVal;
                else
                    maxWorkersCurr = localCluster.NumWorkers;
                end
                parpool(maxWorkersCurr); %Dedicate all available cores to parallel calibration
            end

            if ~isempty(gcp('nocreate'))
               pctRunOnAll warning('off', 'all') %Turn off warnings during parfor
            end

            sOutput = cell(numel(coefAll(:,1)),1);
            nRuns = numel(coefAll(:,1));
            parfor ii = 1 : numel(coefAll(:,1))
                if mod(ii,10) == 0
                   display(['Currently on validation run ' num2str(ii) ' of ' num2str(nRuns)]); 
                end
                %sMeta.coef = coefAll(ii,:);
%                 sOutTemp = struct;

                sOutput{ii} = fullfile(pathValOut, [valRt num2str(ii) '.mat']);
                
%                 sPath.wrtoutput = sOutput{ii};
                [~] = CCHF_engine_v4(sPath, sHydro, sMeta, coefAll(ii,:)', 'path', sOutput{ii});
%                 sOutput{ii} = CCHF_engine_v4(sPath, sHydro, sMeta, coefAll(ii,:)');
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
    if strcmpi(fitType,'kge')
        cellStats = [fitType, statsExtra];
    else
        cellStats = [fitType, 'kge', statsExtra];
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
            mod_v_obs_v2(sObs{mm}, sOutput{mm}, fitType, 'plot', dirModObs, 'scatter','grid','combineType', 'lon', sHydro{mm}.lon, 'lat', sHydro{mm}.lat, 'wrtGridTyp', wrtGrid);
            mod_v_obs_v2(sObs{mm}, sOutput{mm}, fitType, 'plot', dirModObs,'combineType', 'wrtTyp', sMeta.wrtTyp);
        end
        clear mm
    else %This is used when validation is run to assess equifinality
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
                maxWorkersCal = sOpt.maxWorker;
            else
                maxWorkersCal = localCluster.NumWorkers;
            end
            parpool(maxWorkersCal); %Dedicate all available cores to parallel calibration
        end
        
        %Initialize:
        scoreTemp = cell(numel(sOutput(:)), 1);
        obsTypTemp = scoreTemp;
        
        parfor mm = 1 : numel(sOutput(:))
            if strcmpi(typOut, 'struct')
                sAssess = sOutput{mm};
            elseif strcmpi(typOut, 'path')
                
                sAssess = cell(numel(sMeta.region(:)), 1);
                for nn = 1 : numel(sMeta.region(:))
                    [foldLd, fileLd, ~]= fileparts(sOutput{mm}{nn});
                    temp = load(fullfile(foldLd, [fileLd '_' sMeta.region{nn} '.mat']));
                    nmsTemp = fieldnames(temp);
                    if numel(nmsTemp) == 1
                        sAssess{nn} = temp.(nmsTemp{1});
                    else
                        error('CCHFMain:multOutFlds','The output array has an unexpected number of fields.')
                    end
                end
            end

            scoreSameTemp = cell(numel(sAssess),1);
            typSameTemp = scoreSameTemp;
            
            for nn = 1 : numel(sAssess)
                [scoreSameTemp{nn}, typSameTemp{nn}] = report_stats_v2(sObs{nn}, sAssess{nn}, cellStats, '', sMeta, 'no_disp', 'no_write');
            end

            scoreTemp{mm} = scoreSameTemp;
            obsTypTemp{mm} = typSameTemp;
        end
        clear mm
        
            
        if ~isempty(gcp('nocreate'))
           pctRunOnAll warning('on', 'all') %Turn on warnings during parfor
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
if blDispOutput == 1
    strDispOut = '';
elseif blDispOutput == 0
    strDispOut = 'no_disp';
end
if numel(sOutput) == nSites
    for mm = 1 : nSites
        dirModPlots = fullfile(sPath{mm}.output, 'model_output_plots');
            mkdir(dirModPlots);
        pathOutRt = fullfile(dirModPlots,'model');

        %Plot output
        plot_CCHF(sOutput{mm}, {'flow'}, sMeta, [pathOutRt '_flow'], strDispOut);
        % plot_CCHF(sOutput, {'et','pet'}, sMeta, [pathOutRt '_ET'], strDispOut);
        plot_CCHF(sOutput{mm}, {'hfnet','hfrs', 'hfrl', 'hft','hfcp', 'hfgc', 'hfsnc'}, sMeta, [pathOutRt '_heat_components'],'avg', strDispOut);
        % plot_CCHF(sOutput, {'prsn','rain','et','mrro'}, sMeta, [pathOutRt '_inVsOut'],'avg', strDispOut);  
        % plot_CCHF(sOutput, {'prsn','swe','sndwe','rain','mrro','icdwe'}, sMeta, [pathOutRt '_inVsOut'],'avg', strDispOut);
%         plot_CCHF(sOutput{mm}, {'rain','snlr','iclr'}, sMeta, [pathOutRt '_inVsOut'],'avg', strDispOut);
        plot_CCHF(sOutput{mm}, {'snlr','rnrf','iclr'}, sMeta, [pathOutRt '_water_release'],'avg', strDispOut);
        % plot_CCHF(sOutput, {'sndwe','icdwe'}, sMeta, [pathOutRt '_cryoChange'], strDispOut);
        % plot_CCHF(sOutput, {'swe'}, sMeta, [pathOutRt '_dSnowTmp'], strDispOut);
        plot_CCHF(sOutput{mm}, {'sncc'}, sMeta, [pathOutRt '_snowcoldcontent'], strDispOut);
        % plot_CCHF(sOutput, {'icdwe'}, sMeta, [pathOutRt '_dIce'], strDispOut);
        % plot_CCHF(sOutput, {'snw'}, sMeta, [pathOutRt '_snowpack'], strDispOut);
        % plot_CCHF(sOutput, {'swe','lwsnl','snlh','snlr','iclr','lhpme','rnrf'}, sMeta, [pathOutRt '_snowpack'],'avg', strDispOut);
        plot_CCHF(sOutput{mm}, {'lwsnl','snlh','snlr','iclr','lhpme','sndwe'}, sMeta, [pathOutRt '_snowpack'],'avg', strDispOut);
        plot_CCHF(sOutput{mm}, {'tas','tsis','tsn'}, sMeta, [pathOutRt '_temperature'], strDispOut);
        plot_CCHF(sOutput{mm}, {'prsn','rain'}, sMeta, [pathOutRt '_rain'], strDispOut);
        % plot_CCHF(sOutput, {'albedoS', 'albedoI'}, sMeta, [pathOutRt '_albedo'], strDispOut);
        % plot_CCHF(sOutput, 'flow', sMeta, [pathOutRt '_flowrate'], strDispOut);
    end
    clear mm
end


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
