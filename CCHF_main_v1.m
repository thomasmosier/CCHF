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
%Select the name for the region (used for naming output files):
region = 'wolverine';

%Select the time period to downscale:
runType = 'validate';    %Either 'default' (guess a parameter set), 
                            %'calibrate' (for optimizing parameters), 
                            %'calibrate_resume' (use this if a calibration 
                                %routine was interuptted)
                            %'validate' (apply optimized parameter set and 
                                %compare to observation data)
                            %'simulate' (apply optimized parameter set but
                                %do not compare to observations)
                  
nGage = 1;  %Number of gage files to load (if using 'CCHF_gage-*.txt' compiled gagedate written in 
            %previous run, set nGage to 1); 
                    
%REQUIRED MODULE CHOICES (ONE FROM EACH CATEGORY):
%'heat': 'simple', 'Pelli' (Pellicciotti version of ETI), 'Hock', or 'LST'
%'mass': 'step', 'cc', 'enbal', 'Liston'
%'icmlt': 'ratio' (ratio of heat on snow and equivalent ice surface 
    %calculated, then any energy remaining after snowmelt is used for ice 
    %melt), 'time' (determines time left in each time step for ice melt), 
    %'Liston' (glacier melt only occurs during time steps in which there is 
    %no snow and now snow melt has occured), or 'Neumann' (implements a
    %constant temperature Neumann boundary condition).
%'trans': 'Coops' or 'DeWalle'
%'albedo' (not used if 'cryo-simple'): 'Brock' or 'Pelli'
%'toa': 'DeWalle' (works for daily or monthly) or 'Liston' (works for
    %daily)
%'pet': 'Hammon', 'Hargreaves', 'Makkink'
%'runoff': 'bucket' or 'direct'
%'time' (refers to travel time between cells): 'Johnstone' or 'Liston'
%'flow': 'Muskingum', 'lumped', or 'Liston'
%'snlq': 'percent' (snow liquid water holding capacity equal to fitting
    %parameter representing percentage) or 'zero' (all liquid drains
    %immediately)
moduleSet = { ...
    'heat',  'LST'; ... %surface heat flux module representation
    'mass', 'cc'; ... %snow energy and mass module representation
    'icmlt', 'time'; ... %ice melt module
    'runoff', 'bucket'; ... %runoff module representation
    'toa', 'DeWalle'; ... %top-of-atmosphere radiation module representation
    'trans', 'DeWalle'; ... %atmospheric transmissivity module representation
    'albedo', 'Pelli'; ... %albedo module representation
    'pet', 'Hammon'; ... %potential evapotranspiration module representation
    'timelag', 'Johnstone'; ... %timelag (of flow through each grid cell) module representation
    'flow', 'lumped'; ... %flow routing module
    'snlq', 'percent'; ... %snow holding capacity
    'sca', 'max'; ...
	'glacier', 'static'; ... %glacier dynamic processes
    };



%These settings only used in calibration (can to left alone for validation)
optType = 'hybrid';  %optimization method: 
                        %'GA_binary' (Genetic Algorithm using binary chromosomes), 
                        %'GA_real' (genetic algorithm using real numbers), 
                        %'monte_carlo', 
                        %'uniform_sampling', 
                        %'PSO' (Particle Swarm optimization), or 
                        %'hybrid' (combination of PSO, Monte Carlo, and linear sensitivity)
fitType = 'kge_Parajka'; %fitness score options: 
                    %'KGE' (Kling-Gupta Efficiency), 
                    %'NSE' (Nash-Sutcliffe Efficiency), 
                    %'MAE' (mean absolute error), 
                    %'MAPE' (mean absolute percent error),
                    %'WMAPE' (weighted mean absolute percent error), or
                    %others...
nRuns = 8000;   %Maximum number of runs


%Multiple calibration stages: (if commented out, default is to calibrate
%all paremters simultaneously)
    %Categories that can be calibrated seperately: 
        %'input' (precip and temp), 
        %'cryo' (snow and ice, including heat processes), 
        %'land' (ET, soil), 
        %'routing' (velocity/lag time)
nStage = 2; %%Number of calibration stages (Use two stages)
category = {'input','cryo'; ...
    'land','routing'}; %Each row corresponds parameter types to optimize in that stage

useFdr = 1; %Set to one to load and use ArcGIS flow direction grid. If 0, 
            %uses algorithm to computer flow direction, which works well on 
            %large watersheds but not well on small watersheds. 
            
%Select the time-series elements to produce:
startDate = [2000, 9];  %Vector with format '[year, month]' specifying date 
                        %to begin model run (run will include this date).
                        
endDate = [2003, 8];     %Vector with format '[year, month]' specifying date 
                        %to end model run (run will include this date).

monthsSpin = '12 months';	%String defining number of months to 
                            %run model for prior to start date
                        %'startDate'.
                        
timeStep = 'daily';     %String specifying time resolution.  Can be 
                        %'daily', 'monthly', 'hourly'.
                        
dataRes = 'daily';    %String specifying resolution of input time-series 
                        %data 
     
workers = 8; %Specifcy maximum number of CPUs to use in parallell optimization.                        
                        
%For genetic algorithm only:
perUse = []; %If empty, no roulette selection.
nElite = 2;
rateCross = 0.8;
    nCross = 1;
    siteMutate = 1;
%End of genetic algorithm parameters
         
%POINTS WHERE OUTPUT DESIRED:

%For Wolverine use:
testPts = {[-148.921, 60.4197]; [-148.915, 60.3769]; [-148.907, 60.4042]};
reportPt = vertcat('avg', testPts);

%To only report average conditions over the domain:
%reportPt = {'avg'};
%To report at the above points and output at all grid cells:
%reportPt = vertcat('all', testPts);

%All of the model outputs that will be reported:
output = {...
    'flow', reportPt; ... %surface flowrate through cell
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


%Include glaciers? 
%Select type of DEM grid for ice. *In Version 1, leave value as 'none' or 'same'... 
iceGrid = 'same'; %Choices are:
                    %'none' = glaciers not used.
                    %'same' = same grid as main DEM.
                    %'fine' = finer grid than main DEM (requires two grids:
                            %one has boolean ice value and corresponding 
                            %high-res DEM)

useDebris = 0;  %Set to 1 to load debris cover array. If 0, all ice presumed 
                %to be clean ice. Including debris cover adds scalar 
                %fitting parameter in icemelt functions                         
                            
%Clip to non-rectangular region of interest using reference file:
blClip = 0;     %(1 = yes or 0 = no).  If 1, user prompted to 
                %locate reference file in ESRI ACII format.  All downscaled 
                %output will only be produced for cells in reference file 
                %where value is not NaN.
        

%PRINTING OPTIONS: 
%Determine output file type:
writeType = 'ascii';    %'ascii' = individual ascii files (ESRI format)
                        %'netCDF' = A single netCDF file (Analogous to CMIP5 format)
writeAll = 1;       %1 = NetCDF files will be written for any outut parameters recorded at 'all' grid cells.
                    %0 = No files will be written for any outut parameters recorded at 'all' grid cells.

%Defines abbreviated calibration period (if ftiness requirement not met,
%cuts model run for parameter set short. 
%Evaluate model two years after start of year to determine worthiness of
%parameter set:
dateEval = startDate + [2, zeros(1,numel(startDate)-1)];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%IMPLEMENTATION CODE (DO NOT EDIT):    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
cd(pathScript);
addpath(genpath(pathScript));

%Initialize structure arrays for passing variables to downscaling script:
clear global sHydro sAtm sCryo sLand %Possibly simply initialize instead of clearing
% global sHydro
sPath = struct;
    sPath.('dem')    = blanks(0);
    sPath.('pr')     = blanks(0);
    sPath.('tas')    = blanks(0);
    sPath.('tasmin') = blanks(0);
    sPath.('tasmax') = blanks(0);
    sPath.('regclip')= blanks(0);
    sPath.('output') = blanks(0);
    sPath.('coef')   = blanks(0);
sMeta = struct;
    sMeta.('runType')  = runType;
    sMeta.('dateStart')=  startDate; 
    sMeta.('dateEnd')  = endDate;
    sMeta.('spinup')   = monthsSpin;
    sMeta.('dt')       = timeStep;
    sMeta.('dataTsRes')= dataRes;
    sMeta.('region')   = region;
    sMeta.('module')   = moduleSet;
    sMeta.('mode')     = blanks(0);
    sMeta.('coef')     = cell(1,2);
    sMeta.('wrtTyp')   = writeType;
    sMeta.('wrtAll')   = writeAll;
    sMeta.('output')   = output;
    sMeta.('iceGrid')  = iceGrid;
    sMeta.varLd        = {'pr','tas','tasmin','tasmax'};
sHydro = struct;

    
sOpt = struct;
    sOpt.('fitTest') = fitType;

if regexpbl(sMeta.runType,'calibrate')
    sOpt.maxWorker = workers;
    sOpt.('type')    = optType; 

    sOpt.nbest = 100; %Record best _n_ performing parameter sets
    
    %Determine number of generations based on model runs and a
    %fixed population of 30
    nPop = 30;
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
        sOpt.stagnate = 12;
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


%GUI TO LOCATE INPUT FOLDERS:
startPath = pwd;
%Digital Elevation Model (DEM) selection and display:
uiwait(msgbox(sprintf(['Select the Digital Elevation Model (DEM) for ' ...
    sMeta.region '.\n']), '(Click OK to Proceed)','modal'));
[fileDem, foldDem] = uigetfile({'*.asc';'*.txt'},['Select the Digital Elevation Model '...
    '(DEM) for ' sMeta.region], startPath);
sPath.dem = fullfile(foldDem, fileDem);
disp([char(39) sPath.dem char(39) ' has been chosen as the DEM.']);

%Update search path:
[startPath, ~, ~] = fileparts(sPath.dem);

if isempty(fileDem) || isequal(fileDem, 0)
   error('CCHF_main:noDEM','A DEM has not been selected. Therefore, the program is aborting.'); 
end

%Load manual flow direction grid (ESRI ASCII Format):
if useFdr == 1
    %Flow direction selection and display:
    uiwait(msgbox(sprintf(['Select the flow direction grid (typically created using ArcGIS) for ' ...
        sMeta.region '.\n']), '(Click OK to Proceed)','modal'));
    [fileFdr, foldFdr] = uigetfile({'*.asc';'*.txt'},['Select the flow direction grid for ' ...
        sMeta.region], startPath);
    sPath.fdr = fullfile(foldFdr, fileFdr);
    disp([char(39) sPath.fdr char(39) ' has been chosen as the flow direction grid.']);
    
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
        sMeta.region '.\n']), '(Click OK to Proceed)','modal'));
    [fileGlac, foldGlac] = uigetfile({'*.*'},['Select the binary glacier presence grid for ' ...
        sMeta.region], startPath);
    sPath.ice = fullfile(foldGlac, fileGlac);
    disp([char(39) sPath.ice char(39) ' has been chosen as the binary glacier presence grid.']);

    if isempty(fileGlac) || isequal(fileGlac, 0)
       error('CCHF_main:noGlacier',['No glacier presence grid '...
           'has been selected, even though the option was chosen.' ...
           ' Therefore, the program is aborting.']); 
    end
elseif regexpbl(sMeta.iceGrid, 'fine')
    iceAuxDisp = ['on a finer grid that will determine the ' ...
        'spatial scale of glacier process evaluation'];

    uiwait(msgbox(sprintf(['Select the glacier presence geo-referenced grid or shapefile for ' ...
        sMeta.region iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
    [fileGlac, foldGlac] = uigetfile({'*.*'},['Select the glacier presence geo-referenced grid or shapefile for ' ...
        sMeta.region iceAuxDisp '.\n'], startPath);
    sPath.ice = fullfile(foldGlac, fileGlac);
    disp([char(39) sPath.ice char(39) ' has been chosen as the glacier presence grid.']);

    if isempty(fileGlac) || isequal(fileGlac, 0)
       error('CCHF_main:noGlacier',['No glacier presence grid '...
           'has been selected, even though the option was chosen.' ...
           ' Therefore, the program is aborting.']); 
    end

    uiwait(msgbox(sprintf(['Select the glacier surface DEM for ' ...
        sMeta.region iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
    [fileGlacDem, foldGlacDem] = uigetfile({'*.*'},['Select the glacier surface DEM for ' ...
        sMeta.region iceAuxDisp '.\n'], startPath);
    sPath.iceDem = fullfile(foldGlacDem, fileGlacDem);
    disp([char(39) sPath.ice char(39) ' has been chosen as the glacier DEM grid.']);

    if isempty(fileGlacDem) || isequal(fileGlacDem, 0)
       error('CCHF_main:noGlacierDem',['No glacier DEM grid '...
           'has been selected, even though the option was chosen.' ...
           ' Therefore, the program is aborting.']); 
    end
elseif ~regexpbl(sMeta.iceGrid, 'none')
    error('CCHF_main:unknownIceGridType',['The ice grid method is ' sMeta.iceGrid ', which is not known.']);
end
    
if isempty(fileGlac) || isequal(fileGlac, 0)
   error('CCHF_main:noDebris',['No glacier presence file '...
       'has been selected. Therefore, the program is aborting.']); 
end

%Load debris cover grid (if single value, assumed uniform thickness):
if useDebris
    %Flow direction selection and display:
    uiwait(msgbox(sprintf(['Select the debris cover grid for ' ...
        sMeta.region '.\n']), '(Click OK to Proceed)','modal'));
    [fileDeb, foldDeb] = uigetfile({'*.asc';'*.txt'},['Select the debris cover grid for ' ...
        sMeta.region], startPath);
    sPath.debris = fullfile(foldDeb, fileDeb);
    disp([char(39) sPath.debris char(39) ' has been chosen as the debris cover grid.']);

    if isempty(fileDeb) || isequal(fileDeb, 0)
       error('CCHF_main:noDebris',['No glacier surface debris grid '...
           'has been selected, even though the option was chosen.' ...
           ' Therefore, the program is aborting.']); 
    end
end


%If a binary grid is being used to clip the region:
if blClip == 1
    uiwait(msgbox(sprintf(['Select the binary region grid that will be used to clip the ' ...
        sMeta.region ' region.\n']), '(Click OK to Proceed)','modal'));
    [fileRegClip, foldRegClip] = uigetfile({'*.asc';'*.txt'},['Select the binary region clip grid for ' ...
        sMeta.region], startPath);
    sPath.regionClip = fullfile(foldRegClip, fileRegClip);
    disp([char(39) sPath.regionClip char(39) ' has been chosen as the binary region clip grid.']);
    
    if isempty(fileRegClip) || isequal(fileRegClip, 0)
       error('CCHF_main:noClip',['No region clipping file '...
           'has been selected, even though the option was chosen.' ...
           ' Therefore, the program is aborting.']); 
    end
end


%Identify climate variables to load:
[sMeta.varLd, sMeta.varLdDisp] = clm_var_load(sMeta.module);
%Load climate variables:
sPath = clm_path_ui(sPath, sMeta);


%Load DEM and calculate watershed geometry fields:
sDem = read_geodata(sPath.dem, sMeta, 'none');
sHydro.dem = sDem.data;
sHydro.lat = sDem.lat;
sHydro.lon = sDem.lon;
if isfield(sPath, 'fdr')
	sFdr = read_geodata(sPath.fdr , sMeta, 'none');
    sHydro.fdrESRI = sFdr.data;
    clear('sFdr');
end
clear('sDem');

%LOAD PARAMETER FILE (IF VALIDATION OR SIMULATION RUN)
if regexpbl(runType,{'valid','sim'})
    %Load parameters from file:
    uiwait(msgbox(sprintf(['Select the set of parameter coefficients ' ...
        'written during a calibration run of the CCHF model for ' ...
        sMeta.region '.\n']), '(Click OK to Proceed)','modal'));
    [fileCoef, foldCoef] = uigetfile({'*.txt'; '*.csv'}, ['Select the set of parameter coefficients for ' ...
        sMeta.region], startPath);
    sPath.coef = fullfile(foldCoef, fileCoef);
    disp([char(39) sPath.coef char(39) ' has been chosen as the set of parameter coefficients.']);
end


%Load observation data to use for calibration or validation:
if regexpbl(runType,{'calibrat','validat','default'})
    %Calibration/validation data required:
    uiwait(msgbox(sprintf(['Select the file(s) with observation data for ' ...
        sMeta.region '.  If data does not include information such as location, you will '...
        'be required to enter that in the main panel.']), ...
        '(Click OK to Proceed)','modal'));
    
    pathGage = cell(nGage,1);
    for ii = 1 : nGage
        [fileGage, foldGage] = uigetfile({'*.*'}, ...
            ['Select observation data ' num2str(ii) ' of ' ...
            num2str(nGage) ' for ' sMeta.region], startPath);
        pathGage{ii} = fullfile(foldGage, fileGage);
    end

    %Remove any duplicate or empty elements:
    pathGage = unique(pathGage(~cellfun('isempty',pathGage)));  
   
    %Load observation data:
    sObs = read_gagedata(pathGage, sHydro.lon, sHydro.lat, ...
        'time',[sMeta.dateStart;sMeta.dateEnd], ...
        'mask',sHydro.dem);
else
    sObs = struct;
end

if ~isfield(sMeta,'output')
   sMeta.output = cell(0,2); 
end
[sMeta.output, sObs] = output_all_gage(sMeta.output, sObs, sHydro.lon, sHydro.lat);


%Ask to load information from directory of calibration run that was
%interrupted
if regexpbl(sMeta.runType,{'calib','resume'},'and')
    uiwait(msgbox(sprintf(['Select the folder containing containing '...
        'the unfinished calibration for ' sMeta.region '.\n']), ...
        '(Click OK to Proceed)','modal'));
    sPath.resume = uigetdir(startPath,['Select the folder containing '...
        'the unfinished calibration for ' sMeta.region]);
end

%Calculate fdr, fac, flow calculation order, slope, aspect, distance between neighboring centroids
sHydro = watershed(sHydro);

%%Find paremeters needed for current version of model:
%Do this by running CCHF model in 'parameter' mode
sMeta.mode = 'parameter';
sMeta.coefAtt = CCHF_engine_v3(sPath, sHydro, sMeta);
sMeta.global = global_params; %Function that contains global parameter values (albedo of ice, etc.)

%IF VALIDATION, LOAD COEFFICIENT SET AND COMPARE TO PARAMETERS NEEDED FOR CURRENT MODEL:
if regexpbl(runType,{'valid','sim'})
    if regexpbl(sPath.coef,'.txt')
        [cfTemp, hdrParam] = read_CCHF_coef(sPath.coef);

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
                    'for the ' char(39) sMeta.module char(39) ' model '...
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
            end
        end

        sMeta.coef = cfTemp;
    elseif regexpbl(sPath.coef,'.csv')
        prmArray = csvread(sPath.coef,1,0);
        fidPrm = fopen(sPath.coef,'r');
        varPrm = textscan(fidPrm,'%s',numel(prmArray(1,:)),'Delimiter',',');  
            varPrm = varPrm{1};
        fclose(fidPrm);
        
        if numel(varPrm(:)) - 1 == numel(sMeta.coefAtt(:,1))
            varPrm = varPrm{1}(1:end-1)'; %Last column is evaluation metric
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
        [~, ~, extCf] = fileparts(sPath.coef);
        error('CCHF_main:unknownCoefFormat',['The coeffienct file is in a ' ...
            char(39) extCf char(39) ' format, which has not been programmed for.']);
    end
elseif ~regexpbl(sMeta.runType,'calibrate') && numel(sMeta.coefAtt(1,:)) >= 6
    %Reduce dimensions of sMeta.coef (loaded during 'parameter' run)
   sMeta.coef = sMeta.coefAtt(:,[1,4,5]); 
end



%%START PROCESSING LOG AND DISPLAY RELEVENT CONTENT
%Create unique output directory based upon inputs:
if isfield(sPath,'resume')
    sPath.output = sPath.resume;
    [~, sMeta.strModule] = CCHF_out_dir(sPath,sMeta);
else
    [sPath.output, sMeta.strModule] = CCHF_out_dir(sPath,sMeta);
end

%Record all messages displayed on screen to file in output directory
[fileDiary] = downscale_diary(sPath.output);
disp(['All displayed text is being written to ' char(39) fileDiary ...
    char(39) '.']);
%Display modules chosen:
disp('The chosen set of module representations is: ');
for ii = 1 : numel(sMeta.module(:,1))
    disp([sMeta.module{ii,1} ' = ' sMeta.module{ii,2}]);
end
disp(blanks(1));

%Display message with citation information:
[~] = MHS_cite('','CCHF');

%Display modeling package being used and user choices:
[~,dFold1,dFold2] = fileparts(pwd);
disp([char(39) dFold1 dFold2 char(39) ' is being used.' char(10)]);
disp_CCHF_meta(sPath, sMeta)  


if blClip == 1
    sRegionClip = read_geodata(sPath.regionClip, sMeta, 'none');
    sHydro.regionClip = sRegionClip.data;
    clear('sRegionClip');
end

%create time vector to loop over based on start date and time-step
%If this field populated here, it wont be populated inside each model call,
%which saves time
sMeta = dates_run(sMeta, 'spin');
sMeta.progress = 'year';
%Find individual files to load each timestep (saves time later)
sPath = path_find_files(sPath, sMeta);    
%Record root directory in sMeta (for use during simulations)
indOutRt = regexpi(sPath.output,filesep);
sMeta.rtDir = sPath.output(1:indOutRt(end)-1);
% cd(sMeta.rtDir);


%%RUN THE HYDROLOGIC MODEL:
tStrt = now;
sOutput = cell(1,1);
if regexpbl(sMeta.runType,'calib')
    sMeta.mode = 'calibrate';
    sMeta.coef = sMeta.coefAtt;
        sMeta = rmfield_x(sMeta,'progress');
    [sOutput{1}, sMeta.coef] = CCHF_calibrate_v2(sObs, sOpt, sPath, sHydro, sMeta);
elseif regexpbl(sMeta.runType,{'valid','sim'})
    sMeta.mode = 'validate';
    
    if iscell(sMeta.coef) && ~ischar(sMeta.coef{1}) && ~all2d( size(sMeta.coef{1}) == 1 )
        sMeta = rmfield_x(sMeta,'progress');
        
        coefAll = nan(numel(sMeta.coef(:,1)), numel(sMeta.coef{1}));
        for ii = 1 : numel(coefAll(:,1))
            coefAll(ii,:) = sMeta.coef{ii,:};
        end
        
        %Open "Matlab pool" for parallel validation:
        if isempty(gcp('nocreate'))
            localCluster = parcluster('local'); %Check number of possible "workers"
            if isfield(sOpt,'maxWorker') && localCluster.NumWorkers > sOpt.maxWorker
                workers = sOpt.maxWorker;
            else
                workers = localCluster.NumWorkers;
            end
            parpool(workers); %Dedicate all available cores to parallel calibration
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
            sOutput{ii} = CCHF_engine_v3(sPath, sHydro, sMeta, coefAll(ii,:)');
        end
        
        if ~isempty(gcp('nocreate'))
           pctRunOnAll warning('on', 'all') %Turn off warnings during parfor
        end
        
%         %Close dedicated workers used for parallel processing:
%         poolobj = gcp('nocreate');
%         delete(poolobj);
    else
        sOutput{1} = CCHF_engine_v3(sPath, sHydro, sMeta);
    end
elseif regexpbl(sMeta.runType,'default')
    sMeta.mode = 'default';
    sMeta.coef = sMeta.coefAtt(:,[1,4]);
    sOutput{1} = CCHF_engine_v3(sPath, sHydro, sMeta);
else
    error('CCHF_main:runType',[char(39) sMeta.runType char(39) ...
        ' was selected as the type of model run, but this option is not recognized.']);
end
tEnd = now;
disp(['It took ' num2str(round2((tEnd-tStrt)*24*60,1)) ' minutes for the current model run.']);


if ~regexpbl(sMeta.runType,'sim')
    warning('off','all'); %Turn off warning that some time-series not being used.
    %Define additional stats to calculate on data
    statsExtra = {'KGEr', 'KGEs', 'KGEb', 'NSE', 'MAPE', 'MAE'};
    if strcmpi(fitType,'kge')
        cellStats = [fitType, statsExtra];
    else
        cellStats = [fitType, 'kge', statsExtra];
    end
    
    
    if numel(sOutput) == 1
        %Create output directory for assessment plots:
        dirModObs = fullfile(sPath.output, 'mod_v_obs_plots');
        mkdir(dirModObs);
        
        if iscell(sOutput)
            sOutput = sOutput{1};
        end

        %Display all stats and write to file
        report_stats(sObs, sOutput, cellStats, sPath.output, sMeta);

        %Create plots: modeled versus observed
        mod_v_obs_v2(sObs, sOutput, fitType, 'plot', dirModObs, 'scatter','combineType');
        mod_v_obs_v2(sObs, sOutput, fitType, 'plot', dirModObs,'combineType');
    else
        %Loop in reverse order in order to write the best performing
        %calibration parameter set last
        scoreTemp = cell(numel(sOutput),1);
        [scoreTemp{1}, typeTemp] = report_stats(sObs, sOutput{1}, cellStats, sPath.output, sMeta);
        
        scoreOut = cell(numel(scoreTemp{1}{1}),1);
        [scoreOut{:}] = deal(nan(numel(sOutput),numel(cellStats)));
        
        %Open "Matlab pool" for parallel validation:
        if isempty(gcp('nocreate'))
            localCluster = parcluster('local'); %Check number of possible "workers"
            if isfield(sOpt,'maxWorker') && localCluster.NumWorkers > sOpt.maxWorker
                workers = sOpt.maxWorker;
            else
                workers = localCluster.NumWorkers;
            end
            parpool(workers); %Dedicate all available cores to parallel calibration
        end
        
        if ~isempty(gcp('nocreate'))
           pctRunOnAll warning('off', 'all') %Turn off warnings during parfor
        end
        
        parfor mm = 2 : numel(sOutput)
            [scoreTemp{mm}, ~] = report_stats(sObs, sOutput{mm}, cellStats, sPath.output, sMeta, 'no_disp','no_write');
        end
        
        if ~isempty(gcp('nocreate'))
           pctRunOnAll warning('on', 'all') %Turn off warnings during parfor
        end
        
        %Close dedicated workers used for parallel processing:
        poolobj = gcp('nocreate');
        delete(poolobj);
        
        %Reformat statistics so that each obs type is in seperate cell
        %array
        for mm = 1 : numel(sOutput) %loop over model run
            %rows are model runs and columns are metrics
            for nn = 1 : numel(typeTemp{1}) %loop over observation type
                for ll = 1 : numel(typeTemp) %loop over metric
                    scoreOut{nn}(mm,ll) = scoreTemp{mm}{ll}(nn);
                end
            end
        end
        
        %Write all validation results to csv files (one for each obs type):
        %Create output directory for assessment plots:
        dirMult = fullfile(sPath.output, ['mod_stats_' num2str(numel(sOutput)) '_best_cal_runs']);
        mkdir(dirMult);
        
        nCol = numel(cellStats);
        headerFmt = repmat('%s,',1,nCol-1);
        numFmt = repmat('%f,',1,nCol-1);

        hdrStats = cellStats(:)';
        
        for nn = 1 : numel(scoreOut)
            pathStats = fullfile(dirMult, [typeTemp{1}{nn} '_obs.csv']);
            fStats = fopen(pathStats,'w+');
            
            fprintf(fStats, [headerFmt,'%s\n'], hdrStats{:});
            for mm = 1 : numel(scoreOut{nn}(:,1))
                fprintf(fStats, [numFmt,'%f\n'], scoreOut{nn}(mm,:));
            end
            fclose(fStats); 
        end
        
        %Calculate validation statistics and write to file:
        aggStats = {'best_cal_run','mean_val', 'high_val', 'low_val', 'median_val', 'mode_val', 'SD_val'};
        nCol = numel(cellStats) + 1;
        headerFmt = repmat('%s,',1,nCol-1);
        numFmt = repmat('%f,',1,nCol-2);

        hdrStats = [blanks(1), cellStats(:)'];
        
        for nn = 1 : numel(scoreOut)
            pathStats = fullfile(dirMult, [typeTemp{1}{nn} '_stats_4multruns.csv']);
            fStats = fopen(pathStats,'w+');
            
            fprintf(fStats, [headerFmt,'%s\n'], hdrStats{:});
            for mm = 1 : numel(aggStats)
                switch aggStats{mm}
                    case 'best_cal_run'
                        statCurr = scoreOut{nn}(1,:);
                    case 'mean_val'
                        statCurr = mean(scoreOut{nn});
                    case 'high_val'
                        statCurr = max(scoreOut{nn});
                    case 'low_val'
                        statCurr = min(scoreOut{nn});
                    case 'median_val'
                        statCurr = median(scoreOut{nn});    
                    case 'mode_val'
                        statCurr = mode(scoreOut{nn});
                    case 'SD_val'
                        statCurr = std(scoreOut{nn});
                    otherwise
                        display(['The current loop is being skipped '...
                            'because ' aggStats{mm} ' is not a recognized case.']);
                        continue
                end
                
                fprintf(fStats, [numFmt,'%f\n'], statCurr);
            end
            fclose(fStats); 
        end
        
    end
    warning('on','all'); %Turn off warning that some time-series not being used.
end


%MAKE PLOTS OF MODEL OUTPUT
if numel(sOutput) == 1
    dirModPlots = fullfile(sPath.output, 'model_output_plots');
        mkdir(dirModPlots);
    pathOutRt = fullfile(dirModPlots,'model');

    %Plot output
    plot_CCHF(sOutput, {'flow'}, sMeta, [pathOutRt '_flow']);
    % plot_CCHF(sOutput, {'et','pet'}, sMeta, [pathOutRt '_ET']);
    plot_CCHF(sOutput, {'hfnet','hfrs', 'hfrl', 'hft','hfcp', 'hfgc', 'hfsnc'}, sMeta, [pathOutRt '_heat_components'],'avg');
    % plot_CCHF(sOutput, {'prsn','rain','et','mrro'}, sMeta, [pathOutRt '_inVsOut'],'avg');  
    % plot_CCHF(sOutput, {'prsn','swe','sndwe','rain','mrro','icdwe'}, sMeta, [pathOutRt '_inVsOut'],'avg');
    plot_CCHF(sOutput, {'sndwe','icdwe','rain'}, sMeta, [pathOutRt '_inVsOut'],'avg');
    % plot_CCHF(sOutput, {'sndwe','icdwe'}, sMeta, [pathOutRt '_cryoChange']);
    % plot_CCHF(sOutput, {'swe'}, sMeta, [pathOutRt '_dSnowTmp']);
    plot_CCHF(sOutput, {'sncc'}, sMeta, [pathOutRt '_snowcoldcontent']);
    % plot_CCHF(sOutput, {'icdwe'}, sMeta, [pathOutRt '_dIce']);
    % plot_CCHF(sOutput, {'snw'}, sMeta, [pathOutRt '_snowpack']);
    plot_CCHF(sOutput, {'snlr','rnrf','iclr'}, sMeta, [pathOutRt '_water_release'],'avg');
    % plot_CCHF(sOutput, {'swe','lwsnl','snlh','snlr','iclr','lhpme','rnrf'}, sMeta, [pathOutRt '_snowpack'],'avg');
    plot_CCHF(sOutput, {'lwsnl','snlh','snlr','iclr','lhpme','sndwe'}, sMeta, [pathOutRt '_snowpack'],'avg');
    plot_CCHF(sOutput, {'tas','tsis','tsn'}, sMeta, [pathOutRt '_temperature']);
    plot_CCHF(sOutput, {'prsn','rain'}, sMeta, [pathOutRt '_rain']);
    % plot_CCHF(sOutput, {'albedoS', 'albedoI'}, sMeta, [pathOutRt '_albedo']);
    % plot_CCHF(sOutput, 'flow', sMeta, [pathOutRt '_flowrate']);
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
disp(['The ' runTypeDisp ' run for ' sMeta.region ' has finished.' ...
    char(10) 'Results have been written to ' char(39) sPath.output ...
    char(39) char(10)]);
diary off   %Stop recording command history.
