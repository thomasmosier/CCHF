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



%Path for all inputs (file with stored paths is created during model run in
%which user selects inputs through GUI)
path2Inputs = '';

%Path to calibration state file generated during previous, interupted
%calibration run (leave empty unless " runType = 'calibrate_resume' "
pathCalibrateResume = ''; 


%%USER INPUTS:
%Select the type of model run:
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
    'heat',  'TI_Kraaijenbrink'; ... %surface heat flux module representation
    'snmass', 'cc'; ... %snow energy and mass representation
    'icmlt', 'ratio'; ... %ice melt representation
    'runoff', 'bucket'; ... %runoff representation
    'toa', 'DeWalle'; ... %top-of-atmosphere radiation representation
    'atmtrans', 'dem_decay'; ... %atmospheric transmissivity representation
    'snalbedo', 'Pelli'; ... %snow albedo representation
    'icalbedo', 'constant'; ... %ice albedo representation
    'glacier0', 'external'; ... %glacier water equivalent at time 0
    'partition', 'ramp'; ... %Partitioning precip into snow and rain
    'pet', 'Hammon'; ... %potential evapotranspiration representation
    'timelag', 'Johnstone'; ... %timelag (of flow through each grid cell) representation
    'flow', 'lumped'; ... %flow routing representation
    'density', 'Liston'; ... %densification of snow (does not participate in firn to ice processes)
    'sndrain', 'percent'; ... %snow holding capacity representation
    'sca', 'max'; ... %Snow covered area representation
    'sublimate', 'Lutz'; ... %Sublimation representation
    'firn', 'threshold'; ... %firn compaction representation
    'glaciervel', 'sheer'; ... %glacier movement representation
	'glaciermove', 'slide'; ... %glacier movement representation
    'avalanche', 'angle' %snow avalanching (can be 'none')
    };

%Boolean value to toggle on all glacier processes. If off, glaciers can
%still be included, but they will be static (infinite sources of melt water)
blDynamicGlaciers = 0; %0 = static glaciers, 1 = dynamic (use firn/glacier processes selected in module set) 

%These settings only used in calibration (can to left alone for validation)
optType = 'apso';  %optimization method: 
                        %'GA_binary' (Genetic Algorithm using binary chromosomes), 
                        %'GA_real' (genetic algorithm using real numbers), 
                        %'monte_carlo', 
                        %'uniform_sampling', 
                        %'PSO' (Particle Swarm optimization), or 
                        %'hybrid' (combination of PSO, Monte Carlo, and linear sensitivity)
                        %'apso' (Adaptive PSO)
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
%Wolverine:
testPts = {[-148.921, 60.4197]; [-148.915, 60.3769]; [-148.907, 60.4042]};

reportPt = {'avg'};

output = {...
    'flow', {'avg'; testPts{:}}; ... %surface flowrate through cell
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

   
%Values assigned here will overwrite existing values from other sources 
%(and will not be calibrated)
%Format is {name of function, parameter name, value}
%Leave empty if no parameters should be set
knownCfValues = {...
    'atmtrans_dem_decay', 'trans_clear_intercept', 0.7105; ...
    'atmtrans_dem_decay', 'trans_clear_dem', 0.0235; ...
    'atmtrans_dem_decay', 'trans_decay_rate', 3.4318; ...
    'atmtrans_dem_decay', 'trans_decay_power', 1.6142; ...
    'partition_ramp', 'tmp_snow', -1; ...
    'partition_ramp', 'tmp_rain', 3 ...
    }; 
%Values for 'atmtrans' derived from optimization relative to ERA
%Values for 'partition_ramp' derived from visual inspection of figure in 
    %Dai, A. (2008). Temperature and pressure dependence of the rain-snow 
    %phase transition over land and ocean
    

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
                            
useIcePond = 0; %Set to 1 to load glacier lake cover array. If 0, presumed there are no glacier lakes

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
%%%%%ASSIGN INPUTS TO STRUCTURE ARRAYS:    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    sMeta.('iceGrid')  = iceGrid;
    sMeta.('glacierDynamics') = blDynamicGlaciers;
    sMeta.varLd        = {'pr','tas','tasmin','tasmax'};
    sMeta.('blDebris') = useDebris;
    sMeta.('blIcePond') = useIcePond;
    sMeta.('blDispOut') = blDispOutput;
    sMeta.('fixedPrm') = knownCfValues;
    sMeta.('addoutput') = output;
    sMeta.('maxworker') = maxWorkersVal;
    sMeta.('ldFdr') = useFdr;
    sMeta.('nGage') = nGage;
    sMeta.('useprevrun') = blPrevRun;
    sMeta.('pathinputs') = path2Inputs;
    sMeta.('pathresume') = pathCalibrateResume;

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
    nGen = round(nRunsMax / nPop);

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
    if regexpbl(sOpt.type, {'genetic','GA'})
        sOpt.stagnate = 20;
    elseif regexpbl(sOpt.type, {'PSO','hybrid'})
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%IMPLEMENTATION CCHF MODEL:    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sMod, sObs] = CCHF_implement(sMeta, sOpt);