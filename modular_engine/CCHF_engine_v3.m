% Copyright 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the CCHF Hydrologic Modelling Package (referred to 
% as the 'CCHF Package').
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

function varargout = CCHF_backbone_v3(sPath, sHydro, sMeta, varargin)


%'sMeta' is a structure array describing data properties:
    %'currVar' = current variable to load
    %'metVar' = Cell array of all variables requested
    %'metVarDisp' = Cell array of display names for all variables
        %requested
    %'hisTS' = 'GPCC', 'CRU', or 'Willmott'
    %'hisTSDisp' = full version of entry in 'sMeta.hisTS'
    %'yrsOut' = [first year of output data, last year of output data]
    %'yrsClim' = [first year of climatology, last year of climatology]
    %'mnths' = list of months requested
    %'crd' = [lonW, lonE, latS, latN]
    %'dateCurr' = [year, month, day, hour]  - subdivisons depend on time
        %step
%'sHydro' is a structure array describing downscaling characteristics:
    %'sHydro.dem' =  
    %'sHydro.fdr' =  
    %'glacier' =  
    %'lat' = 
    %'lon' = 
%'sPath' is a sctructure array with all the paths to necessary data inputs:
%(Many of these fields are cells, allowing them to store all path strings if 
%multiple variables being downscaled)
    %'hisTS' = path of 0.5 degree monthly time-series input
    %'hisGCM' = path of historical GCM run (only used in downscaling GCM
        %data)
    %'hrClim' = path of 30-arcsecond climatology (i.e. WorldClim)
    %'futTS' = path of projected GCM run (only used in downscaling GCM
        %data)
    %'regClip' = path of boolean data used to clip downscaled outputs (can
        %be empty)
    %'output' = output path

    
%Declare global variables (must clear them because model run in parallel)
if ~regexpbl(sMeta.runType,'sim') %Don't clear land structure array if running simulation (may be dangerous but saves multiple hours for large grids)
    clear global sLand %Possibly simply initialize instead of clearing
end

clear global sAtm sCryo
global sAtm sCryo sLand
sLand.indCurr = 1; %Necessary if sLand not cleared above.
    
    
%%Set default values:
%This value simply determines when to test the fitness of the current
%parameter set. Early convergence is allowed if does not meet fitness
%requirement. Default value is huge, then is edited below.
indEval = 10^6;
sObs = struct;
fitTest = 'KGE';
fitReq = 100;
blEval = 0;
modError = nan;


%If coefficient parameter set given as variable input argument, override
%parameters in sMeta.
if ~regexpbl(sMeta.mode,'parameter')
    if ~isfield(sMeta,'dateRun')
        sMeta = dates_run(sMeta,'spin');
    end
    
    if ~isempty(varargin(:))
        %First variable argument is coefficient set:
        if isnumeric(varargin{1}) && all(~isnan(varargin{1}))
            varargin{1} = num2cell(squeeze(varargin{1}));
        end

        if  isfield(sMeta,'coefAtt') && numel(sMeta.coefAtt(:,1)) == numel(squeeze(varargin{1}))
            sMeta.coef = [sMeta.coefAtt(:,1), varargin{1}(:)];
        elseif isfield(sMeta,'coef') && numel(sMeta.coef(:,1)) == numel(squeeze(varargin{1})) && ischar(sMeta.coef{1})
            sMeta.coef = [sMeta.coef(:,1), varargin{1}(:)];
        else
            error('CCHF_backbone:coefVar',['Variable input argument being ' ...
                'used to assign parameter values, but neither ' char(39) ...
                'sMeta.coef' char(39) ' nor ' char(39) 'sMeta.coefAtt' ...
                char(39) ' appear to contain the necessary parameter name '...
                'fields.']);
        end
        %Second is structure, which will be used for representative analysis over short period. 
        if numel(varargin(:)) > 1
            if numel(varargin(:)) == 3
                if isstruct(varargin{2}) && iscell(varargin{3})
                    sObs = varargin{2};
                    evalAtt = varargin{3};
                elseif isstruct(varargin{3}) && iscell(varargin{2})
                    sObs = varargin{3};
                    evalAtt = varargin{2};
                else
                    error('backbone:varargTyoe',['Variable argument 2 '...
                        'must be structure of observations and 3 must '...
                        'be cell array with evaluation attributes.']);
                end

                fitTest = find_att(evalAtt, 'fitTest');
                dateEval = find_att(evalAtt, 'dateEval');
                fitReq = find_att(evalAtt, 'fitReq');
                if numel(sMeta.dateRun(1,:)) > numel(dateEval)
                    dateEval = [dateEval, ones(1, numel(sMeta.dateRun(1,:)) - numel(dateEval))];
                elseif numel(sMeta.dateRun(1,:)) < numel(dateEval)
                    dateEval = dateEval(1:numel(sMeta.dateRun(1,:)));
                end
                indEval = find(ismember(sMeta.dateRun, dateEval, 'rows') == 1);
                if isempty(indEval)
                    strDateEval = blanks(0);
                    for kk = 1 : numel(dateEval)
                        strDateEval = [strDateEval, num2str(dateEval(kk))];
                        if kk ~= numel(dateEval)
                            strDateEval = [strDateEval, '-'];
                        end
                    end
                   error('backbone:dateMismatch',['The dateEval, ' ...
                       strDateEval ', does not match any of the dates the model will be run for.']); 
                end
            else
               error('backbone:nVarargin',['There are ' ...
                   num2str(numel(varargin{:})) ' extra input arguments but 3 are expected.']);
            end
        end
        
    elseif isfield(sMeta,'coef')
        if numel(sMeta.coef(1,:)) >= 2 && numel(sMeta.coef(1,:)) <= 4
            sMeta.coef = sMeta.coef;
        else
            error('CCHF_backbone:coefHydroDim',[char(39) 'sMeta.coef' char(39)...
                ' has unexpected dimenions.']);
        end
    elseif isfield(sMeta,'coefAtt')
        if numel(sMeta.coefAtt(1,:)) >= 5
            sMeta.coef = sMeta.coefAtt(:,[1,4,5]);
        end
    else
        error('CCHF_backbone:noCoef1',['No input field containing '...
            'the requisite parameters was found.']);
    end
    
    if isempty(sMeta.coef)
        error('CCHF_backbone:noCoef2',['An input field was found but '...
            'its format was not correct.']);
    end
end

if regexpbl(sMeta.mode,'parameter')
    disp(['The CCHF model is being run in ' char(39) sMeta.mode char(39) ' mode.']);
    %Define date vector of time steps to loop over:
    if regexpbl(sMeta.dt,'month')
        sMeta.dateRun = nan(1,2);
    elseif regexpbl(sMeta.dt,{'day','daily'})
        sMeta.dateRun = nan(1,3);
    elseif regexpbl(sMeta.dt,'hour')
        sMeta.dateRun = nan(1,4);
    else
        error('backbone:dtUnknown',[sMeta.dt ' is an unknown time-step']);
    end
    
    %Remove 'progress' field from sMeta if parameter run.
    sMeta = rmfield_x(sMeta, 'progress');
    
    %Initialize coef argument for assigning parameters to:
    coef = cell(0,6);
else    
    %Determine when to start recording output:
    sMeta.indWrite = find(all(sMeta.dateRun == repmat(sMeta.dateStart, numel(sMeta.dateRun(:,1)),1),2) == 1);

    %Create array to hold output(s):
    sModOut = struct;
    
    %Determine if climate input files have been clipped and aligned to
    %input DEM (if so, code field that will expedite loading process):
    for ll = 1 : numel(sMeta.varLd) 
        if ~isfield(sPath,sMeta.varLd{ll}) || isempty(sPath.(sMeta.varLd{ll}))
            continue
        end
        
        if regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},'.nc')
            lonTemp = ncread(sPath.([sMeta.varLd{ll} 'File']){1},'lon')';
            latTemp = ncread(sPath.([sMeta.varLd{ll} 'File']){1},'lat');
        elseif regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},{'.asc','.txt'})
            [~,hdrESRI,hdrMeta] = read_ESRI(sPath.([sMeta.varLd{ll} 'File']){1});
            [latTemp, lonTemp] = ESRI_hdr2geo(hdrESRI,hdrMeta);
        else
            error('backbone:unkownFileType','The present file type has not been coded for.');
        end
            
        if isequal(round2(latTemp,4),round2(sHydro.lat,4)) && isequal(round2(lonTemp,4),round2(sHydro.lon,4 ))
            if ~isfield(sPath,'varNm')
                sPath.varNm = cell(size(sMeta.varLd));
            end
            if isempty(sPath.varNm{ll})
                if regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},'.nc')
                    nmVarTemp = ncinfo(sPath.([sMeta.varLd{ll} 'File']){1});
                    nmVarTemp = squeeze(struct2cell(nmVarTemp.Variables(:)));
                    nmVarTemp = nmVarTemp(1,:)';
                        nmVarTemp(strcmpi(nmVarTemp,'lon')) = [];
                        nmVarTemp(strcmpi(nmVarTemp,'longitude')) = [];
                        nmVarTemp(strcmpi(nmVarTemp,'lat')) = [];
                        nmVarTemp(strcmpi(nmVarTemp,'latitude')) = [];
                        nmVarTemp(strcmpi(nmVarTemp,'time')) = [];
                        nmVarTemp(strcmpi(nmVarTemp,'time_bnds')) = [];
                    if numel(nmVarTemp) == 1
                        sPath.varNm{ll} = nmVarTemp{1};
                    end
                elseif regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},{'.asc','.txt'})
                    if regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},{'pr'})
                        varTmp = 'pr';
                    elseif regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},{'tmin','tmn','tasmin'})
                        varTmp = 'tasmin';
                    elseif regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},{'tmax','tmx','tasmax'})
                        varTmp = 'tasmax';
                    elseif regexpbl(sPath.([sMeta.varLd{ll} 'File']){1},{'tmean','tmp'})
                        varTmp = 'tas';
                    else
                        error('backbone:varDetect',['The variable associated with the file: ' ...
                            sPath.([sMeta.varLd{ll} 'File']){1} ' is unknown.']);
                    end
                    sPath.varNm{ll} = varTmp;
                end
            end
        end
    end
end

%progress update option in sMeta:
if isfield(sMeta,'progress')
    if regexpbl(sMeta.progress,{'year','yr'})
        indEnd = 1;
    elseif regexpbl(sMeta.progress,{'month','mnth'})
        indEnd = 2;
    elseif regexpbl(sMeta.progress,{'day','dy'})
        indEnd = 3;
    elseif regexpbl(sMeta.progress,{'hour','hr'})
        indEnd = 4;
    else
        warning('backbone:progressUnknown',['No progress update will be given because ' ...
            sMeta.progress ' is not recognized.'])
    end
    
    uniqTs = unique(sMeta.dateRun(:,1:indEnd),'rows');
            
    sMeta.indUpdate = nan(numel(uniqTs(:,1)),1);
    for zz = 1 : numel(uniqTs(:,1))
        sMeta.indUpdate(zz) = find(ismember(sMeta.dateRun(:,1:indEnd), uniqTs(zz,1:indEnd), 'rows') == 1, 1, 'first');
    end
end



szDem = size(sHydro.dem);
nTS  = numel(sMeta.dateRun(:,1));


%Enter time loop:
for ii = 1 : nTS
    sMeta.dateCurr = sMeta.dateRun(ii,:);
    sMeta.indCurr = ii;
%     disp(['iteration ' num2str(ii) ' of ' num2str(numel(sMeta.dateRun(:,1)))]);

    if isfield(sMeta, 'indUpdate')
        if ismember(ii, sMeta.indUpdate)
            disp(['The hydrologic model is running time step '...
                num2str(ii) ' of ' num2str(nTS) '.']);
        end
    end
    
    %%LOAD TIME-SERIES CLIMATE DATA:
    if isequal(sMeta.dt,sMeta.dataTsRes) %climate time step and model time step are same
        load_climate_v2(sPath, sHydro, sMeta);
    else %time steps different (in some cases, this may be able to 
        if ii == 1
           sInput = struct; 
        end
        sInput = load_climate_v2(sPath, sHydro, sMeta, sInput);
    end
    

    %Initialize cryosphere and land arrays:
    if ii == 1   
        %Define edges of main grid and record in sMeta (useful for limiting data loaded using 'geodata'):
        edgeLat = box_edg(sHydro.lat);
        edgeLon = box_edg(sHydro.lon);
        
        sMeta.crd = [edgeLon(1), edgeLon(end), edgeLat(end), edgeLat(1)]; %[W, E, S, N]
        
        
        %Initialize cryosphere structure:
        sCryo = rmfield_x(sAtm,{'data','pr','pre','tas','tasmin','tmin','tmn','tmax','tasmax','tmx','lat','lon'});

        %INITIALIZE ICE GRID
        %Load glacier coverage data and convert to depth:
        %Define seperate ice spatial grid if path defined
        sHydro.icbl = ice_grid_init(sHydro, sPath, sMeta);

        %Initialize solid snow depth (water equivalent) grid
        sCryo.snw = zeros(size(sHydro.dem),'single');
        
        %initialize snow internal temperature:
        if regexpbl(sMeta.coef, 'tsn')
            sCryo.tsn  = nan(size(sCryo.snw),'single'); 
            sCryo.tsis = nan(size(sCryo.snw),'single');
        end

        
        
        %%INITIALIZE LAND ARRAY (for vegetation cover - if used, soil
        %moisture, runoff, flow, etc.)
        sLand = rmfield_x(sAtm,{'data','pr','pre','tas','tasmin','tmin','tmn','tmax','tasmax','tmx','lat','lon','time','attTime','indCurr'});

        sLand.mrro = zeros(szDem,'single'); %runoff out of each cell
        sLand.flow   = zeros(szDem,'single'); %flowrate at each cell
        sLand.sm     = zeros(szDem,'single'); %Soil moisture field
        sLand.et     = zeros(szDem,'single'); %evapotranspiration field
    end
    
%%DO EVERY ITERATION:
    %Reset some fields every iteration:
    sLand.mrro = zeros(szDem,'single'); %runoff out of each cell
    
    if isfield(sCryo,'ice')
        sCryo.icdwe = zeros(szDem,'single');
    else
        sCryo.icdwe = nan(szDem,'single');
    end
    sCryo.sndwe = zeros(szDem,'single');

    
    %UPDATE TIME INDICES for sAtm (potentially used when loading climate
    %inputs)
    if ~regexpbl(sMeta.mode,'parameter')
        [sAtm.indCurr, ~] = geodata_time_ind(sAtm, sMeta.dateRun(ii,:));
    end
    
%%FITTING PARAMETERS FOR PRE AND TMP:
%     if regexpbl(sMeta.mode,'parameter')
%         coef = cat(1,coef,{'cfPre', 0.1, 10, 4.0, 'backbone','input'});
%     else
%         cfPre = find_att(sMeta.coef,'cfPre');
%         sAtm.pr = cfPre*sAtm.pr;
%     end
%     if regexpbl(sMeta.mode,'parameter')
%         coef = cat(1,coef,{'cfTmp', -5, 5, 1.50, 'backbone','input'});
%     else
%         cfTmp = find_att(sMeta.coef,'cfTmp');
%         sAtm.tas = cfTmp + sAtm.tas;
%     end
    
    %Allow katabatic scaling of temperature over glaciers:
        %Katabatic wind (cooling effect over glaciers; see Shea and Moore, 
        %2010 or Alaya, Pellicciotti, and Shea 2015 for more complex models):
%     if regexpbl(sMeta.mode,'parameter')
%         coef = cat(1,coef,{'cfKat', -6, 0, -2.75, 'backbone','input'});
%     else
%         cfKat = find_att(sMeta.coef,'cfKat');
% 
%         indIce = find(sCryo.ice > 0);
%         szTmp = size(sAtm.tas);
%         [tIce,~] = meshgrid((1:szTmp(1)),(1:numel(indIce)));
%         [rIce, cIce] = ind2sub(szTmp(2:3),indIce);
%         [~,rIce] = meshgrid((1:szTmp(1)),rIce);
%         [~,cIce] = meshgrid((1:szTmp(1)),cIce);
%         indIce = tIce + (rIce-1)*szTmp(1) + (cIce-1)*szTmp(2)*szTmp(1);
% 
%         sAtm.tas(indIce) = sAtm.tas(indIce) + cfKat;
%     end
    


%%IMPLEMENT MODULES:
    if regexpbl(sMeta.mode,'parameter')
        coefTemp = cat(1,coef, module_implement(sHydro, sMeta));
        %In case there are coefficients with dependencies, need to
        %iterate coefficient search.
        %Currently the only case of this is snowpack temperature, which is
        %dealt with here: 
        if regexpbl(coefTemp(:,1), 'tsn')
            sCryo.tsn = nan(size(sHydro.dem)); %This line is necessary to access some snow temperature dependent fitting parameters
            coef = cat(1,coef, module_implement(sHydro, sMeta));
            indRem = strcmpi(coef(:,1),'tsn');
            coef(indRem,:) = []; 
        else
            coef = coefTemp;
        end
    else
        module_implement(sHydro, sMeta);
    end



    
%%FIND AND WRITE OUTPUT FIELDS:
    if isfield(sMeta,'output') && ~regexpbl(sMeta.mode,'parameter') && ii >= sMeta.indWrite
       sMeta.indCurr = ii;
       sModOut = print_model_state_v2(sModOut, sHydro, sMeta, sPath.output); 
    end
    
%%ALLOW EARLY TERMINATION IF CALIBRATION AND REPRESENTATIVE PERIOD HAS BEEN
%%SURVEYED
    if ii > indEval && ~isempty(fieldnames(sObs)) && blEval == 0
        blEval = 1;
        modErrorTemp = nanmean(mod_v_obs_v2(sObs, sModOut, fitTest, 'combineType', 'rmMod'));
        if modErrorTemp > fitReq %Quit early
            modError = modErrorTemp;
            break
        end
    end
end

if regexpbl(sMeta.mode,'parameter')
    varargout{1} = coef(~cellfun(@isempty, coef(:,1)),:);
else
    varargout{1} = sModOut;
    if nargout > 1 %Will be nan if short analysis not option or if met fitness requirement during short analysis.
        varargout{2} = modError;
    end
end
