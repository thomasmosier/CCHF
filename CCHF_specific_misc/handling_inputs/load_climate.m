function [varargout] = load_climate(sPath, sHydro, sMeta, varargin)
global sAtm

% sMeta.varLd = {'pr','tas','tmn','tmx'};

if isequal(sMeta.dt, sMeta.dataTsRes)
    sAtm = load_ts(sPath, sHydro, sMeta.dateCurr, sMeta);
elseif regexpbl(sMeta.dt,{'day','daily'}) && regexpbl(sMeta.dataTsRes,'month') %Load adjacent months to current time step and then inperpolate monthly data to daily time-steps
    %Find adjacent dates (based on time step):
    dateAdj = adj_dates(sMeta.dataTsRes, sMeta.dateCurr);
    
    %LOAD INPUT TIME-SERIES DATA:
    if ~isempty(varargin) && ~isempty(fieldnames(varargin{1}))
        sInput = varargin{1}; 

        %Find if currently loaded time-series data has current time-step
        [dateRefData, ~] = NC_time_units(sInput.attTime);
        datesData = days_2_date(sInput.time, dateRefData, NC_cal(sInput.attTime));

        %If current time-step not present in already interpolated data,
        %interpolate:
        if sum(ismember(datesData, sMeta.dateCurr, 'rows')) == 0
            sInput = load_ts(sPath, sHydro, dateAdj, sMeta, sInput);
        end
    else
        sInput = load_ts(sPath, sHydro, dateAdj, sMeta);
    end


    
    if sum(~isnan(dateAdj(2,:))) > 1
        datesInterp = [dateAdj(2,1)*ones(eomday(dateAdj(2,1),dateAdj(2,2)),1), dateAdj(2,2)*ones(eomday(dateAdj(2,1),dateAdj(2,2)),1), (1:eomday(dateAdj(2,1),dateAdj(2,2)))'];
    else
        datesInterp = nan(1,3);
    end
    
    [sInput.pre, timeInterp] = quad_fit_exact(sInput, 'pre', datesInterp, 1);
    [sInput.tmp,          ~] = quad_fit_exact(sInput, 'tmp', datesInterp, 0);

    if isfield(sPath,'tmin') && ~isempty(sPath.tmin)
        [sInput.tmn,~] = quad_fit_exact(sInput,'tmn', datesInterp, 0);
    end
    if isfield(sPath,'tmax') && ~isempty(sPath.tmax)
        [sInput.tmx,~] = quad_fit_exact(sInput, 'tmx', datesInterp, 0);
    end 
    
    sInput.time = timeInterp;
    
    sAtm = sInput;
else
    error('backbone:day2XDay',['The case where the model time-step is '...
        sMeta.dt ' and the input data are ' sMeta.dataTsRes ...
        ' has not been coded for.']);
end

%Only create varargout if requested
if nargout > 0
    if ~isequal(sMeta.dt, sMeta.dataTsRes)
        varargout{1} = sInput;
    else
        varargout{1} = sAtm;
    end
end

%Check that sAtm spatial grid aligns with DEM grid:
if ~isempty(fieldnames(sAtm)) && (~isequal(round2(sAtm.lon,4), round2(sHydro.lon,4)) || ~isequal(round2(sAtm.lat,4), round2(sHydro.lat,4)))
    warning('load_climate:diffSpatialGrid',['The spatial grid of the '...
        'DEM being used and the climate inputs are different. The '...
        'cliamte data is being spatially interpolated, which '...
        'significantly slows down the model.']);
    for ii = 1 : numel(sMeta.varLd)
        gridTemp = nan([numel(sAtm.time),numel(sHydro.lat),numel(sHydro.lon)],'single');
        for jj = 1 : numel(sAtm.time)
            gridTemp(jj,:,:) = PCHIP_2D(sAtm.lon, sAtm.lat, squeeze(sAtm.(sMeta.varLd{ii})(jj,:,:)) , sHydro.lon, sHydro.lat);
        end
        
        if regexpbl(sMeta.varLd{ii},'pr')
            gridTemp(sAtm.(sMeta.varLd{ii}) < 0) = 0;
        end
        
        sAtm.(sMeta.varLd{ii}) = gridTemp;
    end
end


