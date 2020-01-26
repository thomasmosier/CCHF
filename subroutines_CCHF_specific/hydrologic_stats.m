function hydrologic_stats(sHydro, sPath, sMod, sMeta)


%Find number of regions:
if iscell(sMeta.region)
    nSites = numel(sMeta.region(:));
elseif ischar(sMeta.region)
    nSites = 1;
    sMeta.region = {sMeta.region};
else
    error('HpatImplement:nmRegions',['The model does not recognize the '...
        'format of the region input variable.']);
end

%Find dates that model was run for:
sMeta = dates_run(sMeta);

%Define variables:
varLon = 'longitude';
varLat = 'latitude';
varPathFlow = 'flow_path';

%Turn off plotting warning
warning('off', 'plot_CCHF:noData');

keyboard
%Do hydrologic calculations only during validation or simulation modes:
for kk = 1 : nSites
    %Load/find flow data written to file during model runs:
    flow = nan([numel(sMeta.dateRun(:,1)), size(sHydro{kk}.dem)], 'single');
    dirFlow = '';
    if iscell(sMod) && numel(sMod(:)) == 1
        if isfield(sMod{kk}.all, 'flow')
            flow = sMod{kk}.all.flow;
        elseif isfield(sMod{kk}.all, varPathFlow)
            dirFlow = sMod{1}.all.(varPathFlow);
        end
    elseif isstruct(sMod)
        if isfield(sMod.all, 'flow')
            flow = sMod.all.flow;
        elseif isfield(sMod{kk}.all, varPathFlow)
            dirFlow = sMod.all.(varPathFlow);
        end
    else
        warning('HpatImplement:unexpectedOutputType','No output plots will be produced because the output are of an unexpected type.');
    end

    %Load flow from files:
    if all(isnan(flow(:)))
        if ~isempty(dirFlow)
            filesFlow = find_files(dirFlow,'nc');
            if ~isempty(filesFlow(:))
                for ii = 1 : numel(sMod{kk}.all.date(:,1))
                    pathCurr = fullfile(dirFlow, filesFlow{ii});
                    flow(ii,:,:) = ncread(pathCurr, 'flow');
                end
            else
                error('HpatImplement:noNcFlow', ['No NetCDF files were found '...
                    'containing the gridded flow data necessary to calculate hydropower potential statistics.']);
            end
        else
            warning('HpatImplement:missingFlowPath', ['No hydropower calculations ' ...
                'can be produced because the path information to modelled '...
                'streamflow information is missing.']);
        end
    end
    
    
    if sum(~isnan(flow(:))) == 0
        error('HpatImplement:noModelledFlow', ['The entire flow array is nan. ' ...
            'Therefore no hydropower calculations can be conducted. This indicates '...
            'that the modelled flow was not stored for the entire grid array either in memory or in files.']);
    end

    flowAvg = nanmean(flow,1);

    %CALCULATE HYDROPOWER POTENTIAL
    if iscell(sMod) && numel(sMod(:)) == 1
        [powerAvg, powerRho, powerSd, powerSm] = power_potential(sHydro{kk}.slopeFdr, sHydro{kk}.dlFdr, flow, sMod{1}.all.date);
    elseif isstruct(sMod)
        [powerAvg, powerRho, powerSd, powerSm] = power_potential(sHydro{kk}.slopeFdr, sHydro{kk}.dlFdr, flow, sMod.all.date);
    else
        warning('HpatImplement:unexpectedOutputType','No output plots will be produced because the output are of an unexpected type.');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%% DEFINE AND IMPLEMENT RANKING METRIC
    %%%%%%%%%%%%%%%%%%%%%%%%
    powerQual = power_quality(powerRho, powerSm);
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%% END OF RANKING METRIC
    %%%%%%%%%%%%%%%%%%%%%%%%


    %Write hydropower potential statistics to files:
    dirHydroStats = fullfile(sPath{kk}.output,'hydropower_stats');
        mkdir(dirHydroStats);
    pathFlowAvg = fullfile(dirHydroStats, 'flow_mean');
    pathPwrAvg  = fullfile(dirHydroStats, 'power_mean');
    pathPwrRho  = fullfile(dirHydroStats, 'power_density');
    pathPwrSd   = fullfile(dirHydroStats, 'power_standard_deviation');
    pathPwrSm   = fullfile(dirHydroStats, 'power_stability_metric');
    pathPwrQual = fullfile(dirHydroStats, 'power_quality');

    if regexpbl(sMeta.wrtTyp, {'nc', 'netcdf'})
        dateVec  = nan(size(sMeta.dateRun(1,:)));
        dateBnds = [sMeta.dateRun(1,:); sMeta.dateRun(end,:)];

        warning('off', 'printGridNc:nanTime');
        print_grid_NC_v2([pathFlowAvg, '.nc'],   flowAvg,      'flow_avg', sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds,    'm^3/s', 1);
        print_grid_NC_v2([ pathPwrAvg, '.nc'],  powerAvg,       'pwr_avg', sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds,        'W', 1);
        print_grid_NC_v2([ pathPwrRho, '.nc'],  powerRho,   'pwr_density', sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds,      'W/m', 2);
        print_grid_NC_v2([  pathPwrSd, '.nc'],   powerSd,        'pwr_sd', sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds,        'W', 2);
        print_grid_NC_v2([  pathPwrSm, '.nc'],   powerSm, 'pwr_stability', sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds, 'unitless', 2);
        print_grid_NC_v2([pathPwrQual, '.nc'], powerQual,   'pwr_quality', sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds,      'W/m', 2);
        warning('on', 'printGridNc:nanTime');
    elseif regexpbl(sMeta.wrtTyp, 'asc')
        hdrWrt = ESRI_hdr(sHydro{kk}.(varLon), sHydro{kk}.(varLat), 'corner');
        write_ESRI_v4(  flowAvg, hdrWrt, [pathFlowAvg, '.asc'], 1);
        write_ESRI_v4( powerAvg, hdrWrt, [ pathPwrAvg, '.asc'], 1);
        write_ESRI_v4( powerRho, hdrWrt, [ pathPwrRho, '.asc'], 2);
        write_ESRI_v4(  powerSd, hdrWrt, [  pathPwrSd, '.asc'], 2);
        write_ESRI_v4(  powerSm, hdrWrt, [  pathPwrSm, '.asc'], 2);
        write_ESRI_v4(powerQual, hdrWrt, [pathPwrQual, '.asc'], 2);
    end
end