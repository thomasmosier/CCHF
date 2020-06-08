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


%Do hydrologic calculations only during validation or simulation modes:
for kk = 1 : nSites
%     %Stats of inerest:
%     Runoff and streamflow
%         Total WY runoff contributions from each grid cell (rain, snow
%         melt, ice melt) (m^3/year). Gridded
%         Date of max runoff contributions
%     SWE 
%         Total WY snowfall (m^3). Gridded
%         Maximum WY SWE (m^3). Gridded
%         Date of max SWE (day of WY). Gridded
%     Glacier mass balance
%         Annual change in glacier mass balance (kg/yr). Gridded
%         Total basin change in glacier mass balance (kg/yr). Lumped.
    %It is contained in: sMod{kk}.all
    
    
    %Find unique water years:
    if isfield(sMod{kk}, 'all') 
        if isfield(sMod{kk}.all, 'date')
            indWYBeg = find(ismember(sMod{kk}.all.date(:,2:3), [10,1], 'rows'));
                indWYBeg = indWYBeg(:); %Ensure column vector
            if ~isempty(indWYBeg)
                indWYEnd = [indWYBeg(2:end) - 1; numel(sMod{kk}.all.date(:,1))];
            else
                warning('hydroStats:waterYearNoWY', ['Hydrologic statistics for ' ...
                sMeta.region{kk} ' are being skipped because no water years have been found.']);
            end
        else
            warning('hydroStats:waterYearNoDate', ['Hydrologic statistics for ' ...
                sMeta.region{kk} ' are being skipped because there is no date field.']);
        end
    else
        warning('hydroStats:waterYearNoDate', ['Hydrologic statistics for ' ...
                sMeta.region{kk} ' are being skipped because there is no gridded field.']);
    end
    
    WYTest = nanmean(indWYEnd(1:end-1) - indWYBeg(1:end-1));
    if indWYEnd(end) - indWYBeg(end) < WYTest - 2 || indWYEnd(end) - indWYBeg(end) > WYTest + 2
        indWYBeg = indWYBeg(1:end-1);
        indWYEnd = indWYEnd(1:end-1);
    end
    
    WYYr = sMod{kk}.all.date(indWYEnd,1);
    
    %Create direcotry for output:
    dirHydroStats = fullfile(sPath{kk}.output,'hydrologic_stats');
        mkdir(dirHydroStats);

        
    %%Runoff and streamflow:
    varLp = {'mrro', 'rnrf', 'snlr', 'iclr', 'snw'};
    varUnit = {'m/s','m/s','m/s','m/s', 'm'};
    
    %Loop over variables:
    for jj = 1 : numel(varLp)
        if isfield(sMod{kk}.all, varLp{jj}) 
            %Loop over water years
            for ii = 1 : numel(indWYBeg)
                [grdWYMax, grdWYDOY] = max(sMod{kk}.all.(varLp{jj})(indWYBeg(ii):indWYEnd(ii),:,:), [], 1);
                    grdWYMax = squeeze(grdWYMax);
                    grdWYDOY = squeeze(grdWYDOY); grdWYDOY(isnan(grdWYMax)) = nan;

                grdWYSum = squeeze(nansum(sMod{kk}.all.(varLp{jj})(indWYBeg(ii):indWYEnd(ii),:,:)));
                    grdWYSum(isnan(grdWYMax)) = nan;

                %Convert from seconds to years:
                if strcmpi(varUnit{jj}, 'm/s')
                    secToDy = 86400;
    %                 secToYr = (indWYEnd(ii) - indWYBeg(ii) + 1) * secToDy;
                    grdWYMax = grdWYMax * secToDy;
                    grdWYSum = grdWYSum * secToDy;

                    unitMax = 'm/day';
                    unitSum = 'm/yr';
                elseif ~strcmpi(varUnit{jj}, 'm')
                   error('hydroStats:unknownUnits', ['Units ' varUnit{jj} ' for variable ' varLp{jj} ' are not recognize.']); 
                else
                    unitMax = varUnit{jj};
                    unitSum = varUnit{jj};
                end


                pathOutSum = fullfile(dirHydroStats, [varLp{jj} '_total_WY' num2str(WYYr(ii))]);
                pathOutMax = fullfile(dirHydroStats, [varLp{jj} '_max_WY' num2str(WYYr(ii))]);
                pathOutDOY = fullfile(dirHydroStats, [varLp{jj} '_dayofyr-max_WY' num2str(WYYr(ii))]);

                if regexpbl(sMeta.wrtTyp, {'nc', 'netcdf'})
                    dateVec  = nan(size(sMod{kk}.all.date(1,:)));
                    dateBnds = [sMod{kk}.all.date(indWYBeg(ii),:); sMod{kk}.all.date(indWYEnd(ii),:)];

                    warning('off', 'printGridNc:nanTime');
                    print_grid_NC_v2([pathOutSum, '.nc'],   grdWYSum, varLp{jj}, sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds, unitSum, 1);
                    print_grid_NC_v2([pathOutMax, '.nc'],   grdWYMax, varLp{jj}, sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds, unitMax, 1);
                    print_grid_NC_v2([pathOutDOY, '.nc'],   grdWYDOY, varLp{jj}, sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds,   'day', 0);
                    warning('on', 'printGridNc:nanTime');
                elseif regexpbl(sMeta.wrtTyp, 'asc')
                    hdrWrt = ESRI_hdr(sHydro{kk}.(varLon), sHydro{kk}.(varLat), 'corner');
                    write_ESRI_v4(grdWYSum, hdrWrt, [pathOutSum, '.asc'], 1);
                    write_ESRI_v4(grdWYMax, hdrWrt, [pathOutMax, '.asc'], 1);
                    write_ESRI_v4(grdWYDOY, hdrWrt, [pathOutDOY, '.asc'], 1);
                end
            end
            clear ii
        else
            warning('hydroStats:unknownVar', ['The variable ' varLp{jj} ...
                ' does not exist in the all field for region ' sMeta.region{kk} ...
                ' and is therefore being skipped.']);
        end
    end
    clear jj
    
    
    %%Ice mass balance:
    if isfield(sMod{kk}.all, 'icmb') && ndims(sMod{kk}.all.icmb) == 3 && numel(sMod{kk}.all.icmb(:,1,1)) > 0
        %Loop over ice mass balance years
        for ii = 1 : numel(indWYBeg)
            icmbWY = squeeze(sMod{kk}.all.icmb(ii,:,:));

            pathOutIcmb = fullfile(dirHydroStats, ['icmb' '_WY' num2str(WYYr(ii))]);

            if regexpbl(sMeta.wrtTyp, {'nc', 'netcdf'})
                dateVec  = nan(size(sMod{kk}.all.date(1,:)));
                dateBnds = [sMod{kk}.all.icmb_dateStart(ii,:); sMod{kk}.all.icmb_dateEnd(ii,:)];

                warning('off', 'printGridNc:nanTime');
                print_grid_NC_v2([pathOutIcmb, '.nc'], icmbWY, 'icmb', sHydro{kk}.(varLon), sHydro{kk}.(varLat), dateVec, dateBnds, 'm/yr', 1);
                warning('on', 'printGridNc:nanTime');
            elseif regexpbl(sMeta.wrtTyp, 'asc')
                hdrWrt = ESRI_hdr(sHydro{kk}.(varLon), sHydro{kk}.(varLat), 'corner');
                write_ESRI_v4(icmbWY, hdrWrt, [pathOutIcmb, '.asc'], 1);
            end
        end
        clear ii
    else
        warning('hydroStats:unknownVar', ['The variable ' 'icmb' ...
                ' does not exist in the all field for region ' sMeta.region{kk} ...
                ' and is therefore being skipped.']);
    end
end