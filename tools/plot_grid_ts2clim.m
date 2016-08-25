clear all

%Read and plot gridded SETI output

% var = 'pr';

foldPlot = uigetdir(pwd,'Select the folder containing gridded NetCDF data to plot.');
indSep = regexpi(foldPlot,filesep);

% if isempty(var)
%     varLd = foldPlot(indSep(end)+1:end);
% else
%     varLd = var;
% end
    
%Find files to plot
filesPlot = dir([foldPlot, filesep, '*', '.nc']);
    filesPlot = extractfield(filesPlot,'name');
    
    
%Open waitbar:
hWait = waitbar(0,'Processing starting.');  
%%READ DATA    
cntr = 0;
for ii = 1 : numel(filesPlot)
    waitbar(ii/numel(filesPlot), hWait, ...
            ['Reading file ' num2str(ii) ' of ' num2str(numel(filesPlot)) '.']);
    
%     if mod(ii,500) == 0
%        disp(['The script is on the ' num2str(ii) ' of ' num2str(numel(filesPlot)) ' iteration.']) 
%     end
    pathCurr = fullfile(foldPlot,filesPlot{ii});
    
    if ii == 1
        latVec = ncread(pathCurr,'lat');
        lonVec = ncread(pathCurr,'lon')';
%         areaGrid = area_geodata(lonVec,latVec);
    end
    
    
    
    sDataCurr.time = ncread(pathCurr,'time');
    
    if ii == 1
        finfo = ncinfo(pathCurr);
        allVar = extractfield(finfo.Variables,'Name');
            allVar(strcmpi(allVar,'lon')) = [];
            allVar(strcmpi(allVar,'longitude')) = [];
            allVar(strcmpi(allVar,'lat')) = [];
            allVar(strcmpi(allVar,'latitude')) = [];
            allVar(strcmpi(allVar,'time')) = [];
            allVar(strcmpi(allVar,'time_bnds')) = [];
        if numel(allVar) == 1
            varLd = allVar{1};
        else
            error('plot_gridded:autoVarFail','Cannot auto-detect variable to load.');
        end
        
        vecDays = nan(numel(filesPlot),numel(sDataCurr.time));
    end
    
    sDataCurr.attTime = ncinfo(pathCurr,'time');
    sDataCurr.attTime = squeeze(struct2cell(sDataCurr.attTime.Attributes))';
    try
        sDataCurr.attData = ncinfo(pathCurr,varLd);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:imagesci:netcdf:unknownLocation'))
%             if ~isempty(var)
%                 varLd = foldPlot(indSep(end)+1:end);
%             end
            
            finfo = ncinfo(pathCurr);
            allVar = extractfield(finfo, 'Variables');
                allVar(strcmpi(allVar,'lon')) = [];
                allVar(strcmpi(allVar,'longitude')) = [];
                allVar(strcmpi(allVar,'lat')) = [];
                allVar(strcmpi(allVar,'latitude')) = [];
                allVar(strcmpi(allVar,'time')) = [];
                allVar(strcmpi(allVar,'time_bnds')) = [];
            if numel(allVar) == 1
                varLd = allVar{1};
            else
                error('plot_gridded:autoVarFail','Cannot auto-detect variable to load.');
            end
        else
            error('plot_gridded:cannotFindVar','Cannot detect variable to load.');
        end
        sDataCurr.attData = ncinfo(pathCurr, varLd);
    end
    
    sDataCurr.attData = squeeze(struct2cell(sDataCurr.attData.Attributes))';
    
%     indUnit = strcmpi(sDataCurr.attData,'units');
%         [rUnit, cUnit] = find(indUnit == 1);
%         dataUnit = sDataCurr.attData{rUnit, cUnit+1};

    sDataCurr.data = ncread(pathCurr, varLd);
    
    vecDays(ii,:) = sDataCurr.time;
    
    for jj = 1 : sum(~isnan(vecDays(ii,:)))
        indTime = find(~isnan(vecDays(ii,:)));
        cntr = cntr + 1;

        [dateRef, Units] = NC_time_units(sDataCurr.attTime);
        cal = find_att(sDataCurr.attTime,'calendar');

        if ii == 1
            nMnth = zeros(12,1);
            sInput.all.(varLd) = zeros([12,size(squeeze(sDataCurr.data(1,:,:)))],'single');
            
            sInput.all.date = nan(numel(vecDays(:)),numel(dateRef));
                sInput.all.lon = NaN;
                sInput.all.lat = NaN;
        end
        
        sInput.all.date(cntr,:) = days_2_date(vecDays(ii,indTime(jj)),dateRef,cal);
        nMnth(sInput.all.date(cntr,2)) = nMnth(sInput.all.date(cntr,2)) + 1;
        sInput.all.(varLd)(sInput.all.date(cntr,2),:,:) = squeeze(sInput.all.(varLd)(sInput.all.date(cntr,2),:,:)) ...
            + squeeze(sDataCurr.data(indTime(jj),:,:));
    end
end 
close(hWait);

for ii = 1 : 12
    sInput.all.(varLd)(ii,:,:) = squeeze(sInput.all.(varLd)(ii,:,:)) / nMnth(ii);
end

foldOut = fullfile(foldPlot,'clim');
mkdir(foldOut);

unqYrs = unique(sInput.all.date(:,1));

for ii = 1 : 12
    pathWrt = fullfile(foldOut,[varLd '_' num2str(unqYrs(1)) 'thru' num2str(unqYrs(end)) '_' num2str(ii) '.nc']);
    print_grid_NC(pathWrt, squeeze(sInput.all.(varLd)(ii,:,:)), varLd, lonVec, latVec, [2000,ii,15]);
%     
%     sData.lon = lonVec;
%     sData.lat = latVec;
%     sData.data = squeeze(sInput.all.(varLd)(ii,:,:));
%     sData.var = varLd;
%     sData.time = nan;
%     sData.attTime = {'type','unknown'};
%     sMeta.region = 'ci';
%     sMeta.currTime = nan(1,2);
%     write_geodata(pathWrt, sData, sMeta, 0, 'nc');
end

% dtTest = diff(vecDays);
% dtTest(dtTest < 0) = [];
% 
% dtTest = nanmean(round(dtTest));
% if dtTest > 0.4 && dtTest < 5
%     sMeta.dt = 'daily';
% elseif  dtTest > 26 && dtTest < 32
%     sMeta.dt = 'monthly';
% else
%     error('plotGridded:tStepUnknown','The time step is not known.');
% end
% 
% %Write Spatially integrated results to text file:
% pathOutputYr = fullfile(foldOut, [varLd 'spatial_integrated_origTS.txt']);
% write_SETI_gagedata(pathOutputYr, sInput);
%     sMeta.dateStart = sInput.all.date(1,:);
%     sMeta.dateEnd   = sInput.all.date(end,:);
% plot_SETI(sInput, {varLd}, sMeta, fullfile(foldOut, [varLd '_origTS']));
% 
% 
% %%SORT BY YEAR
% indUse = find(~isnan(sInput.all.date(:,1)));
% yrs = unique(sInput.all.date(indUse,1));
% sYrSum.all.date = [yrs, 6*ones(numel(yrs),1), 1*ones(numel(yrs),1)];
%     sYrSum.all.lon = NaN;
%     sYrSum.all.lat = NaN;
% sYrSum.all.(varLd) = nan(numel(sYrSum.all.date(:,1)),1);
% 
% for ii = 1 : numel(sYrSum.all.date(:,1))
%     indYr = find(sInput.all.date(:,1) == sYrSum.all.date(ii,1));
%     sYrSum.all.(varLd)(ii) = sum(sInput.all.(varLd)(indYr));
% end
% 
% %Write Spatially integrated results to text file:
% pathOutputYr = fullfile(foldOut, [varLd 'spatial_integrated_yrSum.txt']);
% write_SETI_gagedata(pathOutputYr, sYrSum);
%     sMeta.dateStart = [sYrSum.all.date(1,1), 1, 1];
%     sMeta.dateEnd   = [sYrSum.all.date(end,1), 12, 31];
% plot_SETI(sYrSum, {varLd}, sMeta, fullfile(foldOut, [varLd '_yrSum']));
% 
% 
% %%AVERAGE EACH MONTH:
% sMnthCL.all.date = [ones(12,1), (1:12)', 15*ones(12,1)];
%     sMnthCL.all.lon = NaN;
%     sMnthCL.all.lat = NaN;
% 
% sMnthCL.all.(varLd) = nan(numel(sMnthCL.all.date(:,1)),1);
% for ii = 1 : numel(sMnthCL.all.date(:,1))
%     indMnth = find(sInput.all.date(:,2) == sMnthCL.all.date(ii,2));
%     sMnthCL.all.(varLd)(ii) = sum(sInput.all.(varLd)(indMnth))/numel(sYrSum.all.date(:,1));
% end
% 
% %Write Spatially integrated results to text file:
% pathOutputMnth = fullfile(foldOut, [varLd 'spatial_integrated_mnthCL_' ...
%     num2str(sMnthCL.all.date(1)) '-' num2str(sMnthCL.all.date(end)) '.txt']);
% write_SETI_gagedata(pathOutputMnth, sMnthCL);
%     sMeta.dateStart = [1,1,1];
%     sMeta.dateEnd   = [1,12,31];
% plot_SETI(sMnthCL, {varLd}, sMeta, fullfile(foldOut, [varLd '_mnthCL']));
% 
% %Write annual average:
% pathOutputYrCL = fullfile(foldOut, [varLd 'spatial_integrated_annCL_' ...
%     num2str(yrs(1)) '-' num2str(yrs(end)) '.txt']);
% sYrCL = sMnthCL;
% sYrCL.all.(varLd) = nansum(sMnthCL.all.(varLd));
% sYrCL.all.date = nan(1,3);
% write_SETI_gagedata(pathOutputYrCL, sYrCL);


% for ii = 1 : numel(filesPlot)
%     pathCurr = fullfile(foldPlot,filesPlot{ii});
%     
%     if ii == 1
%         latVec = ncread(pathCurr,'lat');
%         lonVec = ncread(pathCurr,'lon')';
%         areaGrid = area_geodata(lonVec,latVec);
%     end
%     
%     sDataCurr.time = ncread(pathCurr,'time');
%     sDataCurr.attTime = ncinfo(pathCurr,'time');
%         sDataCurr.attTime = squeeze(struct2cell(sDataCurr.attTime.Attributes))';
%     sDataCurr.attData = ncinfo(pathCurr,varLd);
%         sDataCurr.attData = squeeze(struct2cell(sDataCurr.attData.Attributes))';
% %     indUnit = strcmpi(sDataCurr.attData,'units');
% %         [rUnit, cUnit] = find(indUnit == 1);
% %         dataUnit = sDataCurr.attData{rUnit, cUnit+1};
%     sDataCurr.data = ncread(pathCurr,varLd);
%     
%     currSum = sum2d(areaGrid.*squeeze(sDataCurr.data(ii,:,:)));
% end