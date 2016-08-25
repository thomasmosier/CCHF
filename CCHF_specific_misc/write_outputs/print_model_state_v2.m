function sModOut = print_model_state_v2(sModOut, sHydro, sMeta, pathOut)
%persistent sModOut
global sAtm sCryo sLand


szMain = size(sHydro.dem);

%Determine if ice grid different from main grid:
if isfield(sCryo,'iceDem') && ~isequal(szMain,size(sCryo.iceDem))
    iceGrid = 1;
else
    iceGrid = 0;
end

%Don't write during spinup, only write if output requested, and if
%model not bieng run in parameter retrieve mode
if sMeta.indCurr < sMeta.indWrite
    return
else
    indCurrTs = sMeta.indCurr - sMeta.indWrite + 1;
    
    foldGrids = fullfile(pathOut, 'grids');
        
    %On first iteration, create output structure array:
    if indCurrTs == 1 || isempty(sModOut)
        outDate = date_vec_fill(sMeta.dateStart, sMeta.dateEnd, 'Gregorian');
        dataInit = nan(numel(outDate(:,1)),1);
        
        %If ice grid different, but define variables to place output within
        %ice grid:
        if iceGrid == 1
            edgLat = box_edg(sCryo.iceLat);
            edgLon = box_edg(sCryo.iceLon);
            szIce = size(sCryo.iceDem);
        end
        
        %Loop over all variables to output
        for kk = 1 : numel(sMeta.output(:,1)) 
            %loop over points where model output should be written:
            for ll = 1 : numel(sMeta.output{kk,3}(:)) 
                %Switch between types of output points ('avg', 'all', or
                %specific point)
                if strcmpi(sMeta.output{kk,3}{ll},'avg')
                    if ~isfield(sModOut,'avg')
                        sModOut.avg = struct('lon', nan, 'lat', nan, ...
                            'date', outDate, sMeta.output{kk,1}, dataInit, ...
                            'indGage', (1:numel(sHydro.dem)));
                        sModOut.avg.fields = sMeta.output(kk,1);
                    else
                        sModOut.avg.(sMeta.output{kk,1}) = dataInit;
                        sModOut.avg.fields{end+1} = sMeta.output{kk,1};
                    end
                elseif strcmpi(sMeta.output{kk,3}{ll},'all')
                    if ~isfield(sModOut,'all')
                        sModOut.all = struct('lon', sHydro.lon, 'lat', sHydro.lat, ...
                            'date', outDate, sMeta.output{kk,1}, dataInit, ...
                            'indGage', (1:numel(sHydro.dem)));
                        sModOut.all.fields = sMeta.output(kk,1);
                    else
                        sModOut.all.(sMeta.output{kk,1}) = dataInit;
                        sModOut.all.fields{end+1} = sMeta.output{kk,1};
                    end
                    
                    if any(strcmpi(sMeta.output{kk,1},{'casi','snw','snx','icx'}))
                        sModOut.all.(sMeta.output{kk,1}) ...
                            = nan(numel(outDate(:,1)),numel(sHydro.lat),numel(sHydro.lon),'single');
                    end
                    
                    if sMeta.wrtAll == 1 && ~regexpbl(sMeta.runType,'calib')
                        sModOut.all.(['path' sMeta.output{kk,1}]) = fullfile(foldGrids, sMeta.output{kk,1});
                    end
                else
                    ptWrt = ['pt' num2str(sMeta.output{kk,3}{ll})];
                    if ~isfield(sModOut, (ptWrt))
                        [r,c] = ind2sub(szMain,sMeta.output{kk,3}{ll});
                        sModOut.(ptWrt) ...
                            = struct('lon', sHydro.lon(c), 'lat', sHydro.lat(r), ...
                            'date', outDate, sMeta.output{kk,1}, dataInit, ...
                            'indGage', sMeta.output{kk,3}{ll});
                        sModOut.(ptWrt).fields = sMeta.output(kk,1);
                    else
                        sModOut.(ptWrt).(sMeta.output{kk,1}) = dataInit;
                        sModOut.(ptWrt).fields{end+1} = sMeta.output{kk,1};
                    end
                    
                    %The indices for ice points are potentially wrong because ice
                    %can have different grid than main
                    if iceGrid == 1
                        %Find ice grid cell for current output
                        %Lat:
                        rIce = find(sMeta.output{kk,2}{ll}(2) <= edgLat(1:end-1) & sMeta.output{kk,2}{ll}(2) >= edgLat(2:end));
                        %Lon:
                        cIce = find(sMeta.output{kk,2}{ll}(1) >= edgLon(1:end-1) & sMeta.output{kk,2}{ll}(1) <= edgLon(2:end));
                        
                        if isempty(rIce) || isempty(cIce)
                           error('print_model_state:iceOutside',['A model output point for ' sMeta.output{kk,1} ' falls outside the ice grid.']); 
                        end
                        indIce = sub2ind(szIce, rIce, cIce);
                        sModOut.(ptWrt).iceLon = sCryo.iceLon(cIce);
                        sModOut.(ptWrt).iceLat = sCryo.iceLat(rIce);
                        sModOut.(ptWrt).indIce = indIce;
                    end
                end
            end
        end
    end %End of initialization
       
    %Get Names of all points where output data requested:
    namesOut = fieldnames(sModOut);
    
    indDate = sMeta.indCurr - sMeta.indWrite + 1;
    
    %Populate output structure for each of the points:
    for kk = 1 : numel(namesOut(:)) %Loop over grid pts to write to
        for ll = 1 : numel(sModOut.(namesOut{kk}).fields) %Loop over data fields for current grid point
            nmCurr = char(sModOut.(namesOut{kk}).fields{ll});
            switch char(sModOut.(namesOut{kk}).fields{ll})
                case 'flow'
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC(fileNm, sLand.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                    else  
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sLand.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'mrro' %Total runoff from each grid cell
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC(fileNm, sLand.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                    else  
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sLand.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'sm'
                    if isfield(sLand,'sm')
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sLand.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sLand.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'pet'
                    if isfield(sLand,'pet')
                        indPET = find(ismember(sLand.datepet(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
                        szPET = size(sAtm.pr);
                        if isempty(indPET)
                            error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
                        end
                        [gageRow, gageCol] = ind2sub(szMain,sModOut.(namesOut{kk}).indGage);
                        indGage3 = indPET + (gageRow-1)*szPET(1) + (gageCol-1)*szPET(2)*szPET(1);
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            %NOT PROGRAMMED YET
                            if indCurrTs == 1
                                warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                            end
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sLand.(nmCurr)(indGage3).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'et'
                    if isfield(sLand,'et')
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, 'et', [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sLand.(nmCurr), 'et', sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sLand.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'pr'
                    szPre = size(sAtm.pr);
                    [gageRow, gageCol] = ind2sub(szMain,sModOut.(namesOut{kk}).indGage);
                    indGage3 = sAtm.indCurr + (gageRow-1)*szPre(1) + (gageCol-1)*szPre(2)*szPre(1);
                    if strcmpi(namesOut{kk},'all')
                        %NOT PROGRAMMED YET
                        if indCurrTs == 1
                            warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                        end
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sAtm.(nmCurr)(indGage3).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'rain'
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC(fileNm, sAtm.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sAtm.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'rnrf' %Rain that is added to runoff immediately (i.e. no snow beneath)
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC(fileNm, sLand.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sLand.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'prsn'
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC(fileNm, sAtm.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sAtm.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'snw' %snow water equivalent
                    if strcmpi(namesOut{kk},'all') 
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs,:,:) = sCryo.(nmCurr);

                        %Writing to file not programmed yet
                        if isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        end
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'scx' %snow covered area
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC(fileNm, sCryo.nmCurr, nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'icx' %snow covered area
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'casi' %snow/ice covered area
                    if isfield(sCryo,'icdbr')
                        %If debris cover field, only count non-debris covered ice
                        if strcmpi(namesOut{kk},'all') 
                            if isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                                %NOT PROGRAMMED YET
                                if indCurrTs == 1
                                    warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                                end
                            else
                                %NOT PROGRAMMED YET
                                if indCurrTs == 1
                                    warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                                end
                            end
                        else
                            casiTemp = max(sCryo.scx(sModOut.(namesOut{kk}).indGage),sCryo.icx(sModOut.(namesOut{kk}).indGage));
                            casiTemp(sCryo.icdbr > 0) = 0;
                            
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = ...
                                round(100*sum2d(casiTemp.*sHydro.area(sModOut.(namesOut{kk}).indGage)) ...
                            /sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage)));
                        end
                    else %No debris cover field
                        if strcmpi(namesOut{kk},'all') 
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs,:,:) = ...
                                    round(100*max(sCryo.scx, sCryo.icx));
                                
                            %Writing to file not programmed yet
                            if isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                                %NOT PROGRAMMED YET
                                if indCurrTs == 1
                                    warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                                end
                            end
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = ...
                                round(100*sum2d(max(...
                                sCryo.scx(sModOut.(namesOut{kk}).indGage),sCryo.icx(sModOut.(namesOut{kk}).indGage)...
                                ).*sHydro.area(sModOut.(namesOut{kk}).indGage)) ...
                                /sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage)));
                        end
                    end
                case 'rsdt' %Top of atmosphere solar radiation
                    if isfield(sAtm,'rsdt')
                        if isfield(sAtm,'rsdtdate')
                            indRSDT = find(ismember(sAtm.datersdt(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
                            if isempty(indRSDT)
                                error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
                            end
                            if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                                fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                                print_grid_NC(fileNm, squeeze(sAtm.(nmCurr)(indRSDT,:,:)), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                            else
                                temp = squeeze(sAtm.(nmCurr)(indRSDT,:,:));
                                sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(temp(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage)); 
                            end
                        elseif numel(size(sAtm.(nmCurr)))  == 2
                            if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                                fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                                print_grid_NC(fileNm, sAtm.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                            else
                                sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sAtm.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage)); 
                            end
                        else
                            error('print_model_state:rsdt','The field sAtm.rsdt has unexpected dimensions.')
                        end
                    end
                case 'sncc'
                    if isfield(sCryo,'sncc')
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage)); 
                        end
                    end
                case 'tas'
                    szTmp = size(sAtm.(nmCurr));
                    [gageRow, gageCol] = ind2sub(szMain,sModOut.(namesOut{kk}).indGage);
                    indGage3 = sAtm.indCurr + (gageRow-1)*szTmp(1) + (gageCol-1)*szTmp(2)*szTmp(1);
                    if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                        %NOT PROGRAMMED YET
                        if indCurrTs == 1
                            warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                        end
                    else
                        sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sAtm.(nmCurr)(indGage3).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                    end
                case 'tsn' %Snow internal temperature
                    if isfield(sCryo,nmCurr)
                        indUse = find(~isnan(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage)));
                        if ~isempty(indUse)
                            if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                                fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                                print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                            else
                                sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage(indUse)).*sHydro.area(sModOut.(namesOut{kk}).indGage(indUse)))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage(indUse)));
                            end
                        end
                    end
                case 'tsis' %Snow/ice surface temperature
                    if isfield(sCryo, nmCurr)
                        indUse = find(~isnan(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage)));
                        if ~isempty(indUse)
                            if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                                fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                                print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                            else
                                sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage(indUse)).*sHydro.area(sModOut.(namesOut{kk}).indGage(indUse)))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage(indUse)));
                            end
                        end
                    end
                case 'sndt' %Change in snowpack temperature
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'icdt' %Change in ice temperature
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'lwsnl' %Liquid content of snow
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'snlh' %Liquid holding capacity of snow
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'snlr' %runoff from snow melt
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'iclr' %runoff from ice melt
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && ~regexpbl(sMeta.runType,'calib')
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'icdwe'
                    if iceGrid == 1 %CHANGE
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})]) %Print grid to file
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            if numel(sModOut.(namesOut{kk}).indGage) == sum2d(~isnan(sHydro.dem)) %Finding domain average
                                sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area);
                            else %Report ice change at specific point
                                indUse = intersect(sModOut.(namesOut{kk}).indGage, find(sCryo.icx ~= 0));
                                if ~isempty(indUse)
                                    sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage(indUse)).*sHydro.area(sModOut.(namesOut{kk}).indGage(indUse)))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage(indUse)));
                                    sModOut.(namesOut{kk}).icbl(indCurrTs) = 1;
                                else
                                    sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = 0;
                                    sModOut.(namesOut{kk}).icbl(indCurrTs) = 0;
                                end
                            end
                        end 
                    else
                        if isfield(sCryo,'icbl')
                            if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})]) %Print grid to file
                                fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                                print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                            else
                                indUse = intersect(sModOut.(namesOut{kk}).indGage, find(sCryo.icbl == 1));
                                if numel(sModOut.(namesOut{kk}).indGage) == sum2d(~isnan(sHydro.dem)) %Finding domain average
                                    sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr).*sHydro.area)/sum2d(sHydro.area);
                                else %Report ice change at specific point
                                    sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(indUse).*sHydro.area(indUse))/sum2d(sHydro.area(indUse));
                                end
                                
                                %Indicate if ice is present at current grid
                                %cell
                                if ~isempty(indUse)
                                    sModOut.(namesOut{kk}).icbl(indCurrTs) = 1;
                                else
                                    sModOut.(namesOut{kk}).icbl(indCurrTs) = 0;
                                end
                            end 
                        else
                            if indCurrTs == 1
                                warning('print_model_state:noIce','Ice changes are not being recorded because no ice grid is being used in the current model run.')
                            end
                        end
                    end
                case 'sndwe'
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    elseif isfield(sModOut.(namesOut{kk}),'snowpack') %If 'sndwe' is not field, can calculate snowfall from previous record
                        if sMeta.indCurr == sMeta.indWrite
                            prevSnow = 0;
                        else
                            prevSnow = sModOut.(namesOut{kk}).(nmCurr)(indCurrTs-1);
                        end
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.snw - prevSnow, nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.snw(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage)) - prevSnow;
                        end
                    end
                case 'lhpme' %latent heat melt equivalent (units of m SWE)
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hfnet'
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hfrl' %longwave radiation flux
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hfrs' %Shortwave radiation
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hfcs' %Sensible convective
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hfcp' %latent flux from precipitation
                    if isfield(sCryo, nmCurr) 
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hfsnc' %conduction through snow (surface to internal)
                    if isfield(sCryo, nmCurr) 
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hfgc' %Conduction from ground upwards
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.nmCurr, nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'hft' %Heat flux from degree index factor
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end
                case 'snalb'
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end  
                    if isfield(sCryo, nmCurr)
                        if strcmpi(namesOut{kk},'all') && isfield(sModOut.(namesOut{kk}),['path' char(sModOut.(namesOut{kk}).fields{ll})])
                            fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(namesOut{kk}).fields{ll}, sModOut.(namesOut{kk}).date(indDate,:)) '.nc']);
                            print_grid_NC(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, 1);
                        else
                            sModOut.(namesOut{kk}).(nmCurr)(indCurrTs) = sum2d(sCryo.(nmCurr)(sModOut.(namesOut{kk}).indGage).*sHydro.area(sModOut.(namesOut{kk}).indGage))/sum2d(sHydro.area(sModOut.(namesOut{kk}).indGage));
                        end
                    end  
                otherwise
                    if indCurrTs == 1
                        warning('backbone:unknownOutput',[char(39) ...
                            char(sModOut.(namesOut{kk}).fields{ll}) char(39) ...
                            ' is not a known output type.']);
                    end
            end
        end
    end
end