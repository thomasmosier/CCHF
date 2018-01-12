% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function varargout = print_model_state_v4(sHydro, sMeta, pathOut)

%Use global and persistent variles.
%Persistent variable 'sModOut' speeds run time by close to factor of 100
global sAtm sCryo sLand
persistent sModOut

fldsCryo = [fieldnames(sCryo); 'casi'; 'stake'; 'geodetic'];
fldsAtm  = fieldnames(sAtm);
fldsLand = fieldnames(sLand);

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
    indTsPrintCurr = sMeta.indCurr - sMeta.indWrite + 1;
    
    foldGrids = fullfile(pathOut, 'model_output_grids');
        
    %On first iteration, create output structure array:
    if indTsPrintCurr == 1 || isempty(sModOut)
        sModOut = struct;
        outDate = date_vec_fill(sMeta.dateStart, sMeta.dateEnd, 'Gregorian');
        szInitPt = [numel(outDate(:,1)),1];
        szInit3d = [numel(outDate(:,1)),numel(sHydro.lat),numel(sHydro.lon)];
        
        %If ice grid different, but define variables to place output within
        %ice grid:
        if iceGrid == 1
            edgLat = box_edg(sCryo.iceLat);
            edgLon = box_edg(sCryo.iceLon);
            szIce = size(sCryo.iceDem);
        end
        

        %Loop over all variables to output
        for kk = 1 : numel(sMeta.output{sMeta.siteCurr}(:,1)) 
            %loop over points where model output should be written:
            for ll = 1 : numel(sMeta.output{sMeta.siteCurr}{kk,3}(:)) 
                %Switch between types of output points ('avg', 'all', or
                %specific point)
                if strcmpi(sMeta.output{sMeta.siteCurr}{kk,3}{ll},'avg')
                    if ~isfield(sModOut,'avg') || isempty(sModOut.avg)
                        sModOut.avg = struct;
                            sModOut.avg.lon = nan;
                            sModOut.avg.lat = nan;
                            sModOut.avg.date = outDate;
                            sModOut.avg.indGage = (1:numel(sHydro.dem));
                            sModOut.avg.fields = cell(0,1);
%                         sModOut.avg = struct('lon', nan, 'lat', nan, ...
%                             'date', outDate, 'indGage', (1:numel(sHydro.dem)), 'fields', cell(0,1));
                    end
                    sModOut.avg.(sMeta.output{sMeta.siteCurr}{kk,1}) = nan(szInitPt);
                    sModOut.avg.fields{end+1} = sMeta.output{sMeta.siteCurr}{kk,1};

                elseif strcmpi(sMeta.output{sMeta.siteCurr}{kk,3}{ll},'all')
                    if ~isfield(sModOut,'all') || isempty(sModOut.all)
                        sModOut.all = struct;
                            sModOut.all.lon = sHydro.lon;
                            sModOut.all.lat = sHydro.lat;
                            sModOut.all.date = outDate;
                            sModOut.all.indGage = (1:numel(sHydro.dem));
                            sModOut.all.fields = cell(0,1);
%                         sModOut.all = struct('lon', sHydro.lon, 'lat', sHydro.lat, ...
%                             'date', outDate, 'indGage', (1:numel(sHydro.dem)), 'fields', cell(0,1));
                    end
                    sModOut.all.(sMeta.output{sMeta.siteCurr}{kk,1}) = nan(szInit3d, 'single');
                    sModOut.all.fields{end+1} = sMeta.output{sMeta.siteCurr}{kk,1};
                    
%                     if any(strcmpi(sMeta.output{sMeta.siteCurr}{kk,1},{'casi','snw','snx','icx'}))
%                         sModOut.all.(sMeta.output{sMeta.siteCurr}{kk,1}) ...
%                             = nan(numel(outDate(:,1)),numel(sHydro.lat),numel(sHydro.lon),'single');
%                     end
                    
                    if sMeta.wrtGridOut == 1 && ~regexpbl(sMeta.runType,'calib')
                        sModOut.all.([sMeta.output{sMeta.siteCurr}{kk,1} '_path']) = fullfile(foldGrids, sMeta.output{sMeta.siteCurr}{kk,1});
                    end
                else
                    ptWrt = ['pt' num2str(sMeta.output{sMeta.siteCurr}{kk,3}{ll})];
                    
                    %Check that the point is inside the current domain:
                    lonCurr = sMeta.output{sMeta.siteCurr}{kk,2}{ll}(1);
                    latCurr = sMeta.output{sMeta.siteCurr}{kk,2}{ll}(2);

                    lonStep = nanmean(abs(diff(sHydro.lon)));
                    latStep = nanmean(abs(diff(sHydro.lat)));
                    lonEdg = [min(sHydro.lon) - 0.5*lonStep, max(sHydro.lon) + 0.5*lonStep];
                    latEdg = [max(sHydro.lat) + 0.5*latStep, min(sHydro.lat) - 0.5*latStep];

                    if lonCurr >= lonEdg(1) && lonCurr <= lonEdg(2) && latCurr <= latEdg(1) && latCurr >= latEdg(2) 
                        if ~isfield(sModOut, ptWrt) || isempty(sModOut.(ptWrt)) 
                            [r,c] = ind2sub(szMain,sMeta.output{sMeta.siteCurr}{kk,3}{ll});
                            sModOut.(ptWrt) = struct;
                                sModOut.(ptWrt).lon = sHydro.lon(c);
                                sModOut.(ptWrt).lat = sHydro.lat(r);
                                sModOut.(ptWrt).date = outDate;
                                sModOut.(ptWrt).indGage = sMeta.output{sMeta.siteCurr}{kk,3}{ll};
                                sModOut.(ptWrt).fields = cell(0,1);
%                             sModOut.(ptWrt) ...
%                                 = struct('lon', sHydro.lon(c), 'lat', sHydro.lat(r), ...
%                                 'date', outDate, 'indGage', sMeta.output{sMeta.siteCurr}{kk,3}{ll}, 'fields', cell(0,1));

                            %This may not be needed...
                            if isfield(sMeta, 'siteCurr')
                               sModOut.(ptWrt).site = sMeta.siteCurr;
                            end
                        end
                        
                        sModOut.(ptWrt).(sMeta.output{sMeta.siteCurr}{kk,1}) = nan(szInitPt);
                        sModOut.(ptWrt).fields{end+1} = sMeta.output{sMeta.siteCurr}{kk,1};
                    end
                   
                    
                    %The indices for ice points are potentially wrong because ice
                    %can have different grid than main
                    if iceGrid == 1
                        %Find ice grid cell for current output
                        %Lat:
                        rIce = find(sMeta.output{sMeta.siteCurr}{kk,2}{ll}(2) <= edgLat(1:end-1) & sMeta.output{sMeta.siteCurr}{kk,2}{ll}(2) >= edgLat(2:end));
                        %Lon:
                        cIce = find(sMeta.output{sMeta.siteCurr}{kk,2}{ll}(1) >= edgLon(1:end-1) & sMeta.output{sMeta.siteCurr}{kk,2}{ll}(1) <= edgLon(2:end));
                        
                        if isempty(rIce) || isempty(cIce)
                           error('print_model_state:iceOutside',['A model output point for ' sMeta.output{sMeta.siteCurr}{kk,1} ' falls outside the ice grid.']); 
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
       
    
        
    indDate = sMeta.indCurr - sMeta.indWrite + 1;
    
    
    %Get names of all points where output data requested:
    locOut = fieldnames(sModOut);
    fldsOut = cell(numel(locOut), 1);
    for kk = 1 : numel(locOut(:))
        fldsOut{kk} = sModOut.(locOut{kk}).fields;
    end
    
      %NOT NECESSARY WITH PERSISTENT VARIABLE
%     %Use shortcuts for a few key variables to reduce evaluation time
%     indAll = strcmpi(locOut, 'all');
%     indGeo = strcmpi(fldsOut, 'geodetic');
%     indCasi = strcmpi(fldsOut, 'casi');
%     indFlow = strcmpi(fldsOut, 'flow');
%     
%     %Flow:
%     if ~isempty(indAll) && ~isempty(indFlow)
%         sModOut.all.flow(indTsPrintCurr,:,:) = sLand.flow;
% 
%         %Write to file
%         if isfield(sModOut.all, 'flow_path')
%             fileNm = fullfile(foldGrids, 'flow', [file_nm(sMeta.region, '_flow_', sModOut.all.date(indDate,:)) '.nc']);
%             print_grid_NC_v2(fileNm, sLand.flow, 'flow', sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
%         end 
%         
%         %Remove 
%         fldsOut{indAll}(indGeo) = [];
%     end
%     
%     %Casi:
%     if ~isempty(indAll) && ~isempty(indCasi) && ~isfield(sCryo, 'icdbr')
%         sModOut.all.casi(indTsPrintCurr,:,:) = 100*max(sCryo.scx, sCryo.icx);
% 
%         %Write to file
%         if isfield(sModOut.(locOut{kk}), 'casi_path')
%             fileNm = fullfile(foldGrids, 'casi', [file_nm(sMeta.region, 'casi', sModOut.all.date(indDate,:)) '.nc']);
%             print_grid_NC_v2(fileNm, squeeze(sModOut.all.casi(indTsPrintCurr,:,:)), 'casi', sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
%         end 
%         
%         %Remove 
%         fldsOut{indAll}(indCasi) = [];
%     end
%     
%     %Geodetic:
%     if ~isempty(indAll) && ~isempty(indGeo)
%         densBolch = find_att(sMeta.global, 'density_ice_Bolch', 'no_warning');
%         if isfield(sMeta, 'geodetic_density')
%             iceDens = sMeta.geodetic_density;
%         elseif ~isempty(densBolch) && ~isnan(densBolch)
%             iceDens = densBolch;
%         else
%             iceDens = find_att(sMeta.global, 'density_ice');
%         end
% 
%         geodeticTemp = (1000/iceDens)*(sCryo.sndwe + sCryo.icx.*sCryo.icdwe);
% 
%         sModOut.all.geodetic(indTsPrintCurr,:,:) = geodeticTemp;
% 
%         %Write to file
%         if isfield(sModOut.(locOut{kk}), 'geodetic_path')
%             fileNm = fullfile(foldGrids, 'geodetic', [file_nm(sMeta.region, 'geodetic', sModOut.all.date(indDate,:)) '.nc']);
%             print_grid_NC_v2(fileNm, geodeticTemp, 'geodetic', sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
%         end
% 
%         %Remove 
%         fldsOut{indAll}(indGeo) = [];
%     end
    
    %Populate output structure for each of the points:
    for kk = 1 : numel(locOut(:)) %Loop over grid pts to write to
        for ll = 1 : numel(fldsOut{kk}) %Loop over data fields for current grid point
            nmCurr = char(fldsOut{kk}{ll});
                        
            %Indice for 2D area array:
            indAreaCurr = sModOut.(locOut{kk}).indGage;
            
            %CALCULATE INDICES AND EXTRACT MODELED DATA
            if any(strcmpi(fldsLand, nmCurr)) %INFORMATION IN "sLand"
                %Define indices:
                nDimLand = ndims(sLand.(nmCurr));
                if nDimLand == 3 %3D array
                    if isfield(sLand, ['ind' nmCurr])
                        indInputTsCurr = sLand.(['ind' nmCurr]);
                    elseif strcmpi(nmCurr, 'pet')
                        indInputTsCurr = find(ismember(sLand.(['date' nmCurr])(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
                    else
                        indInputTsCurr = find(ismember(sLand.(['date' nmCurr]), sMeta.dateCurr, 'rows') == 1);
                    end
                    
                    if isempty(indInputTsCurr)
                        error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
                    end
                    szCurr = size(sLand.(nmCurr));
                    [gageRow, gageCol] = ind2sub(szMain, sModOut.(locOut{kk}).indGage);
                    indGageCurr = indInputTsCurr + (gageRow-1)*szCurr(1) + (gageCol-1)*szCurr(2)*szCurr(1);
                else %2D array
                    indGageCurr = sModOut.(locOut{kk}).indGage;
                end
                
                %Extract information:
                if strcmpi(locOut{kk}, 'all') 
%                     szCurr = size(sLand.(nmCurr));
%                     [gageRow, gageCol] = ind2sub(szMain, (1:szMain(1)*szMain(2)));
%                     indWrtCurr = indCurrTs + (gageRow-1)*szCurr(1) + (gageCol-1)*szCurr(2)*szCurr(1);

                    if nDimLand == 3
                        sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr,:,:) = squeeze(sLand.(nmCurr)(indInputTsCurr,:,:));
                    else
                        sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr,:,:) = sLand.(nmCurr);
                    end

                    %Write to file
                    if isfield(sModOut.(locOut{kk}), [char(nmCurr) '_path'])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(locOut{kk}).fields{ll}, sModOut.(locOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC_v2(fileNm, sLand.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
                    end
                else
                    sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr) = nansum(sLand.(nmCurr)(indGageCurr).*sHydro.area(indAreaCurr))/nansum(sHydro.area(indAreaCurr));
                end

            elseif any(strcmpi(fldsAtm, nmCurr)) %INFORMATION IN "sAtm"
                %Define indices:
                if ndims(sAtm.(nmCurr)) == 3 %3D array
                    if isfield(sAtm, ['ind' nmCurr])
                        indInputTsCurr = sAtm.(['ind' nmCurr]);
                    else
                        indInputTsCurr = find(ismember(sAtm.(['date' nmCurr]), sMeta.dateCurr, 'rows') == 1);
                    end
                    
                    if isempty(indInputTsCurr)
                        error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
                    end
                    szCurr = size(sAtm.(nmCurr));
                    [gageRow, gageCol] = ind2sub(szMain, sModOut.(locOut{kk}).indGage);
                    indGageCurr = indInputTsCurr + (gageRow-1)*szCurr(1) + (gageCol-1)*szCurr(2)*szCurr(1);
                    
                else %2D array
                    indGageCurr = sModOut.(locOut{kk}).indGage;
                end
                
                %Extract information:
                if strcmpi(locOut{kk}, 'all') 
                    sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr,:,:) = sAtm.(nmCurr);

                    %Write to file
                    if isfield(sModOut.(locOut{kk}), [char(nmCurr) '_path'])
                        fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(locOut{kk}).fields{ll}, sModOut.(locOut{kk}).date(indDate,:)) '.nc']);
                        print_grid_NC_v2(fileNm, sAtm.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
                    end
                else
                    sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr) = nansum(sAtm.(nmCurr)(indGageCurr).*sHydro.area(indAreaCurr))/nansum(sHydro.area(indAreaCurr));
                end
                
            elseif any(strcmpi(fldsCryo, nmCurr)) %INFORMATION IN "sCryo"
                %Potentially use ice grid (can be different from main
                %grid)
                %Also, write boolean parameter to record whether or not ice
                %present at current location
                if regexpbl(nmCurr, 'ic')
                    %Check and record whethere there is ice at requested point
                    indGageCurr = intersect(sModOut.(locOut{kk}).indGage, find(sCryo.icx ~= 0));
                    if ~isempty(indGageCurr)
                        sModOut.(locOut{kk}).icbl(indTsPrintCurr) = 1;
                    else
                        sModOut.(locOut{kk}).icbl(indTsPrintCurr) = 0;
                    end
                    
                    if iceGrid == 1 %Ice grid seperate and current field is ice variable
                        szGrid = szIce;
                        indCurr = sModOut.(locOut{kk}).indIce;
                    else
                        szGrid = szMain;
                        indCurr = sModOut.(locOut{kk}).indGage;
                    end
                else
                    szGrid = szMain;
                    indCurr = sModOut.(locOut{kk}).indGage;
                end
                    
                %'casi' is not actually a field; therefore, need to check
                %corresponding fields
                switch nmCurr
                    case 'casi'
                        fldCheck = 'scx';
                    case 'geodetic'
                        fldCheck = 'icdwe';
                    case 'stake' %fieldObsCurr = {'icdwe','sndwe'}
                        fldCheck = 'sndwe';
                    otherwise
                        fldCheck = nmCurr;
                end
                
                %Define indices (if ice, using special indices defined above):
                if ndims(sCryo.(fldCheck)) == 3 %3D array
                    if isfield(sCryo, ['ind' fldCheck])
                        indInputTsCurr = sCryo.(['ind' fldCheck]);
                    else
                        indInputTsCurr = find(ismember(sCryo.(['date' fldCheck]), sMeta.dateCurr, 'rows') == 1);
                    end
                    
                    if isempty(indInputTsCurr)
                        error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
                    end
                    szCurr = size(sCryo.(fldCheck));
                    [gageRow, gageCol] = ind2sub(szGrid, indCurr);
                    indGageCurr = indInputTsCurr + (gageRow-1)*szCurr(1) + (gageCol-1)*szCurr(2)*szCurr(1);
                    
                else %2D array
                    %Special case for ice grids:
                    if any(strcmpi(fldCheck, {'icdwe', 'iclr', 'icbl'})) && ~strcmpi(locOut{kk},'all')  
                        %Check and record whethere there is ice at requested point
                        indGageCurr = intersect(sModOut.(locOut{kk}).indGage, find(sCryo.icx ~= 0));
                        indAreaCurr = intersect(indAreaCurr, find(sCryo.icx ~= 0));
                        %Create boolean field to display presence of ice at
                        %gage location:
                        if ~isfield(sModOut.(locOut{kk}), 'icbl')
                            sModOut.(locOut{kk}).icbl = zeros(numel(sMeta.date(:,1)) - sMeta.indWrite + 1, 1);
                        end
                        %Record presence of ice at gage location:
                        if ~isempty(indGageCurr)
                            sModOut.(locOut{kk}).icbl(indTsPrintCurr) = 1;
                        else
                            sModOut.(locOut{kk}).icbl(indTsPrintCurr) = 0;
                        end
                    else
                        indGageCurr = indCurr;
                    end
                end
                  
                %Extract information:
                switch nmCurr 
                    case 'stake' %fieldObsCurr = {'icdwe','sndwe'}
                        if strcmpi(locOut{kk}, 'all')
                            warning('print_model_state:stakeAllGrid', 'The glacier stake observations are expected at a single point. An entire grid has not been programmed for.');
                        else
                            sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr) = nansum(...
                                (sCryo.('icdwe')(indGageCurr) + sCryo.('sndwe')(indGageCurr))...
                                .*sHydro.area(indAreaCurr))/nansum(sHydro.area(indAreaCurr));
                        end
                    case 'casi' %Special case for snow covered area:
                        if isfield(sCryo,'icdbr')
                            %If debris cover field, only count non-debris covered ice
                            if strcmpi(locOut{kk},'all') 
                                if isfield(sModOut.(locOut{kk}),['path' char(sModOut.(locOut{kk}).fields{ll})])
                                    %NOT PROGRAMMED YET
                                    if indTsPrintCurr == 1
                                        warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                                    end
                                else
                                    %NOT PROGRAMMED YET
                                    if indTsPrintCurr == 1
                                        warning('print_model_state:varNotProgrammed','Writing a gridded file for this type has not been programmed yet');
                                    end
                                end
                            else
                                casiTemp = max(sCryo.scx(indGageCurr), sCryo.icx(indGageCurr));
                                casiTemp(sCryo.icdbr > 0) = 0;

                                sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr) = ...
                                    100*nansum(casiTemp.*sHydro.area(indAreaCurr)) ...
                                /nansum(sHydro.area(indAreaCurr));
                            end
                        else %No debris cover field
                            if strcmpi(locOut{kk},'all') 
                                casiTemp = 100*max(sCryo.scx, sCryo.icx);
                                sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr,:,:) = casiTemp;

                                %Write to file
                                if isfield(sModOut.(locOut{kk}), [char(nmCurr) '_path'])
                                    fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(locOut{kk}).fields{ll}, sModOut.(locOut{kk}).date(indDate,:)) '.nc']);
                                    print_grid_NC_v2(fileNm, casiTemp, nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
                                end
                            else
                                sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr) = ...
                                    100*nansum(max(sCryo.scx(indGageCurr),sCryo.icx(indGageCurr)...
                                    ).*sHydro.area(indAreaCurr))/nansum(sHydro.area(indAreaCurr));
                            end
                        end
                    case 'geodetic'
                        iceDens = find_att(sMeta.global, 'density_ice_Bolch', 'no_warning');
                        if isfield(sMeta, 'geodetic_density')
                            iceDens = sMeta.geodetic_density;
                        elseif ~isempty(iceDens) && ~isnan(iceDens)
                            %Do nothing
                        else
                            iceDens = find_att(sMeta.global, 'density_ice');
                        end

                        if strcmpi(locOut{kk}, 'all')
                            geodeticTemp = (1000/iceDens)*(sCryo.sndwe + sCryo.icx.*sCryo.icdwe);

                            sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr,:,:) = geodeticTemp;

                            %Write to file
                            if isfield(sModOut.(locOut{kk}), [char(nmCurr) '_path'])
                                fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(locOut{kk}).fields{ll}, sModOut.(locOut{kk}).date(indDate,:)) '.nc']);
                                print_grid_NC_v2(fileNm, geodeticTemp, nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
                            end
                        else
                            geodeticTemp = (1000/iceDens)*(sCryo.sndwe(indGageCurr) + sCryo.icx(indGageCurr).*sCryo.icdwe(indGageCurr));
                            sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr) = nansum(geodeticTemp.*sHydro.area(indAreaCurr))/nansum(sHydro.area(indAreaCurr));
                        end
                    otherwise
                        if strcmpi(locOut{kk}, 'all')
                            sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr,:,:) = sCryo.(nmCurr);

                            %Write to file
                            if isfield(sModOut.(locOut{kk}), [char(nmCurr) '_path'])
                                fileNm = fullfile(foldGrids, nmCurr, [file_nm(sMeta.region, sModOut.(locOut{kk}).fields{ll}, sModOut.(locOut{kk}).date(indDate,:)) '.nc']);
                                print_grid_NC_v2(fileNm, sCryo.(nmCurr), nmCurr, sHydro.lon, sHydro.lat, sMeta.dateCurr, sMeta.dateCurr, 1);
                            end
                        else
                            sModOut.(locOut{kk}).(nmCurr)(indTsPrintCurr) = nansum(sCryo.(nmCurr)(indGageCurr).*sHydro.area(indAreaCurr))/nansum(sHydro.area(indAreaCurr));
                        end
                end
            else %FIELD NOT FOUND
                if indTsPrintCurr == 1
                    warning('backbone:unknownOutput',[char(39) ...
                        char(sModOut.(locOut{kk}).fields{ll}) char(39) ...
                        ' is not a known output type.']);
                end
            end
        end
        
        %Set values to nan at locations where dem is nan (values could be
        %wrong at these locations):
        if strcmpi(locOut{kk}, 'all')
            indNan2d = find(isnan(sHydro.dem));
            [nanRow, nanCol] = ind2sub(szMain, indNan2d);
            
            for ll = 1 : numel(fldsOut{kk})
                sz3d = size(sModOut.(locOut{kk}).(fldsOut{kk}{ll}));
                indNan3d = indTsPrintCurr + (nanRow-1)*sz3d(1) + (nanCol-1)*sz3d(2)*sz3d(1);
                sModOut.(locOut{kk}).(fldsOut{kk}{ll})(indNan3d) = nan;
            end
        end
    end
end

%Write output on select iterations to save time
if nargout > 0
    varargout{1} = sModOut;
end