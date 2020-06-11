function print_model_grid(sData, nmCurr, ptWrtCurr, indTsPrintCurr, sMeta, lon, lat, varargin)

global sModOut


%Indices that should be used (not set to nan):
if ~isempty(varargin(:))
    indInGrid2d = varargin{1};
else
    indInGrid2d = [];
end


%Reinitialize annual grid if it's the first day of the annual time-series
if strcmpi(ptWrtCurr,'annualgrid')
    indDateCurr = find(ismember(sModOut.(ptWrtCurr).([nmCurr '_dateStart']), sMeta.dateCurr, 'rows'));
    if any(indDateCurr)
%         disp(['resetting on ' num2str(sMeta.dateCurr(1)) '/' num2str(sMeta.dateCurr(2))]);
        nDays = days_since(sModOut.(ptWrtCurr).([nmCurr '_dateStart'])(indDateCurr,:), ...
            sModOut.(ptWrtCurr).([nmCurr '_dateEnd'])(indDateCurr,:), 'gregorian') + 1;
        sModOut.(ptWrtCurr).(nmCurr) = nan([nDays, numel(lat), numel(lon)], 'single');
        
        %Reset cntr:
        sModOut.(ptWrtCurr).([nmCurr '_cntr']) = 0;
    end
    
    %Advance counter:
    sModOut.(ptWrtCurr).([nmCurr '_cntr']) = sModOut.(ptWrtCurr).([nmCurr '_cntr']) + 1;
end


if isstruct(sData)
    if ndims(sData.(nmCurr)) == 3
        if isfield(sData, ['ind' nmCurr])
            indInputTsCurr = sData.(['ind' nmCurr]);
        elseif isfield(sData, ['date' nmCurr])
            indInputTsCurr = find(ismember(sData.(['date' nmCurr]), sMeta.dateCurr, 'rows') == 1);
        else
            error('printModelGrid:Array3dNoTime', ['A time field cannot be found for ' nmCurr '.']);
        end

        if ~isempty(indInGrid2d) %The selected indices to be used (not nan)
            if strcmpi(ptWrtCurr, 'all')
                sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)), indInGrid2d, indTsPrintCurr)) = ...
                    single(squeeze(sData.(nmCurr)(ind2_3d(size(sData.(nmCurr)), indInGrid2d, indInputTsCurr))));
            elseif strcmpi(ptWrtCurr,'annualgrid')
                sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)), indInGrid2d, sModOut.(ptWrtCurr).([nmCurr '_cntr']))) = ...
                    single(squeeze(sData.(nmCurr)(ind2_3d(size(sData.(nmCurr)), indInGrid2d, indInputTsCurr))));
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = ...
                    single(squeeze(sData.(nmCurr)(ind2_3d(size(sData.(nmCurr)), indInGrid2d, indInputTsCurr))));
            end
        else
            if strcmpi(ptWrtCurr, 'all')
                sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = single(squeeze(sData.(nmCurr)(indInputTsCurr,:,:)));
            elseif strcmpi(ptWrtCurr,'annualgrid')
                sModOut.(ptWrtCurr).(nmCurr)(sModOut.(ptWrtCurr).([nmCurr '_cntr']),:,:) = single(squeeze(sData.(nmCurr)(indInputTsCurr,:,:)));
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = single(squeeze(sData.(nmCurr)(indInputTsCurr,:,:)));
            end
        end
    else
        if ~isempty(indInGrid2d) %The selected indices to be used (not nan)
            if strcmpi(ptWrtCurr, 'all')
                sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)),indInGrid2d,indTsPrintCurr)) = single(sData.(nmCurr)(indInGrid2d));
            elseif strcmpi(ptWrtCurr,'annualgrid')
                sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)),indInGrid2d,sModOut.(ptWrtCurr).([nmCurr '_cntr']))) = single(sData.(nmCurr)(indInGrid2d));
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = single(sData.(nmCurr)(indInGrid2d));
            end
            
        else
            if strcmpi(ptWrtCurr, 'all')
                sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = single(sData.(nmCurr));
            elseif strcmpi(ptWrtCurr,'annualgrid')
                sModOut.(ptWrtCurr).(nmCurr)(sModOut.(ptWrtCurr).([nmCurr '_cntr']),:,:) = single(sData.(nmCurr));
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = single(sData.(nmCurr));
            end
        end
    end
elseif isnumeric(sData) && ismatrix(sData)
    if ~isempty(indInGrid2d) %The selected indices to be used (not nan)
        if strcmpi(ptWrtCurr, 'all')
            sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)),indInGrid2d,indTsPrintCurr)) = single(sData(indInGrid2d));
        elseif strcmpi(ptWrtCurr,'annualgrid')
        	sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)),indInGrid2d,sModOut.(ptWrtCurr).([nmCurr '_cntr']))) = single(sData(indInGrid2d));
        elseif strcmpi(ptWrtCurr, 'writegrid')
            gridTemp = single(sData(indInGrid2d));
        end
        
    else
        if strcmpi(ptWrtCurr, 'all')
            sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = single(sData);
        elseif strcmpi(ptWrtCurr,'annualgrid')
            sModOut.(ptWrtCurr).(nmCurr)(sModOut.(ptWrtCurr).([nmCurr '_cntr']),:,:) = single(sData);
        elseif strcmpi(ptWrtCurr, 'writegrid')
            gridTemp = single(sData);
        end
        
    end
else
    error('printModelGrid:numericArray3d','The input data is a numeric array with 3 dim. Numeric arrays must only have 2 dim.')
end


%Write to file
if isfield(sModOut.(ptWrtCurr), [char(nmCurr) '_path'])
    if strcmpi(ptWrtCurr,'annualgrid') 
        if any(ismember(sModOut.(ptWrtCurr).([nmCurr '_dateEnd']), sMeta.dateCurr, 'rows'))
            varUnit = {...
                'mrro', 'm/day'; ...
                'rnrf', 'm/day'; ...
                'snlr', 'm/day'; ...
                'iclr', 'm/day'; ...
                'snw', 'm/day'; ...
                'prsn', 'm/day'; ...
                'icdwe', 'm/day'; ...
                'icwe', 'm'};

%             if strcmpi(nmCurr, 'icdwe')
%                 keyboard
%             end
            grdSum = single(squeeze(sum(sModOut.(ptWrtCurr).(nmCurr), 1)));
            signTest = single(sign(grdSum));
            if numel(signTest == 1) > numel(signTest == -1)
                [grdMax, grdDOY] = max(single(sModOut.(ptWrtCurr).(nmCurr)), [], 1);
            else
                [grdMax, grdDOY] = min(single(sModOut.(ptWrtCurr).(nmCurr)), [], 1);
            end
                grdDOY = squeeze(grdDOY);
                grdMax = squeeze(grdMax);
                grdDOY(isnan(grdSum)) = nan;
                grdMax(isnan(grdSum)) = nan;
            
%             %This section of code could be used to find the max negative or
%             %position in the same 3d array
%             temp = sModOut.(ptWrtCurr).(nmCurr);
%             [~,index] = max(abs(temp), [], 1);
%                 index = squeeze(index);
%                 [i,j] = meshgrid((1:size(temp,3)), (1:size(temp,2)));
%                 out = temp(sub2ind(size(temp), index, j, i));
%                 [grdMax(:),out(:)]


            %Convert from seconds to years:
            indUnit = find(ismember(varUnit(:,1), nmCurr) == 1);
            if isempty(indUnit)
                error('printModelState:unitsNotfound',['Units for ' nmCurr ' not found. Add to list in this subroutine.']);
            end

            if strcmpi(varUnit{indUnit,2}, 'm/s')
                secToDy = 86400;
%                 secToYr = (indWYEnd(ii) - indWYBeg(ii) + 1) * secToDy;
                grdMax = grdMax * secToDy;
                grdSum = grdSum * secToDy;

                unitMax = 'm/day';
                unitSum = 'm/yr';
            elseif strcmpi(varUnit{indUnit,2}, 'm/day')
                unitMax = 'm/day';
                unitSum = 'm/yr';
            elseif strcmpi(varUnit{indUnit,2}, 'm')
                unitMax = 'm';
                unitSum = 'm';
            else
                error('hydroStats:unknownUnits', ['Units ' varUnit{indUnit,2} ' for variable ' nmCurr ' are not recognize.']); 
            end

            
            %Create date string for output file
            indDateCurr = find(ismember(sModOut.(ptWrtCurr).([nmCurr '_dateEnd']), sMeta.dateCurr, 'rows'));
            dateBnds = [sModOut.(ptWrtCurr).([nmCurr '_dateStart'])(indDateCurr,:); sModOut.(ptWrtCurr).([nmCurr '_dateEnd'])(indDateCurr,:)];

            if dateBnds(1,2) < 10
                strMnthStrt = ['0' num2str(dateBnds(1,2))];
            else
                strMnthStrt = num2str(dateBnds(1,2));
            end

            if dateBnds(2,2) < 10
                strMnthEnd = ['0' num2str(dateBnds(2,2))];
            else
                strMnthEnd = num2str(dateBnds(2,2));
            end

            strDateWrt = [num2str(dateBnds(1,1)) strMnthStrt '01' '-' ...
                num2str(dateBnds(2,1)) strMnthEnd num2str(eomday(dateBnds(2,1), dateBnds(2,2)))];

            foldWrt = fullfile(sMeta.foldWrtData, nmCurr);
            pathOutSum = fullfile(foldWrt, [nmCurr '_total_' strDateWrt]);
            pathOutMax = fullfile(foldWrt, [nmCurr '_max_' strDateWrt]);
            pathOutDOY = fullfile(foldWrt, [nmCurr '_dayofyr-max_' strDateWrt]);
            
            if regexpbl(sMeta.wrtTyp, {'nc', 'netcdf'})
                warning('off', 'printGridNc:nanTime');
                print_grid_NC_v2([pathOutSum, '.nc'], grdSum, nmCurr, lon, lat, sMeta.dateCurr, dateBnds, unitSum, 2);
                print_grid_NC_v2([pathOutMax, '.nc'], grdMax, nmCurr, lon, lat, sMeta.dateCurr, dateBnds, unitMax, 3);
                print_grid_NC_v2([pathOutDOY, '.nc'], grdDOY, nmCurr, lon, lat, sMeta.dateCurr, dateBnds,   'day', 0);
                warning('on', 'printGridNc:nanTime');
            elseif regexpbl(sMeta.wrtTyp, 'asc')
                hdrWrt = ESRI_hdr(lon, lat, 'corner');
                write_ESRI_v4(grdSum, hdrWrt, [pathOutSum, '.asc'], 2);
                write_ESRI_v4(grdMax, hdrWrt, [pathOutMax, '.asc'], 3);
                write_ESRI_v4(grdDOY, hdrWrt, [pathOutDOY, '.asc'], 0);
            end
        end
    else
        %Create subdirectories for year (if only one year present:
        yrsUniq = unique(sMeta.dateCurr(:,1));
        if numel(yrsUniq) == 1
            foldWrt = fullfile(sMeta.foldWrtData, nmCurr, num2str(yrsUniq));
        else
            foldWrt = fullfile(sMeta.foldWrtData, nmCurr);
        end

        %Create full path:
        fileWrt = [char(file_nm(sMeta.region{sMeta.siteCurr}, nmCurr, sMeta.dateCurr)) '.nc'];
        pathWrt = fullfile(foldWrt, fileWrt);

        %Write file
        if strcmpi(ptWrtCurr, 'all')
            print_grid_NC_v2(pathWrt, squeeze(sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:)), nmCurr, lon, lat, sMeta.dateCurr, sMeta.dateCurr, 1);
        elseif strcmpi(ptWrtCurr, 'writegrid')
            print_grid_NC_v2(pathWrt, gridTemp, nmCurr, lon, lat, sMeta.dateCurr, sMeta.dateCurr, 1);
        else
            error('printModelGrid:ptTypeUnknown',['The point type ' ptWrtCurr ' has not been programmed for.']);
        end
    end
end