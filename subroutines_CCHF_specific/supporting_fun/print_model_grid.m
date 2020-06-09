function print_model_grid(sData, nmCurr, ptWrtCurr, indTsPrintCurr, sMeta, lon, lat, varargin)

global sModOut


%Indices that should be used (not set to nan):
if ~isempty(varargin(:))
    indInGrid2d = varargin{1};
else
    indInGrid2d = [];
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
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = ...
                    single(squeeze(sData.(nmCurr)(ind2_3d(size(sData.(nmCurr)), indInGrid2d, indInputTsCurr))));
            end
        else
            if strcmpi(ptWrtCurr, 'all')
                sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = single(squeeze(sData.(nmCurr)(indInputTsCurr,:,:)));
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = single(squeeze(sData.(nmCurr)(indInputTsCurr,:,:)));
            end
        end
    else
        if ~isempty(indInGrid2d) %The selected indices to be used (not nan)
            if strcmpi(ptWrtCurr, 'all')
                sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)),indInGrid2d,indTsPrintCurr)) = single(sData.(nmCurr)(indInGrid2d));
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = single(sData.(nmCurr)(indInGrid2d));
            end
            
        else
            if strcmpi(ptWrtCurr, 'all')
                sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = single(sData.(nmCurr));
            elseif strcmpi(ptWrtCurr, 'writegrid')
                gridTemp = single(sData.(nmCurr));
            end
        end
    end
elseif isnumeric(sData) && ismatrix(sData)
    if ~isempty(indInGrid2d) %The selected indices to be used (not nan)
        if strcmpi(ptWrtCurr, 'all')
            sModOut.(ptWrtCurr).(nmCurr)(ind2_3d(size(sModOut.(ptWrtCurr).(nmCurr)),indInGrid2d,indTsPrintCurr)) = single(sData(indInGrid2d));
        elseif strcmpi(ptWrtCurr, 'writegrid')
            gridTemp = single(sData(indInGrid2d));
        end
        
    else
        if strcmpi(ptWrtCurr, 'all')
            sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = single(sData);
        elseif strcmpi(ptWrtCurr, 'writegrid')
            gridTemp = single(sData);
        end
        
    end
else
    error('printModelGrid:numericArray3d','The input data is a numeric array with 3 dim. Numeric arrays must only have 2 dim.')
end


%Write to file
if isfield(sModOut.(ptWrtCurr), [char(nmCurr) '_path'])
    %Crete subdirectories for year (if only one year present:
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
    else
        print_grid_NC_v2(pathWrt, gridTemp, nmCurr, lon, lat, sMeta.dateCurr, sMeta.dateCurr, 1);
    end
end