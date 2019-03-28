function sModOut = print_model_grid(sModOut, sData, nmCurr, ptWrtCurr, indTsPrintCurr, sMeta, lon, lat, varargin)

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

        if ~isempty(indInGrid2d)
            szCurr = size(sData.(nmCurr));
            tempOut = nan(szCurr(2:3), 'single');
            tempData = squeeze(sData.(nmCurr)(indInputTsCurr,:,:));
            tempOut(indInGrid2d) = tempData(indInGrid2d);
            sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = tempOut;
        else
            sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = squeeze(sData.(nmCurr)(indInputTsCurr,:,:));
        end
    else
        if ~isempty(indInGrid2d)
            tempOut = nan(size(sData.(nmCurr)), 'single');
            tempOut(indInGrid2d) = sData.(nmCurr)(indInGrid2d);
            sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = tempOut;
        else
            sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = sData.(nmCurr);
        end
    end
elseif isnumeric(sData) && ismatrix(sData)
    if ~isempty(indInGrid2d)
        tempOut = nan(size(sData.(nmCurr)), 'single');
        tempOut(indInGrid2d) = sData(indInGrid2d);
        sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = tempOut;
    else
        sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = sData;
    end
else
    error('printModelGrid:numericArray3d','The input data is a numeric array with 3 dim. Numeric arrays must only have 2 dim.')
end


%Write to file
if isfield(sModOut.(ptWrtCurr), [char(nmCurr) '_path'])
    fileNm = fullfile(sMeta.foldWrtData, nmCurr, [char(file_nm(sMeta.region{sMeta.siteCurr}, nmCurr, sMeta.dateCurr)) '.nc']);
    print_grid_NC_v2(fileNm, squeeze(sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:)), nmCurr, lon, lat, sMeta.dateCurr, sMeta.dateCurr, 1);
end