function sModOut = print_model_grid(sModOut, sData, nmCurr, ptWrtCurr, indTsPrintCurr, sMeta, lon, lat, varargin)

if isstruct(sData)
    if ndims(sData.(nmFld)) == 3

        if isfield(sData, ['ind' nmCurr])
            indInputTsCurr = sData.(['ind' nmCurr]);
        elseif isfield(sData, ['date' nmCurr])
            indInputTsCurr = find(ismember(sData.(['date' nmCurr]), sMeta.dateCurr, 'rows') == 1);
        else
            error('printModelGrid:Array3dNoTime', ['A time field cannot be found for ' nmCurr '.']);
        end

        sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = squeeze(sData.(nmCurr)(indInputTsCurr,:,:));
    else
        sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = sData.(nmCurr);
    end
elseif isnumeric(sData) && ismatrix(sData)
    sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:) = sData;
else
    error('printModelGrid:numericArray3d','The input data is a numeric array with 3 dim. Numeric arrays must only have 2 dim.')
end

if ~isempty(varargin(:))
    indNan2d = varargin{1};
    %Set indices outside domain to nan:
    sModOut.(ptWrtCurr).(nmCurr) = array_time_slice_nan(sModOut.(ptWrtCurr).(nmCurr), indTsPrintCurr, indNan2d);
end

%Write to file
if isfield(sModOut.(ptWrtCurr), [char(nmCurr) '_path'])
    fileNm = fullfile(sMeta.foldWrtData, nmCurr, [char(file_nm(sMeta.region{sMeta.siteCurr}, nmCurr, sMeta.dateCurr)) '.nc']);
    print_grid_NC_v2(fileNm, squeeze(sModOut.(ptWrtCurr).(nmCurr)(indTsPrintCurr,:,:)), nmCurr, lon, lat, sMeta.dateCurr, sMeta.dateCurr, 1);
end