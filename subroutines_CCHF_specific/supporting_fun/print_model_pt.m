function print_model_pt(sData, nmCurr, ptWrtCurr, indTsOutCurr, area, indGage, indArea, sMeta, varargin)

global sModOut

%Set indices outside domain to nan:
if ~isempty(varargin(:))
    indInModel2d = varargin{1};
    indArea = intersect(indArea(:), indInModel2d(:));
    indGage = intersect(indGage(:), indInModel2d(:));
end


if isstruct(sData)
    szInput = size(sData.(nmCurr));
elseif isnumeric(sData)
    szInput = size(nmCurr);
else
    error('printModelPt:unknownType',['The input is of type ' class(sData) ', which has not been programmed for. Struct or numeric is expected.']);
end


if numel(szInput) == 3 %3D array
    if isstruct(sData)
        if isfield(sData, ['ind' nmCurr])
            indTsInCurr = sData.(['ind' nmCurr]);
        elseif isfield(sData, ['date' nmCurr])
            if strcmpi(nmCurr, 'pet')
                indTsInCurr = find(ismember(sData.(['date' nmCurr])(:,2:end), sMeta.dateCurr(:,2:end), 'rows') == 1);
            else   
                indTsInCurr = find(ismember(sData.(['date' nmCurr]), sMeta.dateCurr, 'rows') == 1);
            end
        end

        if isempty(indTsInCurr)
            error('printModelPt:noTsPresent',['No time indice is present in field ' ...
                nmCurr ', but the input array has 3 dimensions.']);
        end
        szCurr = size(sData.(nmCurr));
        [gageRow, gageCol] = ind2sub(size(area), indGage(:));
        
        if numel(indTsInCurr) == 1
            indTsInCurr = repmat(indTsInCurr, numel(gageRow), 1);
        else
            error('printModelPt:multTimeInd',['There are ' num2str(indTsInCurr) ' input time indices. Only one is allowed.']);
        end
        indGage = sub2ind(szCurr, indTsInCurr, gageRow, gageCol);
    else
        error('printModelPt:inputNotStruct3D',['The input array is of class ' ...
            class(sData) ', but struct is required for 3D data.']);
    end
end  


if numel(indGage) == numel(indArea) && ~isempty(indGage)
    if isstruct(sData)
        sModOut.(ptWrtCurr).(nmCurr)(indTsOutCurr) = nansum(sData.(nmCurr)(indGage(:)).*area(indArea(:)))/nansum(area(indArea(:)));
    elseif isnumeric(sData)
        sModOut.(ptWrtCurr).(nmCurr)(indTsOutCurr) = nansum(sData(indGage(:)).*area(indArea(:)))/nansum(area(indArea(:)));
    end
elseif numel(indGage) == numel(indArea) && numel(indGage) == 0
    if isstruct(sData)
        sModOut.(ptWrtCurr).(nmCurr)(indTsOutCurr) = nan;
    elseif isnumeric(sData)
        sModOut.(ptWrtCurr).(nmCurr)(indTsOutCurr) = nan;
    end
else
    error('printModelPt:indicesDiffSz', ['The gage indice has ' ...
        num2str(numel(indGage)) ' elements and the area indice has ' ...
        num2str(numel(indArea)) ' elements.'])
end