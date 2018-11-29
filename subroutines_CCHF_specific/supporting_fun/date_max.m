function [dateEnd, strMult] = date_max(dateIn)

strMult = '';
if iscell(dateIn)
    if numel(dateIn(:)) > 1
        dateAll = nan(numel(dateIn(:)), numel(dateIn{1}(:)));
        for ii = 1 : numel(dateAll(:,1))
            dateAll(ii,:) = dateIn{ii}(:);
        end
        blEqual = nan(1, numel(dateAll(1,:)));
        for ii = 1 : numel(blEqual)
            if all(dateAll(1,ii) == dateAll(:,ii))
                blEqual(ii) = 1;
            else
                blEqual(ii) = 0;
            end
        end
        if all(blEqual == 1)
            dateEnd = dateAll(1,:);
        else
            dateEnd = max(dateAll);
            strMult = 'mult';
        end
    else
        dateEnd = dateIn{1};
    end
else
    dateEnd = dateIn(1);
end