function [dateStart, strMult]= date_min(dateIn)

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
            dateStart = dateAll(1,:);
        else
            dateStart = min(dateAll);
            strMult = 'mult';
        end
    else
        dateStart = dateIn{1};
    end
else
    dateStart = dateIn(1);
end