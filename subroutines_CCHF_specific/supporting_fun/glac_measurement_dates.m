function [datesStart, datesEnd] = glac_measurement_dates(datesUse, datesFld, sMeta)

%Find mb dates (different from other fields because there's a start and end date):
indMbStrt = [];
indMbEnd = [];
unqYrs = unique(datesUse(:,1));
for zz = 1 : numel(unqYrs)
    indMbStrtTemp = find(ismember(datesUse, [unqYrs(zz)-1, sMeta.(datesFld)], 'rows')) + 1;
    indMbEndTemp  = find(ismember(datesUse, [unqYrs(zz)  , sMeta.(datesFld)(1,:)], 'rows'));

    %Colllect dates during years containing both a
    %start and an end for mb:
    if ~isempty(indMbStrtTemp) && ~isempty(indMbEndTemp)
        %Ensure mb dates do not include spin-up:
        if days_since(sMeta.dateCurr, datesUse(indMbStrtTemp,:), 'gregorian') >= 0
            indMbStrt = [indMbStrt; indMbStrtTemp];
            indMbEnd  = [ indMbEnd;  indMbEndTemp];
        end
    end
end

%Assign mb dates:
datesStart = datesUse(indMbStrt,:);
datesEnd = datesUse(indMbEnd,:);