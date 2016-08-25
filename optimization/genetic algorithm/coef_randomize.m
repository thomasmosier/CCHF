function coef = coef_randomize(coef, score, sOpt)
%Puts elite performing coefficient sets first and then randomized the rest:

nElite = find_att(sOpt.att,'elite pop','no_warning');

if isempty(nElite)
   nElite = 0; 
end

[~, coefOrd] = sort(rand( numel(coef(:,1)), 1));

for ii = 1 : nElite
    coefOrd(coefOrd == ii) = [];
end

if nElite > 0
    [~, sortOrd] = sort(score);
    coef = coef(sortOrd,:);
    
    if nElite < numel(coef(:,1))
        coef = [coef(1:nElite,:) ; coef(coefOrd,:)];
    end
else
    coef = coef(coefOrd,:);
end