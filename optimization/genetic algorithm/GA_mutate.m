function coef = GA_mutate(coef, prmBnds, sOpt)


nMutePts = find_att(sOpt.att,'mutation sites');
nMemMut = find_att(sOpt.att,'mutation pop');
nElite = find_att(sOpt.att,'elite pop','no_warning');

if isempty(nElite)
   nElite = 0; 
end

% [~, coefMut] = sort(rand( numel(coef(:,1)), 1));
% coefMut = coefMut(1:nMemMut);

for ii = numel(coef(:,1)) : -1 : numel(coef(:,1)) - nMemMut + 1
    mutCoef = randi([1, numel(coef(1,:))], nMutePts, 1);

    for jj = 1 : nMutePts
        coef(ii, mutCoef(jj)) = rand(1,1)*(prmBnds(mutCoef(jj),2)-prmBnds(mutCoef(jj),1)) + prmBnds(mutCoef(jj),1);
    end
end
