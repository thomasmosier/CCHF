function chromosomes = chrom_mutate(chromosomes, sOpt)


nMutePts = find_att(sOpt.att,'mutation sites');
perMut = find_att(sOpt.att,'mutation rate');

mutTemp = rand( numel(chromosomes), 1);
[~, indMut] = sort(mutTemp);
indMut = indMut(1:ceil(perMut*numel(chromosomes)));

chromLeng = numel(chromosomes{1});
for ii = 1 : numel(indMut)
    mutSite = randi([1, chromLeng], nMutePts, 1);

    chromosomes{indMut(ii)}(mutSite) = num2str(1 - mod(chromosomes{indMut(ii)}(mutSite),2));
end
