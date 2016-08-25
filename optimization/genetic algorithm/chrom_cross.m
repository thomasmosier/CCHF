function chromosomes = chrom_cross(chromosomes, sOpt)


nCrossPts = find_att(sOpt.att,'crossover sites');
perCross = find_att(sOpt.att,'crossover rate');

%Randomly select indices for chromosomes to cross
xTemp = rand(numel(chromosomes), 1);

%Using indices of the sort function ensures no repreated indices
[~, indCross] = sort(xTemp);
%Keep only the percent that are supposed to be crossed:
indCross = indCross(1:ceil(perCross*numel(chromosomes)));
%Ensure pairs:
if mod(numel(indCross),2) == 1    
    indCross(end) = []; 
end


%Loop over chromosomes, crossing each consecutive pair:
for ii = 1 : length(indCross)/2
    %Select sites where crossing occurs:
    indsites = randi([1, numel(chromosomes{1})], nCrossPts, 1);
    %Perform crossover for current pair:
    for jj = 1 : nCrossPts
        chromosomes{indCross(2*(ii-1)+1), : } = [chromosomes{indCross(2*(ii-1)+1)}(1:indsites(jj)), chromosomes{indCross(2*(ii-1)+2)}(indsites(jj)+1:end)];
        chromosomes{indCross(2*(ii-1)+2), : } = [chromosomes{indCross(2*(ii-1)+2)}(1:indsites(jj)), chromosomes{indCross(2*(ii-1)+1)}(indsites(jj)+1:end)];
    end
end