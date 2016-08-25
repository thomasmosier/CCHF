function chrom = chrom_dec2bin(parameters, prmBnds, paramLeng)

chrom = cell(length(parameters(:,1)),1);
[ chrom{:} ]= deal( blanks( sum(paramLeng)));

if iscell(prmBnds)
    prmBnds = cell2mat(prmBnds);
end

nPts = 2.^paramLeng - 1;

nBit = 1;
for ii = 1 : numel(paramLeng)  %Loop over each parameter
    
    if prmBnds(ii,2) - prmBnds(ii,1) == 0
        continue;
    else
        
        for jj = 1 : length( parameters(:,1) )  %Loop over each member of parameter population
            chrom{jj}(nBit:nBit + paramLeng(ii) - 1) = dec2bin((parameters(jj,ii) - prmBnds(ii,1))*nPts(ii)/(prmBnds(ii,2) - prmBnds(ii,1)),paramLeng(ii));
        end
        nBit = nBit + paramLeng(ii);
    end
end