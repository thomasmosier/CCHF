function parameters = chrom_bin2dec(chromosomes, prmBnds, paramLeng)

if iscell(prmBnds)
    prmBnds = cell2mat(prmBnds);
end

%1st indice is parameter set (i.e. individual), 2nd indice is parameter
parameters = zeros( numel(chromosomes(:,1)), length(prmBnds(:,1)));

nPts = 2.^paramLeng - 1;

ind = 1;
for ii = 1 : length(prmBnds(:,1))
    if isequal(prmBnds(ii,1), prmBnds(ii,2))
       parameters(:,ii) = prmBnds(ii,1);
    else
        for jj = 1 : numel(chromosomes(:))
            parameters(jj,ii) = prmBnds(ii,1) ...
                + (prmBnds(ii,2) - prmBnds(ii,1))/nPts(ii)*bin2dec( chromosomes{jj}(ind : ind + paramLeng(ii)-1) ); 
        end
        ind = ind + paramLeng(ii); 
    end
    
end

