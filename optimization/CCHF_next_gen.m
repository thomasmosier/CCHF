function coef = CCHF_next_gen(coef, fitScore, ii, prmBnds, sOpt)

%Indices of coef are: [iteration, member, parameter]

% %%PENALIZE MODELLED FLOWS CONTAINING UNACCEPTABLE RESULTS:
% %Defines penalty for flowrates resulting in ?NaN? values.
% penNaN = -99999999;
% 
% %If flowrate is 'nan', changes to very large negative number:
% gageMod( isnan(gageMod) ) = penNaN;
% %Coefficients incur penalty if more than half of flow values are zero
% if sum(gageMod == 0) > 0.5*length(gageMod)
%     gageMod(gageMod == 0) = penNaN;
% end

 
if regexpbl(sOpt.type,{'GA','genetic','evolution'}) && ~regexpbl(sOpt.type,{'binary'}) && regexpbl(sOpt.type,{'real'})
    %Reduce dimensionality of inputs:
    coef = squeeze(coef(ii,:,:));
    fitScore = squeeze(fitScore(ii,:));
    
    if ~isempty(find_att(sOpt.att,'fitness requirement'))
        coef = GA_roulette(coef, fitScore, sOpt);
    end
    coef = coef_randomize(coef, fitScore, sOpt);
    coef = GA_cross(coef, sOpt);
    coef = GA_mutate(coef, prmBnds, sOpt);
elseif regexpbl(sOpt.type,{'GA','genetic','evolution'}) && regexpbl(sOpt.type,{'binary'})
    %Reduce dimensionality of inputs:
    coef = squeeze(coef(ii,:,:));
    fitScore = squeeze(fitScore(ii,:));
    
    %%CONVERT COEFFICIENTS TO BINARY CHROMOSOMES AND MUTATE
    %Find desirable range of parameter values
    paramLeng = chrom_scaling(prmBnds);
    %Convert coefficients to binary chromosomes
    chromosomes = chrom_dec2bin(coef,prmBnds,paramLeng);
    %Select the chromosomes to "reproduce"
    chromosomes = chrom_roulette(chromosomes, fitScore, sOpt);
    %Perform crossover on selected chromosomes
    chromosomes = chrom_cross(chromosomes, sOpt);
    %Randomly mutate set percentage of chromosomes at determined number of
    %points
    chromosomes = chrom_mutate(chromosomes, sOpt);
    %Convert from binary chromosomes back to coefficients:
    coef = chrom_bin2dec(chromosomes, prmBnds, paramLeng); 
elseif regexpbl(sOpt.type,{'hybrid'})
    coef = opt_hybrid(coef, fitScore, prmBnds, ii);
elseif regexpbl(sOpt.type,{'PSO','particle swarm'})
    coef = PSO(coef, fitScore, prmBnds, ii);
elseif regexpbl(sOpt.type,{'Uniform','sampling'},'and')
    %Do nothing because all parameter sets have already been defined
elseif regexpbl(sOpt.type,{'Monte','Carlo'},'and')
    coef = Monte_Carlo(prmBnds,numel(fitScore));
else
    error('CCHF_next_gen:unknownType',['The optimization scheme ' ...
        char(39) sOpt.type char(39) ' has not been coded for.'])
end

%Ensure no coefficients are out of bounds:
for ii = 1 : numel(coef(1,:))
    coef(coef(:,ii) < prmBnds(ii,1), ii) = prmBnds(ii, 1);
    coef(coef(:,ii) > prmBnds(ii,2), ii) = prmBnds(ii, 2);
end