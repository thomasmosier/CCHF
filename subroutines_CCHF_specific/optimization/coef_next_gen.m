% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


function coef = coef_next_gen(coef, fitScore, ii, prmBnds, sOpt)

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
elseif strcmpi(sOpt.type, 'pso') || regexpbl(sOpt.type, 'particle swarm')
    coef = pso(coef, fitScore, prmBnds, ii);
elseif strcmpi(sOpt.type, 'apso')
    coef = apso(coef, fitScore, prmBnds, ii);
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