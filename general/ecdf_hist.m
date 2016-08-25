function [n, c] = ecdf_hist( data , nBin)
%n = number of elements in each bin
%c = nanmean of all values in each bin

if iscell(data)
   error('e_cdf:cell','This function cannot handle cell input.'); 
end

[cdf, val] = e_cdf( data );
cdf = cdf(2:end);
val = val(2:end);

cdfCent = linspace(0,1,nBin+1);

c = nan(nBin,1);
n = nan(nBin,1);
for ii = 1 : nBin
    ind = find(cdf > cdfCent(ii) & cdf < cdfCent(ii+1));
    c(ii) = nanmean(val(c));
    n(ii) = numel(ind);
end