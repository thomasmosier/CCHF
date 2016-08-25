function out = stream(fac,thresh)

out = nan(size(fac));

out(fac > thresh) = 1;