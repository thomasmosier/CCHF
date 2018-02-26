function trans = transmissivity_Coops(tasRng, elev, a, b, c)
%Coops, N. C., Waring, R. H., & Moncrieff, J. B. (2000). Estimating mean 
%monthly incident solar radiation on horizontal and inclined slopes from 
%mean monthly temperatures extremes. International Journal of 
%Biometeorology, 44(4), 204–211.


tClear = a*(0.65 + 0.008*elev);

%If elevation is 2d and tasrng is 3d, make 3d version of tClear
if ismatrix(tClear) && ndims(tasRng) == 3 && isequal(size(tClear), size(squeeze(tasRng(1,:,:))))
    temp = nan([1,size(tClear)], 'single');
    temp(1,:,:) = tClear;
    tClear = repmat(temp, numel(tasRng(:,1,1)), 1, 1);
elseif ndims(tasRng) == 3 && ndims(tClear) == 3 && ~isequal(size(tClear), size(tasRng))
    error('transmissivityCoops:arraySzDiff', 'Both input arrays are 3d, but their sizes are different.');
end

tDecay = b*(0.031+0.201*exp(-0.185*tasRng));

trans = tClear.*(1-exp(-tDecay.*tasRng.^c));

trans(trans > 1) = 1;
trans(trans < 0) = 0;
