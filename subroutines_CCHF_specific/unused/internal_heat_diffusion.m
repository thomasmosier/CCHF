function tmpArray = internal_heat_diffusion(tmpArray, thick2d, nLyr2d, cf)
keyboard
%tmp3d = [indices (corresponding to land surface), layers]
%layer = 1 = surface; layer N = bottom
for ii = 1 : numel(nLyr2d)
    tmpArray(ii,1:nLyr2d(ii)) = heat_diffusion_CN(tmpArray(ii,1:nLyr2d(ii)), tmpArray(ii,1), tmpArray(ii,nLyr2d(ii)), cf/(thick2d(ii)^2));
end