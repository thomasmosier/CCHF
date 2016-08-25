function slope = slope_Liston(dem,fdr,lon,lat)
%Slope is the angle of inclination along the flowpath


slope = nan(size(dem));

if sum(size(lon) == 1) > 0 && sum(size(lon) > 1) > 0
    [lon, lat] = meshgrid(lon, lat);
end

tic
for jj = 2 : numel(dem(:,1)) - 1
   
%    disp([num2str(jj) ', ' num2str(toc)]);
   
   for ii = 2 : numel(dem(1,:)) - 1
       crd = fdr_inv_mult(fdr(jj,ii),2);
       
        if numel(crd) == 2
            slope(jj,ii) = atand(slope_2pts(dem(jj,ii) - dem(jj+crd(1),ii+crd(2)), [lon(jj,ii), lat(jj,ii)], [lon(jj+crd(1),ii+crd(2)), lat(jj+crd(1),ii+crd(2))]));
        elseif numel(crd) == 4
            if sum(abs(crd(1,:) - crd(2,:)) > 1)
                slope(jj,ii) = 0;
            else
                slope(jj,ii) = mean([ ...
                    atand(slope_2pts(dem(jj,ii) - dem(jj+crd(1,1),ii+crd(1,2)), [lon(jj,ii), lat(jj,ii)], [lon(jj+crd(1,1),ii+crd(1,2)), lat(jj+crd(1,1),ii+crd(1,2))])), ...
                    atand(slope_2pts(dem(jj,ii) - dem(jj+crd(2,1),ii+crd(2,2)), [lon(jj,ii), lat(jj,ii)], [lon(jj+crd(2,1),ii+crd(2,2)), lat(jj+crd(2,1),ii+crd(2,2))])) ...
                    ]);
            end
        else
            slope(jj,ii) = 0;
        end
           
   end
end