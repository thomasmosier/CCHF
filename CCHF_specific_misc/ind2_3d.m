function ind3 = ind2_3d(sz3,ind2,tSlice)
%Given indices for 2d array and a 3d array where the 2d array is the 2nd
%and 3rd dimension, create the indices for the 3d array.

if numel(ind2 > 1) && numel(ind2(1,:)) > numel(ind2(:,1))
   ind2 = ind2'; 
end

[r, c] = ind2sub(sz3(2:3),ind2);

ind3 = tSlice*ones(size(r)) + (r-1)*sz3(1) + (c-1)*sz3(2)*sz3(1);