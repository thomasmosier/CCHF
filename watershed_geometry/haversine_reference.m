function d = haversine_reference(lon0,lat0,lonVec,latVec,r)

szGrid = size(lonVec);
    
if ~any(size(lonVec) == 1) && ~any(size(latVec) == 1) && isequal(size(latVec),size(lonVec))
    lonVec = reshape(lonVec,[],1);
    latVec = reshape(latVec,[],1);
elseif any(size(lonVec) == 1) && ~any(size(latVec) == 1) || ~any(size(lonVec) == 1) && any(size(latVec) == 1)
    error('haversine:matrixAndVec','One input appears to be a vector and the other a matrix.  Both inputs must have same form.')
end

lon0 = repmat(lon0,numel(lonVec),1);
lat0 = repmat(lat0,numel(latVec),1);

d = 2*r*asin(sqrt(sind(0.5*(lat0-latVec)).^2 + cosd(lat0).*cosd(latVec).*sind(0.5*(lonVec-lon0)).^2));

if all(szGrid ~= 1)
   d = reshape(d,szGrid); 
end
