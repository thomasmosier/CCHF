function d = haversine_neighbors(lon,lat,r,varargin)


if any(size(lon) == 1) && any(size(lon) == 1)
    [lon, lat] = meshgrid(lon,lat);
elseif any(size(lon) == 1) && ~any(size(lon) == 1) || ~any(size(lon) == 1) && any(size(lon) == 1)
    error('haversine:matrixAndVec','One input appears to be a vector and the other a matrix.  Both inputs must have same form.')
end

[ic, id] = ixneighbors(lon);

d = 2*r*asin(sqrt(sind(0.5*(lat(ic)-lat(id))).^2 + cosd(lat(ic)).*cosd(lat(id)).*sind(0.5*(lon(id)-lon(ic))).^2));

if ~isempty(varargin(:)) 
    if regexpbl(varargin{1},'sparse')
        d = sparse(ic, id, d, numel(lon), numel(lon));
    elseif regexpbl(varargin{1},'full')
        d = full(sparse(ic, id, d, numel(lon), numel(lon)));
    else
        error('haversine:outputType',['Output format ' varargin{1} ' is not recognized.'])
    end
end