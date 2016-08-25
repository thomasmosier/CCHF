function [pJoint, xC, yC] = joint_hist(x,y,varargin)


if numel(varargin) == 1
    nBin = varargin{1};
    
    [~, valX] = ecdf_data(x, nBin);
    [~, valY] = ecdf_data(y, nBin);
elseif numel(varargin) == 2
    valX =  varargin{1};
    valY =  varargin{2};
    nBin = numel(valX);
else
    error('joint_hist:arguments','The number of inputs is not acceptable.')
end

if valX(1) == valX(2) && valY(1) == valY(2)
    valX = valX(2:end);
    valY = valY(2:end);
    nBin = nBin - 1;
end

pJoint = zeros(nBin);
xC = nan(nBin);
yC = nan(nBin);

for ii = 1 : nBin %Loop over x
    if ii == 1
        indX = find(x <= valX(ii));
    else
        indX = find(x > valX(ii-1) & x <= valX(ii));
    end
    
    for jj = 1 : nBin %Loop over y
        if jj == 1
            indY = find(y(indX) <= valY(jj));
        else
            indY = find(y(indX) > valY(jj-1) & y(indX) <= valY(jj));
        end
        
        indJ = indX(indY);
        if ~isempty(indJ)
            pJoint(ii,jj) = numel(indJ);
            xC(ii,jj) = nanmean(x(indJ));
            yC(ii,jj) = nanmean(y(indJ));
        end
    end
end

%Normalize nJoint:
pJoint = 100*pJoint/sum2d(pJoint);

