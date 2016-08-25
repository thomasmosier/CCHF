function [sinkBl, dem] = is_sink(dem)


dimDem = size(dem);
sinkBl = nan(dimDem); 

[~, indDem] = sort(dem(:),'ascend');
%Remove all nan elements:
indDem(isnan(dem(indDem))) = [];

for ii = 1 : numel(indDem) %Loop over columns
    [r,c] = ind2sub(dimDem,indDem(ii));
    
    if r == 1 || r == dimDem(1) || c == 1 || c == dimDem(2)
       continue
    end
    
    demCurr = reshape(dem(r-1:r+1, c-1:c+1),[],1);
    [val,~] = min(demCurr);

    ind = find(demCurr == val);
    if numel(ind) == 1 && ind == 5 %If center is lowest, it may be a sink
        sinkBl(r,c) = 1;
        demCurr(5) = [];
        dem(r,c) = min(demCurr);
    end
end