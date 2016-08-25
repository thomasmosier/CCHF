function upslope = upstream(fdr, lon, lat, pt)

keyboard

%fdr = sparse matrix where row is flowFrom cell and column is the cell that
%flow goes to
%pt = [lat, lon]

%upslope = boolean with all upstream cells marked as 1 and others nan.
szUp = [numel(lat) numel(lon)];
upslope = nan(szUp);

colPt = min(abs(pt(1)-lon));
rowPt = min(abs(pt(2)-lat));

indPt = sub2ind(szUp, rowPt, colPt);

indFlowF = find(fdr(:,indPt) == 1);

upslope(indFlowF) = 1;

iiLoop = 1;
while ~isempty(indFlowF)
    iiLoop = iiLoop + 1;
    
    flowTemp = fdr^ii;
    indFlowF = find(flowTemp(:,indFlowF) == 1);

    upslope(indFlowF) = 1;
end