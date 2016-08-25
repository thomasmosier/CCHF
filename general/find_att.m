function val = find_att(listAtt,strAtt, varargin)

flagWarning = 1;
if ~isempty(varargin)
    for ii = 1 : numel(varargin)
        if regexpbl(varargin{ii},{'no','warning'},'and')
            flagWarning = 0;
        end
    end
end
    
ind = strcmpi(listAtt,strAtt);

if sum2d(ind) == 1
    [row, col] = find(ind == 1);
    val = listAtt{row, col+1};
elseif sum2d(ind) > 1
    error('find_att:mult',['Multiple attributes labeled ' char(39) strAtt char(39) ' were found.']);
elseif sum2d(ind) == 0
    if flagWarning
        warning('find_att:none',['No attributes labeled ' char(39) strAtt char(39) ' were found.']);
    end
    
    val = [];
else
    val = [];
end