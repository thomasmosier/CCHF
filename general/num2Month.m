function str = num2Month(n,varargin)

str = cell(numel(n),1);

for ii = 1 : numel(n)
    switch n(ii)
        case 1
            str{ii} = 'January';
        case 2
            str{ii} = 'February';
        case 3
            str{ii} = 'March';
        case 4
            str{ii} = 'April';
        case 5
            str{ii} = 'May';
        case 6
            str{ii} = 'June';
        case 7
            str{ii} = 'July';
        case 8
            str{ii} = 'August';
        case 9
            str{ii} = 'September';
        case 10
            str{ii} = 'October';
        case 11
            str{ii} = 'November';
        case 12
            str{ii} = 'December';
        otherwise
            str{ii} = 'Unknown';
    end

    if ~isempty(varargin)
       str{ii} = str{ii}(1:varargin{1}); 
    end
end    