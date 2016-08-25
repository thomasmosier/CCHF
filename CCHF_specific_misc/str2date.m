function arrayDate = str2date(strDate)

%Input string format must be year/month/day/hour (with any symbol
%seperating)

strSymb = '[-\/.,;:|_]';

nDate = numel(regexpi(strDate{1},strSymb)) + 1;

arrayDate = nan(numel(strDate(:)),nDate);


for ii = 1 : numel(strDate(:))
    indCurr = regexpi(strDate{ii},strSymb);
    
    switch nDate
        case 1
           arrayDate(ii) = str2double(strDate{ii});
        case 2
            arrayDate(ii,:) = [...
                str2double(strDate{ii}(1:indCurr(1)-1)), ...
                str2double(strDate{ii}(indCurr(1)+1:end))];
        case 3
            arrayDate(ii,:) = [...
                str2double(strDate{ii}(1:indCurr(1)-1)), ...
                str2double(strDate{ii}(indCurr(1)+1:indCurr(2)-1)), ...
                str2double(strDate{ii}(indCurr(2)+1:end))];
        case 4
            arrayDate(ii,:) = [...
                str2double(strDate{ii}(1:indCurr(1)-1)), ...
                str2double(strDate{ii}(indCurr(1)+1:indCurr(2)-1)), ...
                str2double(strDate{ii}(indCurr(2)+1:indCurr(3)-1)), ...
                str2double(strDate{ii}(indCurr(3)+1:end))];
        otherwise
            error('str2date:nElements',['The date string has ' ...
                num2str(nDate) ' elements.  This case has not been programmed for.'])
    end
            
    
end