function strDate = date2str(vecDate,varargin)

szDate = size(vecDate);
strDate = cell(szDate(1),1);

%define seperation symbol
sep = blanks(0);
form = blanks(0);

strSymb = '[-\/.,;:|_]';

if ~isempty(varargin(:))
    for jj = 1 : numel(varargin(:))
        if ~isempty(regexpi(varargin{jj},'[a-z]'))
            form = varargin{jj};
        elseif ~isempty(regexpi(varargin{jj},strSymb)) && isempty(regexpi(varargin{jj},'[a-z]'))
           sep = varargin{jj}; 
        end
    end
end
if isempty(sep) && ~isempty(form)
   sep = form(regexpi(form, strSymb, 'once'));
elseif isempty(sep)
    sep = '/';
end
if isempty(form)
    form = 'm/d/y'; 
end


indSep = regexpi(form,strSymb);
indSwap = ones(szDate(2),1);
for jj = 1 : szDate(2) 
    if jj == 1
        strComp = form(1:indSep(1)-1);
    elseif jj == szDate(2) 
        strComp = form(indSep(end)+1:end);
    else
        strComp = form(indSep(jj-1)+1: indSep(jj)-1);
    end
    
    if regexpbl(strComp, 'y')
        indSwap(jj) = 1;
    elseif regexpbl(strComp, 'm')
        indSwap(jj) = 2;
    elseif regexpbl(strComp, 'd')
        indSwap(jj) = 3;
    elseif regexpbl(strComp, 'h')
        indSwap(jj) = 4;
    end
end

for ii = 1 : szDate(1)
    switch szDate(2)
        case 3
            strDate{ii} = [num2str(vecDate(ii,indSwap(1))) sep num2str(vecDate(ii,indSwap(2))) sep num2str(vecDate(ii,indSwap(3)))];
        case 2
            strDate{ii} = [num2str(vecDate(ii,indSwap(1))) sep num2str(vecDate(ii,indSwap(2)))];
        case 4
            strDate{ii} = [num2str(vecDate(ii,indSwap(1))) sep num2str(vecDate(ii,indSwap(2))) sep num2str(vecDate(ii,indSwap(3))) sep num2str(vecDate(ii,indSwap(4)))];
        otherwise
            error('date2str:unknownFormat',['The input date vector has ' szDate(2) ' elements.  This has not been coded for.']);
    end
end

if numel(strDate(:)) == 1
   strDate = char(strDate); 
end