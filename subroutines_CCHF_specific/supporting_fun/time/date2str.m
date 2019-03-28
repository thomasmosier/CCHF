% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

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