function [dataGage, lon, lat, date] = read_Arc_data(path)
%Read Arc data that has been copied into Excel file.  
%Must have lon and lat attributes, date, and type


sepDate = {'\','/','-',':'};

%Get raw data:
[data, txtRaw] = xlsread(path);

%Skip all lines before column labels:
flagHdr = 1; cntrHdr = 0;
while flagHdr == 1
    cntrHdr = cntrHdr + 1;
    
    cntrCurrEmpt = 0;
    for ii = 1 : numel(txtRaw(1,:))
        if ~isempty(txtRaw{cntrHdr,ii})
            cntrCurrEmpt = cntrCurrEmpt + 1;
        end
    end
    
    if cntrCurrEmpt == numel(txtRaw(1,:))
       flagHdr = 0; 
    end
end
%Record column labels:
hdrLabels =  txtRaw(cntrHdr,:);

% %Find which columns have txt and which are numbers:
% colIsNum = nan(numel(txtRaw(cntrHdr+1,:)),1);
% for ii = 1 : numel(txtRaw(cntrHdr+1,:))
%     if isempty(txtRaw{cntrHdr+1,ii})
%         colIsNum(ii) = 1;
%   	end
% end

%Initialize column numbers:
colLon = [];
colLat = [];
colData = [];
colDate = [];
colDateStrt = [];
colDateEnd = [];

%Find columns corresponding to data, lat, lon, date
for ii = 1 : numel(hdrLabels(:))
   if regexpbl(hdrLabels{ii},{'lon','x_point','point_x'})
       colLon = ii;
   elseif regexpbl(hdrLabels{ii},{'lat','y_point','point_y'})
       colLat = ii;
   elseif regexpbl(hdrLabels{ii},{'balance','swe'})
       colData = ii;
   elseif regexpbl(hdrLabels{ii},{'date','time'})
       if regexpbl(hdrLabels{ii},'start')
            colDateStrt = ii;
       elseif regexpbl(hdrLabels{ii},'end')
            colDateEnd = ii;
       else
           colDate = ii;
       end
   else
       %Do nothing: just means it's a non-relevant column
   end
end

%Error messages if column header not found:
if isempty(colLon)
    error('read_Arc_data:colLon','No longitude column label was found.');
end
if isempty(colLat)
    error('read_Arc_data:colLat','No latitude column label was found.');
end
if isempty(colData)
    error('read_Arc_data:colData','No data column label was found.');
end
if isempty(colDate) && (isempty(colDateStrt) || isempty(colDateEnd))
    error('read_Arc_data:colDate','No date column label was found.');
end

%Put all data into text cell array (to ensure columns lineup)
cntrData = 1;
for ii = 1 : numel(txtRaw(1,:))
   if isempty(txtRaw{2,ii}) 
        while isnan(data(2,cntrData))
            cntrData = cntrData + 1;
        end
        txtRaw(2:end,ii) = num2cell(data(:,cntrData));
        cntrData = cntrData + 1; 
   end
end

%Get lat, lon, and data:
lon = cell2mat(txtRaw(2:end,colLon));
lat = cell2mat(txtRaw(2:end,colLat));
dataGage = cell2mat(txtRaw(2:end,colData));


%%Convert dates in text format to numeric array format:

%Retrieve either one date or two:
if ~isempty(colDateStrt) && ~isempty(colDateEnd)
    colDate = [colDateStrt, colDateEnd];
end
   
%Identify seperator used:
for ll = 1 : numel(sepDate(:))
    if regexpbl(txtRaw{cntrHdr+1,colDate(1)},sepDate{ll});
        sepUse = sepDate{ll}; 
    end
end

nColDate = numel(regexpi(txtRaw{cntrHdr+1,colDate(1)},sepUse)) + 1;
%Initialize output array:
date = cell(numel(colDate),1);
[date{:}] = deal(nan(numel(lon), nColDate));


for ii = 1 : numel(colDate)
    %Convert each date:
    for kk = 1 : numel(lon)
        rowDate = cntrHdr + kk;
        indSep = regexpi(txtRaw{rowDate, colDate(ii)}, sepUse);

        if nColDate == 2
            date{ii}(kk,:) = [ ...
                str2double(txtRaw{rowDate, colDate(ii)}(1:indSep(1)-1)), ...
                str2double(txtRaw{rowDate, colDate(ii)}(indSep(1)+1:end))];
        elseif nColDate == 3
            date{ii}(kk,:) = [ ...
                str2double(txtRaw{rowDate, colDate(ii)}(1:indSep(1)-1)), ...
                str2double(txtRaw{rowDate, colDate(ii)}(indSep(1)+1:indSep(2)-1)), ...
                str2double(txtRaw{rowDate, colDate(ii)}(indSep(2)+1:end))];
        else
            error('read_Arc_data:DateNumber',['There are ' ...
                num2str(nColDate) ' date values.  This has not '...
                'been programmed for.']);
        end
    end
    
    colYr   = find(all(date{ii} > 31));
    colMnth = find(all(date{ii} < 13) == 1);
    
    %Swap date columns:
    if nColDate == 2
        date{ii} = [date{ii}(:,colYr), date{ii}(:,colMnth)];
    elseif nColDate == 3
        colDay = find(all(date{ii} < 33) == 1 & any(date{ii} > 13) == 1);
        date{ii} = [date{ii}(:,colYr), date{ii}(:,colMnth), date{ii}(:,colDay)];
    else
        error('read_Arc_data:DateNumber',['There are ' ...
            num2str(numel(date{ii}(1,:))) ' date values.  This has not '...
            'been programmed for.']);
    end
    
    if numel(txtRaw(2:end,1)) ~= numel(date{ii}(:,1))
        error('read_Arc_data:diffLength',['There are ' ...
            num2str(numel(txtRaw(2:end,1))) ' entries in the data array and '...
            num2str(numel(date{ii}(:,1))) ' entries in the date array.  '...
            'The cause of this mismatch needs to be investigated.']);
    end
end

if numel(colDate) == 1
   date = date{1}; 
end
    





