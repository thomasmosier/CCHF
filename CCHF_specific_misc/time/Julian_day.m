function day = Julian_day(currTime, varargin)

%currTime = [year, month, day, ...]
%day = Julian day

%Test relative to http://quasar.as.utexas.edu/BillInfo/JulianDateCalc.html

% %Start of Julian Calendar in Julian is Jan 1, 4713 BC:
% day = days_since([-4713, 1, 1],currTime(:,1:3),'gregorian'); 



%If inputs have hour and/or minutes colum, convert to day decimal
if numel(currTime(1,:)) == 4
    currTime = [currTime(:,1:2), currTime(:,3) + currTime(:,4)/24];
elseif numel(currTime(1,:)) == 5
    currTime = [currTime(:,1:2), currTime(:,3) + currTime(:,4)/24 + currTime(:,5)/(60*24)];
end

day = nan(numel(currTime(:,1)),1);

for ii = 1 : numel(day)
    if ~isempty(varargin(:)) && regexpbl(varargin{1},'ofyear')
       day(ii) = ...
           datenum([num2str(currTime(ii,2)) '/' num2str(currTime(ii,3)) '/' num2str(currTime(ii,1))], 'mm/dd/yyyy') ...
           - datenum(['01/01/' num2str(currTime(ii,1))], 'mm/dd/yyyy') ...
           + 1;
    else
        %Difference between date number and Julian date is 1721058.5
        if numel(currTime(1,:)) == 3
            day(ii) = 1721058.5 + ...
                datenum([num2str(currTime(ii,2)) '/' num2str(currTime(ii,3)) '/' num2str(currTime(ii,1))], 'mm/dd/yyyy');
        else
           error('Julian_date:inputLength',['The indices in the current input is ' num2str(numel(currTime(ii,:))) ', which has not been programmed for.']); 
        end
    end
end

