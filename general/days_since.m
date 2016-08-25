function nDays = days_since(startVec, endVec, strCal)
% Copyright 2013, 2014 Thomas M. Mosier, Kendra V. Sharp, and David F. Hill
% This file is part of the GlobalClimateData Downscaling Package.
% 
% The GlobalClimateData Downscaling Package is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The Downscaling Package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

%Input:
%Starting and ending date vectors with the format '[year, month, day]' (can
%also have [year, month, day, hours, minutes]
%strCal = 'julian', 'gregorian' (also 'standard'), '360', '365' (also 'noleap')
        
%Leap year calculations in Gregorian calendar:
%A year is a Leap Year if:
    %the year can be evenly divided by 100, it is NOT a leap year, unless;
    %The year is also evenly divisible by 400. Then it is a leap year.
%In Julian calendar, leap year happens every year where mod(yr,4) == 0
    %also, additional 10 days October 1582 relative to Gregorian


%Output:
%'nDays' = number of days between dates (i.e. difference between dates)

nDec = nan; %Used for rounding at end
%If inputs nan, return nan
if all(all(isnan(endVec))) || all(all(isnan(startVec)))
   nDays = nan(numel(endVec(:,1)));
   return
end

if sum(isnan(startVec)) > 0
    error('daysSince:refVec','The date reference vector contains a NaN value.')
end

if length(startVec(1,:)) == 2
    startVec = [startVec, ones(length(startVec(:,1)),1)];
end
if length(endVec(1,:)) == 2
    endVec = [endVec, ones(length(endVec(:,1)),1)];
    nDec = 0;
end


%If inputs have hour and/or minutes colum, convert to day decimal
if numel(startVec(1,:)) == 4
    startVec = [startVec(:,1:2), startVec(:,3) + startVec(:,4)/24];
    nDec = 4;
elseif numel(startVec(1,:)) == 5
    startVec = [startVec(:,1:2), startVec(:,3) + startVec(:,4)/24 + startVec(:,5)/(60*24)];
    nDec = 6;
end
if numel(endVec(1,:)) == 4
    endVec = [endVec(:,1:2), endVec(:,3) + endVec(:,4)/24];
    nDec = 4;
elseif numel(endVec(1,:)) == 5
    endVec = [endVec(:,1:2), endVec(:,3) + endVec(:,4)/24 + endVec(:,5)/(60*24)];
    nDec = 6;
end

nDays = NaN(length(endVec(:,1)),1);

%Converts calendar name into number corresponding to class:
blCal = cal_class(strCal);

strFormat = 'mm/dd/yyyy';

if length(startVec(:,1)) > 1 && length(startVec(1,:)) == 1
    startVec = startVec';
end
if length(endVec(:,1)) > 1 && length(endVec(1,:)) == 1
    endVec = endVec';
end

if blCal == 2
    nDays = 360*(endVec(:,1) - startVec(:,1));
    nDays = nDays + 30*(endVec(:,2) - startVec(:,2));
    nDays = nDays + endVec(:,3) - startVec(:,3);
    return;
end

if length(startVec(:,1)) ~= 1
   error('days_between:incorrectStartDate', ...
       ['The date vector defining the beginning date for the '...
       'interval cannot have more than one date.']);
end

for ii = 1 : length(nDays(:))
    %If a year or month element of the endVec contains a nan value, resulting
    %value is nan
    if sum(isnan(endVec(:,1:2))) > 0
        nDays = NaN;
        continue;
    end
    
    intDayStart = floor(startVec(3));
    decDayStart = startVec(3) - intDayStart;
    intDayEnd = floor(endVec(ii,3));
    decDayEnd = endVec(ii,3) - intDayEnd;
    
    strStart = [num2str( startVec(2)) '/' num2str( intDayStart) '/' num2str( startVec(1))];
    strEnd   = [num2str(endVec(ii,2)) '/' num2str( intDayEnd) '/' num2str(endVec(ii,1))];
    
    nDays(ii) = datenum(strEnd,strFormat) - datenum(strStart,strFormat) ...
        + decDayEnd - decDayStart;
    
    %Julian and Gregorian calendars do not include year 0 but datenum does, so subtract
    if startVec(1) < 0 && endVec(1) > 0
        nDays(ii) = nDays(ii) - (datenum('12/31/0000',strFormat) - datenum('1/1/0000',strFormat)+1);
    end

    %Convert from Serial number (essentially Julian date) to other calendars
    if blCal == 0
        %Subtract ten days if time interval covers October 1582:
        sub10 = 0;
        if sum(startVec(1) < [1582,10,4]) == 3 && sum(endVec(ii,1) > [1582,10,15]) == 3
            sub10 = 1;
        end
        
        if sub10 == 1
            nDays(ii) = nDays(ii) - 10;
        end
    elseif blCal == 3 %For Gregorian, have to add 10 days and then subtract a few leap years:
        if startVec(2) <= 2 && endVec(ii,2) >= 2 %Two extra potential leap years (first and ending year):
           nYrs = endVec(ii,1) - startVec(1) + 1;
           strtLoop = startVec(1);
        elseif startVec(2) > 2 && endVec(ii,2) > 2  %One extra potential leap years (ending year):
           nYrs = endVec(ii,1) - startVec(1);
           strtLoop = startVec(1) + 1;
        elseif startVec(2) <= 2 && endVec(ii,2) < 2  %One extra potential leap years (first year):
           nYrs = endVec(ii,1) - startVec(1);
           strtLoop = startVec(1);
        elseif startVec(2) > 2 && endVec(ii,2) <= 2 %One extra potential leap years
           nYrs = endVec(ii,1) - startVec(1) - 1;
           strtLoop = startVec(1) + 1;
        else %Case not written
            error('days_since:LeapYrs',['No function has been written '...
                'for the case when the starting month is ' ...
                num2str(startVec(2)) ' and the ending month is ' ...
                num2str(endVec(2))]);
        end
        
        if nYrs > 0
            nLeap = 0;
            for jj = 1 : nYrs
                if mod(strtLoop + jj - 1, 4) == 0 && eomday(strtLoop + jj - 1, 2) == 28
                    nLeap = nLeap + 1;
                end 
            end

            nDays(ii) = nDays(ii) + nLeap;
        end
        

    elseif blCal == 1 %Subtract number of leap years between dates in Julian Calendar:
        if startVec(2) <= 2 && endVec(ii,2) >= 2 %Two extra potential leap years (first and ending year):
           nYrs = endVec(ii,1) - startVec(1) + 1;
           strtLoop = startVec(1);
        elseif startVec(2) > 2 && endVec(ii,2) > 2  %One extra potential leap years (ending year):
           nYrs = endVec(ii,1) - startVec(1);
           strtLoop = startVec(1) + 1;
        elseif startVec(2) <= 2 && endVec(ii,2) < 2  %One extra potential leap years (first year):
           nYrs = endVec(ii,1) - startVec(1);
           strtLoop = startVec(1);
        elseif startVec(2) > 2 && endVec(ii,2) <= 2 %One extra potential leap years
           nYrs = endVec(ii,1) - startVec(1) - 1;
           strtLoop = startVec(1) + 1;
        else %Case not written
            error('days_since:LeapYrs',['No function has been written '...
                'for the case when the starting month is ' ...
                num2str(startVec(2)) ' and the ending month is ' ...
                num2str(endVec(2))]);
        end

        nLeap = 0;
        for jj = 1 : nYrs
            if eomday(strtLoop + jj - 1,2) == 29
                nLeap = nLeap + 1;
            end
        end

        nDays(ii) = nDays(ii) - nLeap;
        %'360' day calendar dealt with above
%     elseif regexpbl(strCal,'360')
%         daysMnth = 30;
%         
%         nTs = (endVec(ii,1) - startVec(1))*12 + endVec(ii,2) - startVec(2);
%         
%         for kk = 1 : nTs
%             dateCurr = startVec + [0,1,0]*(kk - 1);
%             if dateCurr(2) == 13
%                 dateCurr(2) = 1;
%                 dateCurr(1) = dateCurr(1) + 1;
%             end
%             
%             %Terminate if current month and year are same as ending vec
%             if isequal(dateCurr(1:2), endVec(ii,1:2))
%                break 
%             end
%             
%                nDays(ii) = nDays(ii) - eomday(dateCurr(kk,1),dateCurr(kk,2)) + 30;  
%         end       
    else
        error('days_since:unkownCal', [char(39) strCal char(39) 'is an unknown calendar type.']);
    end
    
    if isnan(nDec)
        %Round output based on time resolution of inputs
        switch numel(endVec(1,:))
            case 1 || 2
                nDec = 0;
            case 3
                nDec = 1;
            case 4
                nDec = 4;
        end
    end

    nDays(ii) = round2(nDays(ii),nDec);
            
end