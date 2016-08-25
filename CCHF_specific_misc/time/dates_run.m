function sMeta = dates_run(sMeta,varargin)



% dateRef = [0,1,1,1,1]; %year, month, day, hour, minute
cal = 'Gregorian';
%Adjust start and end date vector to have same temporal resolution as timestep
if regexpbl(sMeta.dt,{'daily','day'})
    if length(sMeta.dateStart) == 2
        sMeta.dateStart = [sMeta.dateStart, 1];
    elseif length(sMeta.dateStart) > 3
        sMeta.dateStart = sMeta.dateStart(1:3);
    elseif length(sMeta.dateStart) == 3
        
    else
        error('dates_run:shortDate',['The start date has a length that '...
            'has not been programmed for.']);
    end
    
    if length(sMeta.dateEnd) == 2
        sMeta.dateEnd = [sMeta.dateEnd, eomday(sMeta.dateEnd(1), sMeta.dateEnd(2))];
    elseif length(sMeta.dateEnd) > 3
        sMeta.dateEnd = sMeta.dateEnd(1:3);
    elseif length(sMeta.dateStart) == 3
            
    else
        error('dates_run:shortDate',['The end date has a length that '...
            'has not been programmed for.']);
    end
elseif regexpbl(sMeta.dt,'month')
    if length(sMeta.dateStart) == 1
        sMeta.dateStart = [sMeta.dateStart, 7];
    elseif length(sMeta.dateStart) == 2
           
    elseif length(sMeta.dateStart) > 2
        sMeta.dateStart = sMeta.dateStart(1:2);
    else
        error('dates_run:shortDate',['The start date has a length that'...
            ' has not been programmed for.']);
    end
    
    if length(sMeta.dateEnd) == 1
        sMeta.dateEnd = [sMeta.dateEnd, 6];
    elseif length(sMeta.dateStart) == 2
           
    elseif length(sMeta.dateEnd) > 2
        sMeta.dateEnd = sMeta.dateEnd(1:2);
    else
        error('dates_run:shortDate',['The end date has a length that'...
            ' has not been programmed for.']);
    end
elseif regexpbl(sMeta.dt,'hour')
    if length(sMeta.dateStart) == 3
        sMeta.dateStart = [sMeta.dateStart, 1];
    elseif length(sMeta.dateStart) == 2
        sMeta.dateStart = [sMeta.dateStart, 1, 1];
    elseif length(sMeta.dateStart) == 4
           
    elseif length(sMeta.dateStart) > 4
        sMeta.dateStart = sMeta.dateStart(1:4);
    else
        error('dates_run:shortDate',['The start date has a length that '...
            'has not been programmed for.']);
    end
    
    if length(sMeta.dateEnd) == 3
        sMeta.dateEnd = [sMeta.dateEnd, 24];
    elseif length(sMeta.dateEnd) == 2
        sMeta.dateEnd = [sMeta.dateEnd, eomday(sMeta.dateEnd(1), sMeta.dateEnd(2)), 24];
    elseif length(sMeta.dateStart) == 4
           
    elseif length(sMeta.dateEnd) > 4
        sMeta.dateEnd = sMeta.dateEnd(1:3);
    else
        error('dates_run:shortDate',['The end date has a length that '...
            'has not been programmed for.']);
    end
end

%Add spin to time:
if ~isempty(varargin(:)) && regexpbl(varargin{1},'spin')
    indSpin = regexpi(sMeta.spinup,'[0-9]');
    
    if isempty(indSpin)
        error('dates_run:spinupNumber',['No numbers detected in spinup string' char(39) sMeta.spinup char(39)]);
    elseif any(diff(indSpin)> 1)
        error('dates_run:spinupNumber',['Multiple numbers detected in spinup string' char(39) sMeta.spinup char(39)]);
    else
        nSpin = str2double(sMeta.spinup(indSpin(1):indSpin(end)));
    end
    
	if  regexpbl(sMeta.spinup,'month')
        if length(sMeta.dateStart) == 2
            dateStart = sMeta.dateStart - [0, nSpin];
        elseif length(sMeta.dateStart) == 3
            dateStart = sMeta.dateStart - [0, nSpin, 0];
        elseif length(sMeta.dateStart) == 4
            dateStart = sMeta.dateStart - [0, nSpin, 0, 0];
        else
            error('dates_run:spinupMonth',[length(sMeta.dateStart) ' for time unit of month is an unrecognized length.']);
        end
        
        %Correct for negative months:
        if dateStart(2) < 1
           yrExtraFrac = -dateStart(2)/12; 
           dateStart(1) = dateStart(1) - ceil(yrExtraFrac);
           dateStart(2) = mod(dateStart(2)-1,12)+1;
        end
	elseif regexpbl(sMeta.spinup,'year')
        if length(sMeta.dateStart) == 1
            dateStart = sMeta.dateStart - nSpin;
        elseif length(sMeta.dateStart) == 2
            dateStart = sMeta.dateStart - [nSpin, 0];
        elseif length(sMeta.dateStart) == 3
            dateStart = sMeta.dateStart - [nSpin, 0, 0];
        elseif length(sMeta.dateStart) == 4
            dateStart = sMeta.dateStart - [nSpin, 0, 0, 0];
        else
            error('dates_run:spinupMonth',[length(sMeta.dateStart) ' for time unit of year is an unrecognized length.']);
        end
	elseif regexpbl(sMeta.spinup,{'day','daily'})
        if length(sMeta.dateStart) == 3
            dateStart = sMeta.dateStart - [sMeta.spinup, 0];
        else
            error('dates_run:spinupMonth',[length(sMeta.dateStart) ' for time unit of year is an unrecognized length.']);
        end
        
        if dateStart(3) < 1 
            if eomday(sMeta.dateStart(1), mod(sMeta.dateStart(2)-1,12)+1) - dateStart(3) > 0
                sMeta.dateStart(2) = sMeta.dateStart(2) - 1;
                sMeta.dateStart(3) = eomday(sMeta.dateStart(1), mod(sMeta.dateStart(2)-1,12)+1) - dateStart(3);
            else
                error('dates_run:spinupDay','The spinup period is in units of days but greater than one month.  This has not been programmed for.')
            end
        end
	else
       error('dates_run:spinupUnit',[sMeta.spinup ' is an unrecognized time unit.']);
	end
    
else
    dateStart = sMeta.dateStart;
end

sMeta.dateRun = date_vec_fill(dateStart,sMeta.dateEnd,cal);

% if length(sMeta.dateStart) < 5
%     dateRef = dateRef(1:length(sMeta.dateStart));
% end