function dateFilled = date_vec_fill(dateStart,dateEnd,cal)
%Assumes Gregorian calendar

if ~regexpbl(cal,'greg')
   error('date_vec_fill:unkownCal',[cal ' has not been programmed for.  '...
       'Currently this function only works with the Gregorian Calendar.']); 
end

if numel(dateStart(1,:)) ~= numel(dateEnd(1,:))
    if numel(dateStart(1,:)) < numel(dateEnd(1,:))
        strCrop = 'starting';
        dateStart = dateStart(1:numel(dateEnd(1,:))); 
    else
        strCrop = 'ending';
        dateEnd = dateEnd(1:numel(dateStart(1,:))); 
    end
    
    warning('date_vec_fill:diffLength',['The time resolution of the '...
        'input vectors is different.  The vectors are being cropped '...
        'to the resolution of the ' strCrop ' vector.']);
end

if numel(dateStart(1,:)) == 1
    nSteps = dateEnd - dateStart + 1;
elseif numel(dateStart(1,:)) == 2
    nSteps = 12*(dateEnd(1) - dateStart(1) + 1);
elseif numel(dateStart(1,:)) == 3
    nSteps = 366*(dateEnd(1) - dateStart(1) + 1);
elseif numel(dateStart(1,:)) == 4
    nSteps = 8784*(dateEnd(1) - dateStart(1) + 1);
end


dateFilled = nan(nSteps,numel(dateStart));
dateFilled(1,:) = dateStart;

if numel(dateStart(1,:)) == 1
    for ii = 1 : dateEnd - dateStart 
        dateFilled(ii+1,:) = dateStart + ii;
    end
elseif numel(dateStart(1,:)) == 2
    ii = 1;
    while ~isequal(dateFilled(ii,:), dateEnd)
        if dateFilled(ii,2) ~= 12
           dateFilled(ii+1,:) = dateFilled(ii,:) + [0,1];
        else
            dateFilled(ii+1,:) = [dateFilled(ii,1) + 1, 1];
        end

        ii = ii + 1;
    end
elseif numel(dateStart(1,:)) == 3
    ii = 1;
    while ~isequal(dateFilled(ii,:), dateEnd)
        if dateFilled(ii,3) ~= eomday(dateFilled(ii,1),dateFilled(ii,2))
           dateFilled(ii+1,:) = dateFilled(ii,:) + [0,0,1];
        else
            if dateFilled(ii,2) == 12
                dateFilled(ii+1,:) = [dateFilled(ii,1) + 1, 1, 1];
            else
                dateFilled(ii+1,:) = [dateFilled(ii,1), dateFilled(ii,2) + 1, 1];
            end
        end

        ii = ii + 1;
    end
elseif numel(dateStart(1,:)) == 4
    ii = 1;
    while ~isequal(dateFilled(ii,:), dateEnd)
        if dateFilled(ii,4) ~= 24;
            dateFilled(ii+1,:) = dateFilled(ii,:) + [0,0,0,1];
        else
            if dateFilled(ii,3) ~= eomday(dateFilled(ii,1),dateFilled(ii,2))
               dateFilled(ii+1,:) = [dateFilled(ii,1:3) + [0,0,1], 1];
            else
                if dateFilled(ii,2) == 12
                    dateFilled(ii+1,:) = [dateFilled(ii,1) + 1, 1, 1, 1];
                else
                    dateFilled(ii+1,:) = [dateFilled(ii,1), dateFilled(ii,2) + 1, 1, 1];
                end
            end
        end
        
        ii = ii + 1;
    end
end

dateFilled( isnan(dateFilled(:,1)), :) = [];