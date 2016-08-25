function dataOut = interp_month2day(datesMonth, dataIn, datesDays, method, blNorm)


daysInMonthIn = nan(numel(datesMonth(:,1)),1);
for ii = 1 : numel(daysInMonthIn)
   daysInMonthIn(ii) = eomday(datesMonth(ii,1), datesMonth(ii,2));
end
daysIn = [0; cumsum(daysInMonthIn(1:end-1))] + 0.5*daysInMonthIn;

daysInMonthOut = nan(numel(datesDays(:,1)),1);
cumDays = 0;
for ii = 1 : numel(daysInMonthOut)
    if ii ~= 1 && diff([datesDays(ii,2),datesDays(ii-1,2)]) ~= 0
        cumDays = cumDays + eomday(datesDays(ii-1,1), datesDays(ii-1,2));
    end
   daysInMonthOut(ii) = cumDays + datesDays(ii,3);
end

compareFirst = datesDays(1,:) < [datesMonth(1,:), 0.5*eomday(datesMonth(1,1), datesMonth(1,2))];

if isequal(compareFirst,[0,0,1])
    daysOffset = 0;
else
    error('interp_month2day:diffStartDates', ['The input data date and '...
        'interp2 date are not the same.  This has not been programmed for.']);
end

daysIn = daysIn + daysOffset;

if sum(~isnan(dataIn)) > 2
    dataOut = interp1(daysIn, dataIn, daysInMonthOut, method);
else
    dataOut = nan(numel(daysInMonthOut),1);
end

%Normalize if specified:
if blNorm == 1
    for ii = 1 : numel(dataIn)
        indCurr = ismember(datesDays(:,1:2), datesMonth(ii,:),'rows');
        
        scale = dataIn(ii)/nansum(dataOut(indCurr));
        scale(isnan(scale)) = 0;
        dataOut(indCurr) = scale*dataOut(indCurr);
    end
elseif blNorm ~= 0
    error('quad_fit_exact:normalize',['Normalization value is ' ...
        num2str(blNorm) ' but must be either 0 or 1.']);
end

% plot(daysInMonthOut, dataOut, '.', daysIn, dataIn, '*');