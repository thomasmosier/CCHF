function dateStr = date_range_str(dateStart, dateEnd)

[dateStart, strStartMult] = date_min(dateStart);
[dateEnd  ,   strEndMult] = date_max(dateEnd);

dateStr = [strStartMult num2str(dateStart(1)) 'thru' strEndMult num2str(dateEnd(1))];