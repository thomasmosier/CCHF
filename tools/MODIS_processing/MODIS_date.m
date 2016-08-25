function dateVec = MODIS_date(fileNm)

%Dates are between the 2nd and 3rd periods:
indPer = regexpi(fileNm, '\.');

dateVec = [...
    str2double(fileNm(indPer(1)+2:indPer(1)+5)), ...
    str2double(fileNm(indPer(2)-3:indPer(2)-1))];