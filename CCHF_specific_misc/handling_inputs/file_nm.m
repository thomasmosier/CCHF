function fileNm = file_nm(region,var,dateVec)


dateStr = blanks(0);

for ii = 1 : numel(dateVec)
    dateStr = [dateStr, num2str(dateVec(ii))];
    if ii ~= numel(dateVec)
        dateStr = [dateStr, '.'];
    end
end

fileNm = [region '_' var '_' dateStr];