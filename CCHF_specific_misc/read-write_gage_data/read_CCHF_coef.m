function [coef, hdr] = read_CCHF_coef(path)

[~, flNm, ~ ] = fileparts(path);
%Check to see if filename appears consistent with expected format:
if ~regexpbl(flNm, {'coefficients','coef'})
    warning('read_CCHF_coef:unexpectedNm',['The file ' char(39) flNm char(39) ' does not appear to be a coefficient file (based on naming).'])
end

fid = fopen(path,'rt');


hdr = cell(0,1);
for ii = 1 : 10
    hdr{ii,1}   = fgets(fid);
    if strcmpi(hdr{ii},char(10))
        rHdr = ii;
        break 
    end 
end
hdr = hdr(1:rHdr-1,:);

labels = fgets(fid);
nCol = numel(regexpi(labels, ',')) + 1;

dataStr = '%s';
for jj = 1 : nCol -1
    dataStr = [dataStr, '\t%s'];
end

cellRaw = textscan(fid, dataStr);

fclose(fid);

coef = cell(numel(cellRaw{1}), numel(cellRaw(1,:)));

for kk = 1 : numel(cellRaw{1})
    for ll = 1 : nCol
        coef{kk,ll} = cellRaw{ll}(kk);
        if iscell(coef{kk,ll})
            if ischar(coef{kk,ll}{1})
                coef{kk,ll} = coef{kk,ll}{1};
                if sum(isletter(coef{kk,ll})) == 0
                    coef{kk,ll} = str2double(coef{kk,ll});
                end
            end

        end
    end
    
    indComma = regexpi(coef{kk,1},',');
    if ~isempty(indComma)
        
        coef{kk,1} = coef{kk,1}(1:indComma(1)-1);
    end
end