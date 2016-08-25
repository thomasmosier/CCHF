function [dirAvg, dirStd] = climatology(foldTs,climYrs,mnths,calcTyp,varargin)

if ~isempty(varargin(:))
    scl = varargin{1};  %Should be 1 except for PRISM temperature
else
    scl = 1;
end
filesTemp = dir(fullfile(foldTs,'*.asc'));
gridFiles = struct2cell(filesTemp);
gridFiles = gridFiles(1,:);

%Find indices seperating year and month (should be in format '*_year_mth.asc')
undInd = regexpi(gridFiles{1},'_');
if length(undInd) > 2
   undInd = undInd(end-1:end);
elseif length(undInd) < 2
    error('climFromTs:und',['File format does not contain expected number '...
        'of underscores.']);
end

if scl ~= 1
    disp(['Warning: A scale factor of ' num2str(scl) ' is being used.'])
end


%Create Output directory:
nmClimFold = ['clim' num2str(climYrs(1)) 'thru' num2str(climYrs(2))];
dirAvg = mk_output_dir(fullfile(foldTs, nmClimFold, gridFiles{1}), 'avg');
dirStd = mk_output_dir(fullfile(foldTs, nmClimFold, gridFiles{1}), 'std');

nYrs = (climYrs(2) - climYrs(1) + 1);
for ii = 1 : length(mnths) 
    filesCurrMnth = cellfun(@(x) ...
        str2double(x(undInd(end) + 1 : end-4)) == mnths(ii) ...
        && str2double(x(undInd(1) + 1 : undInd(end) - 1)) >= climYrs(1) ...
        && str2double(x(undInd(1) + 1 : undInd(end) - 1)) <= climYrs(2) ...
        , gridFiles);

    %Check that correct number of files selected:
    if sum(filesCurrMnth) ~= nYrs
       error('climFromTs:incomYrs',['An incorrect number of time ' ...
           'series files were selected from ' foldTs ' for use in ' ...
           'producing the current climatology.' ]); 
    end

    filesUse = gridFiles(filesCurrMnth);

    %Load files and calculate climatology:
    for mm = 1 : length(filesUse)
        [dataTsCurr, hdrCurr] = read_ESRI(fullfile(foldTs, filesUse{mm}));

        if mm == 1
            dataTs = zeros([nYrs,size(dataTsCurr)]);
            hdrFir = hdrCurr;
        else
            if sum(hdrFir == hdrCurr) ~= 6
                error('climFromTs:misalignHdr',['The header ' ...
                    'corresponding to ' filesUse{mm} ' does not align ' ...
                    'with that of ' filesUse{1} '.']);
            end
        end

        dataTs(mm,:,:) = dataTsCurr * scl;
    end

    if regexpbl(calcTyp,{'mean','average','avg'})
        dataClim = squeeze(mean(dataTs,1));
    elseif regexpbl(calcTyp,'sum')
        dataClim = squeeze(sum(dataTs,1));
    else
        error('climatology:calcType',[char(39) calcTyp char(39) ' is not a valid option.']);
    end
    dataStd = squeeze(std(dataTs,0,1));

    %Write climatology to file:
    nmClim = [char(filesUse{1}(1:undInd(1))) 'avg_' ...
        num2str(climYrs(1)) 'thru' num2str(climYrs(2)) '_' ...
        num2str(mnths(ii)) '.asc'];
    fullPathClim = fullfile(dirAvg, nmClim);
    disp([nmClim ' is being written.']);
    write_ESRI_v4(dataClim, hdrFir, fullPathClim, 3);
    disp([nmClim ' has been successfully written.']);

%     warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
%     fracComplt = ((jj-1)*(length(mnths)*nClim) + (kk-1)*length(mnths)+ii)/(nFold*nClim*length(mnths));
%     if fracComplt > 1
%        disp(['The progress fraction is ' num2str(fracComplt) '.']);
%        fracComplt = 1;
%     end
%     waitbar(fracComplt, hWait, ...
%         [strrep(nmClim,'_','\_') ' has been written.']);
%     warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');

    %Write standard deviation to file:
    nmClim = [char(filesUse{1}(1:undInd(1))) 'std_' ...
        num2str(climYrs(1)) 'thru' num2str(climYrs(2)) '_' ...
        num2str(mnths(ii)) '.asc'];
    fullPathStd = fullfile(dirStd, nmClim);
    disp([nmClim ' is being written.']);
    write_ESRI_v4(dataStd, hdrFir, fullPathStd, 3);
    disp([nmClim ' has been successfully written.']);
end