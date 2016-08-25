function varargout = report_stats(sObs, sOutput, cellStats, dirStats, sMeta,varargin)


score = cell(numel(cellStats),1);
type = cell(numel(cellStats),1);
strScoreOut = cell(numel(cellStats),1);
strScoreDisp = cell(numel(cellStats),1);

blDisp = 1;
blWrt = 1;
if ~isempty(varargin)
    for ii = 1 : numel(varargin)
        if regexpbl(varargin{ii}, {'no','disp'}, 'and')
            blDisp = 0; 
        elseif regexpbl(varargin{ii}, {'no','write'}, 'and')
            blWrt = 0; 
        end
    end
end

% [score{1}, type{1}] = mod_v_obs_v2(sObs, sOutput, fitType,'combineType');
%     %If NSE used, convert to regular expression (1 - error); in 
%     %'fitness' NSE is calculate only as error)
%     if regexpbl(fitType, {'NSE','KGE'}) 
%         score{1} = 1 - score{1};
%     end
%     strScoreOut{1} = blanks(0);
% 
%     for ii = 1 : numel(score{1})
%         strScoreOut{1} = [strScoreOut{1}, num2str(round2(score{1}(ii),2)) ' for ' type{1}{ii}];
%         if ii ~= numel(score{1})
%             strScoreOut{1} = [strScoreOut{1}, ', '];
%         end
%     end
% strScoreDisp{1} = ['The ' fitType ' scores between modeled and observed '...
%     'values for the current run are: ' strScoreOut{1}];
% 
% mod_v_obs_v2(sObs, sOutput, fitType, 'plot', dirModObs,'combineType');

for yy = 1 : numel(cellStats)
    [score{yy}, type{yy}] = mod_v_obs_v2(sObs, sOutput, cellStats{yy},'combineType');
    %If NSE used, convert to regular expression (1 - error); in 
    %'fitness' NSE is calculate only as error)
    if regexpbl(cellStats{yy}, {'NSE','KGE'}) 
        score{yy} = 1 - score{yy};
    end
    strScoreOut{yy} = blanks(0);

    for ii = 1 : numel(score{yy})
        strScoreOut{yy} = [strScoreOut{yy}, num2str(round2(score{yy}(ii),2)) ...
            ' for ' type{yy}{ii}];
        if ii ~= numel(score{yy})
            strScoreOut{yy} = [strScoreOut{yy}, ', '];
        end
    end

    strScoreDisp{yy} = ['The ' cellStats{yy} ' scores between modeled '...
        'and observed values for the current run are: ' strScoreOut{yy}];
end
    warning('on','all')

  
if nargout > 0
    varargout{1} = score;

    if nargout > 1
        varargout{2} = type;
    end
end
    
% simple_text(fullfile(dirStats,'model_performance.txt'), ...
%     ['The model was calibrated using the ' cellStats{yy} '.' char(10)]);
for yy = 1 : numel(score)
    if blDisp == 1
        disp(strScoreDisp{yy});
    end
    if blWrt == 1
        simple_text(fullfile(dirStats,'model_performance.txt'), [strScoreDisp{yy} char(10)]);
    end
end


%Write output to CCHF formatted text file:
if blWrt
    nmOutputRoot = [sMeta.region, '_' sMeta.runType, '_' date2str(sMeta.dateStart,'_'), 'thru' date2str(sMeta.dateEnd,'_')];
    outputRoot = fullfile(dirStats, nmOutputRoot);
    pathOutput = [outputRoot '_modelled.txt'];
    write_CCHF_gagedata(pathOutput, sOutput);
end