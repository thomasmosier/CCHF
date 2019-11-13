% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function varargout = report_stats_v2(sObs, sOutput, cellStats, dirStats, sMeta, varargin)


scoreAll = cell(numel(cellStats),1);
type = cell(numel(cellStats),1);
strScoreOut = cell(numel(cellStats),1);
strScoreDisp = cell(numel(cellStats),1);

blDisp = 1;

if isfield(sMeta, 'wrtGridEval')
    blWrt = sMeta.wrtGridEval;
else
    blWrt = 1;
end

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

% mod_v_obs_v2(sObs, sOutput, fitType, 'plot', dirModObs,'combineType');
[scoreAll, typeAll] = mod_v_obs_v2(sObs, sOutput, cellStats,'combineType');

for yy = 1 : numel(scoreAll(:))
%     [score{yy}, type{yy}] = mod_v_obs_v2(sObs, sOutput, cellStats{yy},'combineType');
    %If NSE used, convert to regular expression (1 - error); in 
    %'fitness' NSE is calculate only as error)
    if regexpbl(cellStats{yy}, {'NSE','KGE'}) 
        scoreAll{yy} = 1 - scoreAll{yy};
    end
    strScoreOut{yy} = blanks(0);

    for ii = 1 : numel(scoreAll{yy})
        strScoreOut{yy} = [strScoreOut{yy}, num2str(round2(scoreAll{yy}(ii),2)) ...
            ' for ' typeAll{yy}{ii}];
        if ii ~= numel(scoreAll{yy})
            strScoreOut{yy} = [strScoreOut{yy}, ', '];
        end
    end

    strScoreDisp{yy} = ['The ' cellStats{yy} ' scores between modeled '...
        'and observed values for the current run are: ' strScoreOut{yy}];
end
    warning('on','all')

  
if nargout > 0
    varargout{1} = scoreAll;

    if nargout > 1
        varargout{2} = typeAll;
    end
end
    
% simple_text(fullfile(dirStats,'model_performance.txt'), ...
%     ['The model was calibrated using the ' cellStats{yy} '.' char(10)]);
for yy = 1 : numel(scoreAll(:))
    if blDisp == 1
        disp(strScoreDisp{yy});
    end
    if blWrt == 1
        simple_text(fullfile(dirStats,'model_performance.txt'), [strScoreDisp{yy} char(10)]);
    end
end


%Write output to CCHF formatted text file:
if blWrt
    %Set current model run dates
    dateStr = date_range_str(sMeta.dateStart, sMeta.dateEnd);

    nmOutputRoot = [sMeta.runType, '_' dateStr];
    outputRoot = fullfile(dirStats, nmOutputRoot);
    pathOutput = [outputRoot '_modelled.txt'];
    write_CCHF_gagedata(pathOutput, sOutput);
end