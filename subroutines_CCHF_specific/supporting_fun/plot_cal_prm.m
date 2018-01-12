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

function plot_cal_prm(prmIn, indStage, prmBnds, strPrm, varargin)


%%SET PROPERTIES:
%Set font size:
ftSz = 14;
    ftSzT = ftSz + 2;
    ftSzAx = ftSz - 2;
%Set line width:
lnWd = 2;
    lnWdA = lnWd - 0.5;
strFont = 'Arial';
%Use these colors instead?
%From color Brewer):
clrBrewer = [228,26,28; ...
    55,126,184; ...
    77,175,74; ...
    152,78,163]/255;
%Nice colors for two tone comparisons:
clrCust = [0.1,0.6,0.1; 0.1,0.1,0.7]; %Custom green and purple
% szFrame = [3.5,3];
szFrame = [12,7];

   
% %Set colors:
% if numel(numel(prm)) <= 4
%     colorsUse = clrBrewer;
% else
%     colorsUse = distinguishable_colors( numel(prm) );
% end

if ismatrix(prmIn) && ~iscell(prmIn)
    temp = prmIn;
    prmIn = cell(1,1);
    prmIn{1} = temp;
    indStage = cell(1,1);
    indStage{1} = (1:numel(temp(1,:)));
end

nStage = numel(prmIn(:));
nPrm = numel(prmIn{1}(1,:));
bestPrm = nan(nPrm, 1);
avgPrm = bestPrm;
sdPrm = avgPrm;
for ii = 1 : nStage
    indCurr = indStage{ii};
    
    bestPrm(indCurr) = prmIn{ii}(1,indCurr);
    avgPrm(indCurr) = nanmean(prmIn{ii}(:,indCurr));
    sdPrm(indCurr) = nanstd(prmIn{ii}(:,indCurr));
end
clear ii

perBestPrm = 100*(bestPrm - prmBnds(:,1))./(prmBnds(:,2) - prmBnds(:,1));
perAvgPrm = 100*(avgPrm - prmBnds(:,1))./(prmBnds(:,2) - prmBnds(:,1));

hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame]);
hold on
hBar = barh((1:nPrm)',perAvgPrm,'FaceColor',clrBrewer(2,:),'LineWidth',lnWd);
hPoint = scatter(perBestPrm, (1:nPrm)', 100, 'black', 'filled');
hErr = herrorbar(perAvgPrm, (1:nPrm)', sdPrm, sdPrm, 'oblack');
set(hErr, 'markerFaceColor', 'black','markerSize', 0.01, 'lineWidth', lnWd)
% hErr = errorbar(perBestPrm(:)', (1:numel(avgPrm)), sdPrm(:)','o','horizontal')
hold off
xlim([0,125]);
ylim([0.5, nPrm+ 0.5]);

set(gca, 'xtick', (0:20:100))
%Add text showing absolute value to plot:
hText = nan(numel(avgPrm),1);
for ii = 1 : numel(avgPrm)
    x1 = 100;
%     x1 = perAvgPrm(ii) + 1;
    y1 = ii;
    str1 = [num2str(round2(avgPrm(ii), 3)) '  [' num2str(round2(prmBnds(ii,1), 3)) '\rightarrow' num2str(round2(prmBnds(ii,2), 3)) ']' ]; %'\leftarrow sin(\pi) = 0';
    hText(ii) = text(x1,y1,str1);
end
    

%Change X labels:
if numel(strPrm(1,:)) > 1
    strPrm = strPrm(:,1);
end
strPrm = strrep(strPrm,'_',' ');
set(gca,'YTick',(1:numel(avgPrm)));
set(gca,'YTickLabel', strPrm);

hXLab = xlabel('Parameter Bounds (%)');
hYLab = ylabel('Parameters (units vary)');

%Set figure data properties:
set(hBar, ...
    'LineWidth', lnWd);
set([hXLab, hYLab, hText']  , ...
    'FontSize'   , ftSz, ...
    'Color', 'black', ...
    'FontName'   , strFont);
% set(hTtl, ...
%     'FontWeight' , 'bold', ...
%     'FontSize'   , ftSzT);
set(gca, ...
    'Box'         , 'off', ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'YMinorTick'  , 'off'      , ...
    'fontSize'    , ftSz, ...
    'LineWidth'   , lnWdA, ...
    'FontName'   , strFont);


%%WRITE FIGURE:
if ~isempty(varargin)
    pathPlot = varargin{1};
    
    [pathPlot, name, ext] = fileparts(pathPlot);

    if ~isempty(ext)
        pathPlot = fullfile(pathPlot, [name '_parameter_plot']);
    else
        pathPlot = fullfile(pathPlot, name, 'parameter_plot');
    end 

    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savefig(hFig, [pathPlot '.fig']);
    print(hFig,    [pathPlot '.eps'],'-depsc2');
    print(hFig,    [pathPlot '.png'],'-dpng','-r600');
end

