function plot_parameters(prm, prmBnds, strPrm, varargin)


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

if iscell(prm)
   prm = cell2mat(prm);
end

perPrm = 100*(prm - prmBnds(:,1))./(prmBnds(:,2) - prmBnds(:,1));

hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame]);
hBar = barh((1:numel(prm))',perPrm,'FaceColor',clrBrewer(2,:),'LineWidth',lnWd);

xlim([0,100]);

%Add text showing absolute value to plot:
hText = nan(numel(prm),1);
for ii = 1 : numel(prm)
    x1 = perPrm(ii) + 1;
    y1 = ii;
    str1 = [num2str(prm(ii)) '  [' num2str(prmBnds(ii,1)) '\rightarrow' num2str(prmBnds(ii,2)) ']' ]; %'\leftarrow sin(\pi) = 0';
    hText(ii) = text(x1,y1,str1);
end
    

%Change X labels:
if numel(strPrm(1,:)) > 1
    strPrm = strPrm(:,1);
end
strPrm = strrep(strPrm,'_',' ');
set(gca,'YTick',(1:numel(prm)));
set(gca,'YTickLabel', strPrm);

hXLab = xlabel('Percent of Designated Parameter Bounds');
hYLab = ylabel('Model Parameters (units vary)');

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
    'fontSize'    , ftSzAx, ...
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

