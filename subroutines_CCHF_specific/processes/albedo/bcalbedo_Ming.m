function varargout = bcalbedo_Ming(varargin)

global sCryo


if isempty(varargin(:))
    varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'albedo_fresh', 0.5, 1, 0.71, 'albedo_Pellicciotti','cryo'});

    return
else
	%sMeta = varargin{1};
end

%Accounts for black carbon impacts according to:
%Ming, J., Xiao, C., Cachier, H., Qin, D., Qin, X., Li, Z., & Pu, J. (2009). 
%Black Carbon (BC) in the snow of glaciers in west China and its potential 
%effects on albedos. Atmospheric Research, 92(1), 114?123. 
%https://doi.org/10.1016/j.atmosres.2008.09.007

%Fig 7: Reduced albedo (%) = 0.0757 x black carbon (micrograms per kilogram) + 0.0575
%where "Reduced albedo (%)" = -100*(A1-A0)/A0 = (100/A0)*(A0-A1) 
%with A0 = original albedo and A1 = new albedo
%Therefore, A1 = A0 - (A0/100)*(0.0757 x black carbon + 0.0575)
%This is equivalent to A1 = A0*(1 - (1/100)*(0.0757 x black carbon + 0.0575)

%See also discussion (for example Table 4) in:
%Yasunari, T. J., Bonasoni, P., Laj, P., Fujita, K., Vuillermoz, E., 
%Marinoni, A., ? Lau, K.-M. (2010). Estimated impact of black carbon 
%deposition during pre-monsoon season from Nepal Climate Observatory ? 
%Pyramid data and snow albedo changes over Himalayan glaciers. Atmospheric 
%Chemistry and Physics, 10(14), 6603?6615. https://doi.org/10.5194/acp-10-6603-2010



%Get albedo of ice
aFresh = find_att(varargin{1}.global,'albedo_snow_fresh');

%From Ming et al. (2009):
m = 0.0757;
b = 0.0575;

albRedTop = 1 - (m*sCryo.snbctop + b)/100;
albRedTot = 1 - (m*sCryo.bctot   + b)/100;

%%Snow albedo calculation:
sCryo.snalbbc = albRedTop*sCryo.snalb;

%Set albedo limits
sCryo.snalbbc(sCryo.snalbbc > aFresh) = aFresh;
sCryo.snalbbc(sCryo.snalbbc < 0) = 0;
%sCryo.snalbbc(isnan(sCryo.snalbbc)) = aUnder; 

%%Ice albedo calculation:
sCryo.icalbbc = albRedTot*sCryo.icalb;

%Set albedo limits
sCryo.snalbbc(sCryo.snalbbc > aFresh) = aFresh;
sCryo.snalbbc(sCryo.snalbbc < 0) = 0;