function varargout = bcalbedo_Ming(varargin)

global sCryo sAtm


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


%The basic premise of the formulation here is:
%(1) BC deposition since last snowfall affects snow albedo
%(2) Total accumulated BC deposition affects ice albedo. Of note, ice
    %albedo only matters once all snow has melted. In reality
    %there is an intermediary step of BC trapped inside the snowpack that
    %eventually melts and becomes deposited on the glacier. But, tracking 
    %this intermediary doesn't have any impact on the surface properties 
    %since ice is buried when it is snow covered. 


%Get global albedo values (fresh snow and typical ice)
aSnowFresh = find_att(varargin{1}.global,'albedo_snow_fresh');
aIce       = find_att(varargin{1}.global,       'albedo_ice');
% aDebris    = find_att(varargin{1}.global,    'albedo_debris');

%Assume density of fresh snow is ~350 kg/m^3
dSnow = 350; %Based on mid-point density of fresh snow on Yala Glacier
%Assume density of ice is ~850 kg/m^3
dIce = 850; %Based on recommendation from Tobias Bolch

%convert from micrograms per m^3 to micrograms per m^2 (based on Ming Fig
%3A)
sclBc = 0.0120; %=24/2000;

%From Ming et al. (2009): (snow albedo reduction) = m*x + b, where 
m = 0.0757; %Units = kg/micrograms
b = 0.0575; %Unitless

Ming = @(x) m*x+b;

%Set Snow depth limit:
dSnwLmt = 0.02; %2 cm (taken from Yasunari et al. 2010)]
dIceLmt = 0.1; %10 cm (this is hypothesis based on the holes that bc creates in ice)
%Set maximum BC value:
bcMax = 500; %Units = micrograms per kg
%This threshold is maximum BC deposition recorded in Table 3 of:
%Ming, J., Xiao, C., Du, Z., & Yang, X. (2013). An overview of black carbon 
%deposition in High Asia glaciers and its impacts on radiation balance. 
%Advances in Water Resources, 55, 80-87.

%Ming empirical formulation assumes "x" has units of micrograms per kg. The
%input is new BC deposition in units of micrograms per day per meter
%squared
%Ming's "x" is obtained by:
%      x        =     (new BC dep)    * (new snow density)^-1 * (new snow density)^-1
%[microgram/kg] = [microgram/day/m^2] *       [m^3/kg]        *        [day/m]

varTopSnow = 'sntop';
varBcSnow = 'snbctop';
varSnAlb = 'snalb';
varSnAlbNoBC = 'snalbnobc';
varIcAlb = 'icalb';
varIcAlbNoBC = 'icalbnobc';
varBcIce = 'icbc';
varAlbRedTop = 'albRedTop';
varAlbRedIce = 'albRedIce';

%Set indices of entire spatial grid
indAll = (1:numel(sAtm.prsn));

%New black carbon:
newBC = squeeze(sAtm.bcdep(sAtm.indbcdep,:,:));

%%Calculate albedo reduction of snow:
%Find indices of new snow (this will be the snow depth used in the
%calculation)
indNew = find(sAtm.prsn > 0);

%Initialize top snow layer 
if ~isfield(sCryo, varTopSnow) 
    indNew = indAll;
    sCryo.(varTopSnow) = zeros(size(sAtm.prsn));
    sCryo.(varBcSnow)  = zeros(size(sAtm.prsn));
end

%Set depth of new snow layer
sCryo.(varTopSnow)(indNew) = sAtm.prsn(indNew);

%Set minimum to mixing/optical limit 
sCryo.(varTopSnow)(sCryo.(varTopSnow) < dSnwLmt) = dSnwLmt;

%Set black carbon deposition for new snow:
sCryo.(varBcSnow)(indNew) = newBC(indNew);

%Find locations without new snow:
indOld = intersect(indAll, indNew);
%Add black carbon to these locations:
if ~isempty(indOld)
    sCryo.(varBcSnow)(indOld) = sCryo.(varBcSnow)(indOld) + newBC(indOld);
end

%Save calculated snow albedo that doesn't account for BC deposition
%(This prevents the bc aledo reduction from being applied again and again)
sCryo.(varSnAlbNoBC) = sCryo.(varSnAlb);

%Calculate "x" for Ming's equation:
%THINK ABOUT SCALING FACTOR
xSn = sclBc*sCryo.(varBcSnow)./(dSnow*sCryo.(varTopSnow));
%x = sclBc*sCryo.(varBcSnow)(1:10,1:10)./(dSnow*sCryo.(varTopSnow)(1:10,1:10));
%Range should be about 0 to 500 micrograms per kg (based on figures in
%papers)
xSn(xSn > bcMax) = bcMax;

%Calculate snow layer albedo reduction based on Ming's equation:
sCryo.(varAlbRedTop) = Ming(xSn)/100; %Units are fraction (same as albedo)
%Set limits:
sCryo.(varAlbRedTop)(sCryo.(varAlbRedTop) < 0) = 0;
sCryo.(varAlbRedTop)(sCryo.(varBcSnow) == 0) = 0;
sCryo.(varAlbRedTop)(isinf(sCryo.(varAlbRedTop))) = 0;

%For diagnostics:
%xSn(1:10,1:10)
%max(xSn(:))
%sCryo.(varAlbRedTop)(1:10,1:10)

%Apply albedo reduction (fraction of original albedo):
sCryo.(varSnAlb) = sCryo.(varSnAlbNoBC).*(1-sCryo.(varAlbRedTop));

%Set max ice albedo for threshold
mxSnAlb = max(max(sCryo.(varSnAlbNoBC)(:)), aSnowFresh);

%Set albedo limits:
sCryo.(varSnAlb)(sCryo.(varSnAlb) > mxSnAlb) = mxSnAlb;
sCryo.(varSnAlb)(sCryo.(varSnAlb) < 0) = 0; 
%THINK ABOUT LOWER SATURATION LIMIT

%sCryo.snalbbc(isnan(sCryo.snalbbc)) = aUnder; 


%%Calculate albedo reduction of ice:
%Iniitialize ice BC grid
if ~isfield(sCryo, varBcIce) 
    sCryo.(varBcIce) = zeros(size(sAtm.prsn));
end

sCryo.(varBcIce) = sCryo.(varBcIce) + newBC;

%Save calculated ice albedo that doesn't account for BC deposition
%(This prevents the bc aledo reduction from being applied again and again)
sCryo.(varIcAlbNoBC) = sCryo.(varIcAlb);

%Calculate "x" for Ming's equation:
xIc = sclBc*sCryo.(varBcIce)./(dIceLmt*dIce);
%x = sclBc*sCryo.(varBcSnow)(1:10,1:10)./(dSnow*sCryo.(varTopSnow)(1:10,1:10));
xIc(xIc > bcMax) = bcMax;

%Calculate ice layer albedo reduction:
sCryo.(varAlbRedIce) = Ming(xIc)/100; %Units are fraction (same as albedo)
%Set limits:
sCryo.(varAlbRedIce)(sCryo.(varAlbRedIce) < 0) = 0;
sCryo.(varAlbRedIce)(sCryo.(varBcIce) == 0) = 0;
sCryo.(varAlbRedIce)(isinf(sCryo.(varAlbRedIce))) = 0;

%For diagnostics:
%xIc(1:10,1:10)
%max(xIc(:))
%sCryo.(varAlbRedIce)(1:10,1:10)

%Apply albedo reduction (fraction of original albedo):
sCryo.(varIcAlb) = sCryo.(varIcAlbNoBC).*(1-sCryo.(varAlbRedIce));

%Set max ice albedo for threshold
mxIcAlb = max(max(sCryo.(varIcAlbNoBC)(:)), aIce);

%Set albedo limits:
sCryo.(varIcAlb)(sCryo.(varIcAlb) > mxIcAlb) = mxIcAlb;
sCryo.(varIcAlb)(sCryo.(varIcAlb) < 0) = 0; 