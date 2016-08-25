function varargout = runoff_bucket(sHydro, varargin)


%The processes of seperating between ground water, direct runoff, etc. are 
%taken from:
%Moore, R. D., Trubilowicz, J. W., & Buttle, J. M. (2012). Prediction of 
%Streamflow Regime and Annual Runoff for Ungauged Basins Using a 
%Distributed Monthly Water Balance Model1

%Consider alternative groundwater models, for example:
%Maxwell, R. M., & Miller, N. L. (2005). Development of a coupled land 
%surface and groundwater model. Journal of Hydrometeorology, 6(3), 233-247
%or
%Schaake, J. C., Koren, V. I., Duan, Q. Y., Mitchell, K., & Chen, F. 
%(1996). Simple water balance model for estimating runoff at different 
%spatial and temporal scales. Journal of Geophysical Research: Atmospheres 
%(1984�2012), 101(D3), 7461-7475

global sLand sCryo
%To evaulate function, varargin{1} = varargin{1}

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    varargout{1} = cat(1, varargout{1}, {'soil_moisture_capacity', 0, 4, 3.44, 'runoff_bucket','land'}); %Units of meters
    varargout{1} = cat(1, varargout{1}, {'drain_rate', 0, 0.034, 0.016, 'runoff_bucket','land'}); %Fraction drained per day (max = smc)
    return
else
	smc  = find_att(varargin{1}.coef,'soil_moisture_capacity');
	drain= find_att(varargin{1}.coef,            'drain_rate');  
end

%Drain rate is standardized for daily time-step; if dt is other, change:
if regexpbl(varargin{1}.dt,'month')
    drain = drain*eomday(varargin{1}.dateCurr(1), varargin{1}.dateCurr(2));
elseif regexpbl(varargin{1}.dt,'hour')
    drain = drain / 24;
elseif regexpbl(varargin{1}.dt,{'day','daily'})
    %Do nothing
else
    error('groundwater_bucket:dtUnknown',[varargin{1}.dt ' is an unknown time-step.'])
end

if drain > 1
   drain = 1; 
end


%Reset runoff field:
sLand.mrro = zeros(size(sLand.mrro));

indX = [];
if isfield(sCryo, 'ice')
   indX = find(sCryo.icx > 0);
end

%Create soil moisture field:
if ~isfield(sLand,'sm')
    sLand.sm = zeros(size(sLand.mrro));
end
%Set soil moisture holding capacity:
if ~isfield(sLand,'smc')
    sLand.smc = smc*ones(size(sLand.mrro));
    if ~isempty(indX) %Set moisture holding capacity to zero at glaciated cells
        sLand.smc(indX) = smc.*sCryo.icx(indX);
    end
end


%Calculate additions to soil moisture (rain, snowmelt, and icemelt):
sLand.sm = sLand.sm + sLand.rnrf + sCryo.snlr + sCryo.iclr;


%%Account for evapotranspiration losses from soil moisture 
if isfield(sLand,'pet')
    %Start by setting ET to PET
    indPET = find(ismember(sLand.datepet(:,2:end),varargin{1}.dateCurr(2:end), 'rows') == 1);
    if isempty(indPET)
        error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
    end
    sLand.et = squeeze(sLand.pet(indPET,:,:));
else
    error('direct_runoff:noPET','There should be a potential evapotranspiration field.');
end

%Set ET to 0 when snowpack or glaciers exist:
sLand.et(sCryo.snw > 0) = 0;
if ~isempty(indX)
    sLand.et(indX) = sCryo.icx(indX).*sLand.et(indX);
end

%Max ET is water stored in soil:
sLand.et(sLand.et > sLand.sm) = sLand.sm(sLand.et > sLand.sm); 
%Remove ET from ground storage:
sLand.sm = sLand.sm - sLand.et;


%%During each iteration, add percent of Soil moisture (equal to 'drain' %) to runoff:
%First add excess runoff
indExcess = find(sLand.sm > sLand.smc);
if ~isempty(indExcess)
   sLand.mrro(indExcess) = sHydro.area(indExcess).*(sLand.sm(indExcess) - sLand.smc(indExcess));
   sLand.sm(indExcess) = sLand.smc(indExcess);
end

%Second, add runoff from drain rate
sLand.mrro = sLand.mrro + drain*sHydro.area.*sLand.sm;
    sLand.mrro(sLand.mrro < 0 ) = 0;
sLand.sm = (1 - drain)*sLand.sm;
