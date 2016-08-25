function varargout = tLag_Liston(sHydro, varargin)

global sCryo sLand
%%Calculate velocity through each cell (from Liston's 'timecoefs.f')

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    varargout{1} = cat(1, varargout{1}, {'time_cf_slow', 0.2, 5, 2.778, 'tLag_Liston','routing'}); %Unitless scalar
    varargout{1} = cat(1, varargout{1}, {'time_cf_fast', 0.5, 5, 3.838, 'tLag_Liston','routing'}); %Unitless scalar
    
    return
else
    sMeta = varargin{1};
    alphaS = find_att(sMeta.coef,'time_cf_slow'); %Unitless scalar
    alphaF = find_att(sMeta.coef,'time_cf_fast'); %Unitless scalar
end



vSlowScIce = 0.12;
vSlowSfIce = 0.20;
vSlowScLand = 0.10;
vSlowSfLand = 0.08;


tFast = sparse(sHydro.fdr);
    tFast(sHydro.dl ~= 0) = nan;
tSlow = sparse(sHydro.fdr);
    tSlow(sHydro.dl ~= 0) = nan;



%Slope as degree:
if ~isfield(sHydro,'slopeLis') || ~isfield(sHydro,'sHydro.slopeCorr')
    sHydro.slopeLis = atand(sHydro.slope); %Slope is in rise/run, but Liston uses degrees
    sHydro.slopeLis(sHydro.slope < 5) = 5; 
    sHydro.slopeLis(sHydro.slope > 45) = 45;
    sHydro.slopeLis = sHydro.slopeLis.^0.8;
    sHydro.slopeLis(sHydro.dl == 0) = 0;
    %
    sHydro.slopeCorr = sHydro.slopeLis/(15^0.8);
end


%Calculate velocities at each point (m/s):
indDl = find(sHydro.dl ~= 0);
[rDl, ~] = ind2sub(size(sHydro.dl), indDl);
for ii = 1 : numel(indDl)
    [rC, cC] = ind2sub(size(sHydro.dem), rDl(ii));
    if sCryo.snw(rC,cC) > 0 && sCryo.icx(rC,cC) > 0 %Snow covered ice
        vFast = max(vSlowScIce, vSlowScIce*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowScIce;
    elseif sCryo.snw(rC,cC) == 0 && sCryo.icx(rC,cC) > 0 %Snow free ice
        vFast = max(vSlowSfIce, vSlowSfIce*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowSfIce;
    elseif sCryo.snw(rC,cC) > 0 && sCryo.icx(rC,cC) == 0 %Snow covered land
        vFast = max(vSlowScLand, vSlowScLand*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowScLand;
    else %Snow free land
        vFast = max(vSlowSfLand, vSlowSfLand*sHydro.slopeCorr(indDl(ii)));
        
        vSlow = vSlowSfLand;
    end
    
    %Populate residence time, 'tLag', array:
    tFast(indDl(ii)) = sHydro.dl(indDl(ii)) ./ vFast; %units of seconds
    tSlow(indDl(ii)) = sHydro.dl(indDl(ii)) ./ vSlow; %units of seconds
end
    
%Set cap for slow flow at one month: 
dtMonth = 24*3600*eomday(sMeta.datecurr(1), sMeta.dateCurr(2));
tSlow(tSlow > dtMonth) = dtMonth - 1;



%Adjust slow and fast flow times with fitting parameters: 
sLand.tlags = alphaS*tSlow; 
sLand.tLag = alphaF*tFast;