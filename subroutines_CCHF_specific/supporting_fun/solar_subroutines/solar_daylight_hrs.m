function [hrRiseSlope, hrSetSlope, hrRiseNoShift, hrSetNoShift] = solar_daylight_hrs(angleDec, latFlat, latSlope, lonDiff)


%Ensure lat is grid rather than vector
if ~isequal(size(latSlope), size(lonDiff))
    [~, latFlat] = meshgrid(lonDiff, latFlat);
    [~, latSlope] = meshgrid(lonDiff, latSlope);
end

latSame = isequal(latFlat, latSlope);

%Initialize output
hrRiseSlope = nan([numel(angleDec), size(latSlope)], 'single');
hrSetSlope  = nan([numel(angleDec), size(latSlope)], 'single');

hrRiseFlat = nan([numel(angleDec), size(latSlope)], 'single');
hrSetFlat  = nan([numel(angleDec), size(latSlope)], 'single');

angVel = 15; %degrees per hour
for ii = 1 : numel(angleDec)
    hrRiseSlope(ii,:,:) = real(-acosd(-tand(latSlope)*tand(angleDec(ii))) - lonDiff)/angVel;
    hrSetSlope(ii,:,:)  = real(+acosd(-tand(latSlope)*tand(angleDec(ii))) - lonDiff)/angVel;
    
    hrRiseFlat(ii,:,:) = real(-acosd(-tand(latFlat)*tand(angleDec(ii))))/15;
    hrSetFlat(ii,:,:)  = real(+acosd(-tand(latFlat)*tand(angleDec(ii))))/15;
end

%Ensure that the sloped times have equal or lower magnitudes to flat times
if ~latSame
    indRiseFlat = find(hrRiseSlope < hrRiseFlat);
    hrRiseSlope(indRiseFlat) = hrRiseFlat(indRiseFlat);

    indSetFlat = find(hrSetSlope > hrSetFlat);
    hrSetSlope(indSetFlat) = hrSetFlat(indSetFlat);
end

%Ensure no impossible values:
hrRiseSlope(hrRiseSlope < -12) = -12;
hrRiseSlope(hrRiseSlope >   0) =   0;

hrSetSlope(hrSetSlope > 12) = 12;
hrSetSlope(hrSetSlope <  0) =  0;

hrRiseNoShift = nan([numel(angleDec), size(latSlope)], 'single');
hrSetNoShift = nan([numel(angleDec), size(latSlope)], 'single');
for ii = 1 : numel(angleDec)
    hrRiseNoShift(ii,:,:) =  squeeze(hrRiseSlope(ii,:,:)) + lonDiff/angVel;
    hrSetNoShift(ii,:,:)  =  squeeze( hrSetSlope(ii,:,:)) + lonDiff/angVel;
end
% hrRiseFlat(hrRiseFlat < -12) = -12;
% hrRiseFlat(hrRiseFlat >   0) =   0;
% 
% hrSetFlat(hrSetFlat > 12) = 12;
% hrSetFlat(hrSetFlat <  0) =  0;
