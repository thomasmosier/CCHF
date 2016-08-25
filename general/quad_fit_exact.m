function [arrayOut, daysInterp] = quad_fit_exact(sData, var, datesInterp, blNorm)
%INPUTS:
%sData = time-series data for three consecutive time-steps in geodata format
%DatesInterp = vector of matrix of dates to interpolate to in form [year1, month1, day1; year2, month2, day2; ...] 
%blNorm = if 1, normalizes output to have total value of middle input (used
%for precipitation, so that e.g. daily values add to monthly values);
%     
% keyboard

if ~isfield(sData,var)
    arrayOut   = nan;
    daysInterp = nan;
    
    return
end
    
if numel(sData.time) < 3 %Improper input date format
    error('quad_fit_exact:fewTime',['The current data structure has ' ...
        num2str(numel(sData.time)) ' time-series elements but at least '...
        '3 are needed.']);
elseif sum2d(isnan(datesInterp)) > 0 || sum(isnan(sData.time)) > 0 %If input times are nan, return input data as output and exit function.
    arrayOut = sData.(var);
%     sDataOut.time = nan;
    return
else
    %Convert date vectors to days, using sData calendar and reference date:
    [tRef, ~] = NC_time_units(sData.attTime);
    strCal =  NC_cal(sData.attTime);
    daysInterp = days_since(tRef, datesInterp, strCal);
    
    if ~issorted(daysInterp)
        daysInterp = sort(daysInterp);
    end
    
    %Choose 3 time-series elements from sData to use:
    %Want mid-date and dates bracketing requested outputs.
    indUse = nan(3,1);
    diffDays = indUse;
    indUse(1) = find(sData.time < daysInterp(1),1,'last');
    diffDays(1) = abs(sData.time(indUse(1)) - daysInterp(1));
    [diffDays(2), indUse(2)] = min(abs(sData.time - nanmean(daysInterp)));
    indUse(3) = find(sData.time > daysInterp(end),1,'first');
    diffDays(3) = abs(sData.time(indUse(3)) - daysInterp(end));
    
        diffThresh = 20;
    if sum(diffDays > diffThresh) > 0
       warning('quad_fit_exact:diffDates',[num2str(sum(diffDays > 20)) ...
           ' data points are more than ' num2str(diffThresh) ...
           ' time units apart from the expected date.']); 
    end
%     daysIn = sData.time(indUse);
end


%Initialize output array
if numel(size(sData.(var))) == 3
    arrayOut = nan([numel(daysInterp), size(squeeze(sData.(var)(1,:,:)))], 'single');
else
    arrayOut = nan([numel(daysInterp), size(sData.(var)(:,:))], 'single');
end


%USE EXACT-FIT QUADRATIC TO INTERPOLATE DAYS:
%These fitting parameters are matrices the same size as the
%climate grids
aFit = ( (sData.(var)(indUse(2),:,:) ...
        - sData.(var)(indUse(1),:,:))*(sData.time(indUse(1),:,:) ...
        - sData.time(indUse(3),:,:)) + (sData.(var)(indUse(3),:,:) ...
        - sData.(var)(indUse(1),:,:))*(sData.time(indUse(2),:,:) ...
        - sData.time(indUse(1),:,:)) ) ...
    / ( (sData.time(indUse(1),:,:) ...
        - sData.time(indUse(3),:,:))*(sData.time(indUse(2),:,:)^2 ...
        - sData.time(indUse(1),:,:)^2) + (sData.time(indUse(2),:,:) ...
        - sData.time(indUse(1),:,:))*(sData.time(indUse(3),:,:)^2 ...
        - sData.time(indUse(1),:,:)^2) );

bFit = ( (sData.(var)(indUse(2),:,:) ...
        - sData.(var)(indUse(1),:,:)) ...
        - aFit*(sData.time(indUse(2),:,:)^2 ...
        - sData.time(indUse(1),:,:)^2) ) ...
    / (sData.time(indUse(2),:,:) ...
        - sData.time(indUse(1),:,:));

cFit = (sData.(var)(indUse(1),:,:) ...
    - aFit*sData.time(indUse(1),:,:)^2 ...
    - bFit*sData.time(indUse(1),:,:));

%Apply coeffifient matrices to desired output dates:
for ii = 1 : numel(daysInterp)
    arrayOut(ii,:,:) = aFit*daysInterp(ii)^2 ...
        + bFit*daysInterp(ii) + cFit;
end
    
%Set negative precipitation to 0;
if regexpbl(sData.var,'pr')
        arrayOut(arrayOut < 0 ) = 0;   %Prevents negative precipitation
end

%Normalize if specified:
if blNorm == 1
    %Create scaling matrix for output data, which is middle time element divided by sum of time-series values at each grid point
    norm = squeeze(sData.(var)(indUse(2),:,:)) ./ squeeze(nansum(arrayOut(:,:,:),1));
    
    %Scale by middle time element of input grid: 
    for ii = 1 : numel(daysInterp)
        arrayOut(ii,:,:) = squeeze(arrayOut(ii,:,:)).*norm;
    end
elseif blNorm ~= 0
    error('quad_fit_exact:normalize',['Normalization value is ' ...
        num2str(blNorm) ' but must be either 0 or 1.']);
end
    

end