function raad = solar_daily_Dewalle(latGrid, aspect, slope, typeTopo, date, tRes, varargin)


% %Example from pg. 401 of DeWalle and Rango:
% latGrid = 55;
% dateMod = [2009, 2, 7];
% aspect = 135;
% slope = -30;
% typeTopo = 'slope';

I0 = nan;
if ~isempty(varargin(:))
    for ii = 1 : numel(varargin(:))
        if regexpbl(varargin{ii}, {'solar','constant'}, 'and')
            I0 = varargin{ii+1};
        end
    end
end

%Solar Constant:
if isnan(I0) %If not input, use default from DeWalle and Rango
    I0 =  (4.896*10^6)/86400; %Units = (W/hr/m^2)
elseif I0 > 1000 && I0 < 2000 %if units = Watts/m^2
    I0 = I0 / 24; %Convert to Watts/hr (because integrated over sunlight hours in day below)
end
%Dewalle and Rango use a value of 1360 W; Liston uses value of 1370 W

%Get unique days of the year to loop over (calculations same regardless of year)
%Values are assigned to input dates at the end of function
[~, indDate] = unique(date(:,2:end),'rows');
dateDayJ = date(indDate,:);

%Calculate Julian dates:
[dayJ, dateDayJ] = dates_2_Julian_days(dateDayJ, tRes);
    
%Declination angle:
angleDec = solar_dec_angle(dayJ);

%Eccentricity of Earth's orbit:
radVec = solar_rad_sq(dayJ);

%Correct latitude/lon to account for topography:
[latGridCorr, lonDiff] = solar_topo_slope(latGrid, aspect, slope, typeTopo);

%Calculate daylight hrs based on corrected topography
[hrRiseSlope, hrSetSlope, hrRiseAct, hrSetAct] ...
    = solar_daylight_hrs(angleDec, latGrid, latGridCorr, lonDiff);
%Calculate total hours:
hrs = hrSetSlope - hrRiseSlope;
hrs(hrs > 24) = 24;
hrs(hrs < 0) = 0;

angFreqRad = 3.8197; %Angular velocity of Earth (units = hrs per radian)
angVelDeg = 15; %Angular velocity of Earth (units = degrees per hr)
%Calculate daily radiation:
tempRaad = nan([numel(dayJ), size(latGridCorr)]);
for ii = 1 : numel(dayJ)
%     %Example from DeWalle and Rango Example:
%     cZTest = (3.76 - (-4.435))*sind(30.43)*sind(-15.583) + 3.8197*cosd(30.43)*cosd(-15.583)*(sind(80.57) - sind(-42.32));
%     raadTest = (I0./radVec(ii))*cZTest;
    %Integrated cos(Z) (Eq. B.4)
    cZenith = squeeze(hrs(ii,:,:)).*sind(latGridCorr).*sind(angleDec(ii)) ...
        + angFreqRad*cosd(latGridCorr).*cosd(angleDec(ii)).*squeeze(sind(angVelDeg*hrSetAct(ii,:,:))-sind(angVelDeg*hrRiseAct(ii,:,:))); %Hours
    %Daily Radiation (Units of Watts)
    %tempRaad(ii,:,:) = I0*radVec(ii).*cZenith/24;
    tempRaad(ii,:,:) = (I0./radVec(ii))*cZenith;
end

%Ensure no negative values
tempRaad(tempRaad < 0) = 0;


%Aggregate to desired time-step:
nTimeOut = numel(date(:,1));
raad = nan([nTimeOut, size(latGrid)], 'single');

datePrec = numel(date(1,:));
if regexpbl(tRes,{'day','daily'})
    if datePrec > 3
        warning('raadDailyDeWalle:dayMorePrecise', ...
            'Input date vector more precise than daily. Output has daily resolution.');
    elseif datePrec < 3
        error('raadDailyDeWalle:dayLessPrecise', ...
            'Input date vector less precise than daily, but daily output requested.')
    end
    
    for ii = 1 : nTimeOut
        raad(ii,:,:) = tempRaad(ismember(dateDayJ(:,2:3), date(ii,2:3), 'rows'), :, :);
    end
elseif regexpbl(tRes,'month')
    if datePrec > 2
        warning('raadDailyDeWalle:monthMorePrecise', ...
            'Input date vector more precise than monthly. Output has monthly resolution.')
    elseif datePrec < 2
        error('raadDailyDeWalle:monthLessPrecise', ...
            'Input date vector less precise than monthly, but monthly output requested.')
    end
    
    for ii = 1 : nTimeOut
        raad(ii,:,:) = squeeze(nanmean(tempRaad(ismember(dateDayJ(:, 2), date(ii, 2), 'rows'),:,:), 1)); 
    end
else 
    error('raadDailyDeWalle:unknownTimeStep', ['The time step ' tRes ' has not been programmed for.']);
end

