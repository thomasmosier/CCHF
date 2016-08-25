function varargout = toa_rad_DeWalle(date, sHydro, varargin)

%Formulation comes from
%DeWalle, D. R., & Rango, A. (2008). Principles of snow hydrology. 
%Cambridge University Press. (Appendix B)

global sAtm
%OUTPUT SCALING FITTING PARAMETER:
if isempty(varargin(:))
	varargout{1} = cell(0,6);
 
    return
else
    sMeta = varargin{1};
end


%NO FITTING PARAMETERS:
% if isempty(varargin(:))
% 	argout = cell(0,5);
%     return
% else
%     sMeta = varargin{1};  
% end



%Solar Constant (W/m^2):
I0 =  find_att(varargin{1}.global, 'solar_constant'); 


[~, indDate] = unique(date(:,2:end),'rows');
date = date(indDate,:);

%Calculate Julian date:
if regexpbl(sMeta.dt,{'day','daily'}) %Calculate for day:
    if numel(date(1,:)) == 3
        dayJ = Julian_day(date,'ofyear');
        datersdt = date;
    else
        error('toa_rad_DeWalle:dateNumel',['The date vector has ' num2str(numel(date(1,:))) ' elements but 3 are needed.']);
    end
elseif regexpbl(sMeta.dt,'month') %if month time-step, calculate PET for each day and average to monthly value
    dayJ = nan(31*numel(date(:,1)));
    datersdt = nan([31*numel(date(:,1)),3]);
    cntr = 0;
    for ii = 1 : numel(date(:,1))
        for jj = 1 : eomday(date(ii,1),date(ii,2))
            cntr = cntr + 1;
            datersdt(cntr,:) = [date(ii,1:2),jj];
            dayJ(cntr) = Julian_day(datersdt(cntr,:),'ofyear');
        end
    end
    dayJ = dayJ(1:cntr);
    datersdt = datersdt(1:cntr,:);
else
    error('toa_rad_DeWalle:tUnit',['Current time units are ' sMeta.dt ' but days or months are required.  Program for this.'])
end


if ~isfield(sHydro,'latGridC')
    indFlat = find(sHydro.aspectFdr == -1); %-1 is value assigned to flat surfaces.
    %Correct latitude grid to account for slope and aspect:
        %Use negative of slope angle because equation calls for angle of
        %incline (see example on pg. 399) 
    sHydro.latGridC = asind(sind(-sHydro.slopeFdrA).*cosd(sHydro.aspectFdr).*cosd(sHydro.latGrid) + cosd(-sHydro.slopeFdrA).*sind(sHydro.latGrid));
    indLatP = find(isnan(sHydro.latGridC));
    %Set cells that are nan and cells that are flat to the original latitude
    sHydro.latGridC(indLatP) = sHydro.latGrid(indLatP);
    sHydro.latGridC(indFlat) = sHydro.latGrid(indFlat);
end


%Declination angle:
%See Holbert_ASU-solarCalcs.pdf in SETI folder
angleDec = 23.45*sind((dayJ+284)*360/365);

%Eccentricity of Earth's orbit:
%radVec = (1 - e^2)/(1+e*cosd(theta))
radVec = 1 + 0.033*cosd(dayJ*360/365);

%Difference in longitude between sloped surface and equivalent horizontal
%surface:
lonDiff = atand((sind(sHydro.aspectFdr).*sind(-sHydro.slopeFdrA))./(cosd(-sHydro.slopeFdrA).*cosd(sHydro.latGrid) - cosd(sHydro.aspectFdr).*sind(-sHydro.slopeFdrA).*sind(sHydro.latGrid)));
lonDiff(isnan(lonDiff)) = 0;


%Calculate daily radiation:
tempRaad = nan([numel(angleDec),size(sHydro.latGridC)]);
for ii = 1 : numel(angleDec)
    hrs = daylight_hrs(sHydro,dayJ(ii)); %Accounts for slope and aspect
    
%     hourAngle = acosd(-tand(sHydro.latGrid)*tand(angleDec(ii))); %This is hour angle of an equivalent horizontal surface
    hourAngle1 = -7.5*hrs + lonDiff;
    hourAngle2 = 7.5*hrs + lonDiff;
    
    %Integrated cos(Z)
    cZenith = hrs.*sind(sHydro.latGridC).*sind(angleDec(ii)) ...
        + (1/15)*cosd(sHydro.latGridC).*cosd(angleDec(ii)).*(sind(hourAngle2)-sind(hourAngle1)); %Unitless
    %Daily Radiation (Units of Watts)
    tempRaad(ii,:,:) = (I0./radVec(ii))*cZenith/24;
end

%Set negative values to zero:
tempRaad(tempRaad < 0) = 0;


%Aggregate to desired time-step:
if regexpbl(sMeta.dt,{'day','daily'})
    if numel(datersdt(1,:)) == 3
        sAtm.datersdt = datersdt;
        sAtm.rsdt = tempRaad;
    else
        error('toa_rad_DeWalle:tUnit',['Current time units are ' ...
            sMeta.dt ' but date resolution is ' ...
            num2str(numel(sAtm.datersdt(1,:))) ', which is not expected.']);
    end
elseif regexpbl(sMeta.dt,'month')
    if numel(datersdt(1,:)) == 3
        sAtm.datersdt = date;
        sAtm.rsdt = nan([numel(date(:,1)), size(sHydro.dem)],'single');
        for ii = 1 : numel(date(:,1))
            sAtm.rsdt(ii,:,:) = squeeze(nanmean(tempRaad(ismember(datersdt(:,1:2), date(ii,1:2),'rows'),:,:),1)); 
        end
    elseif numel(datersdt(1,:)) == 2
        sAtm.datersdt = datersdt;
        sAtm.rsdt = tempRaad;
    else
        error('toa_rad_DeWalle:tUnit',['Current time units are ' ...
            sMeta.dt ' but date resolution is ' ...
            num2str(numel(sAtm.datersdt(1,:))) ', which is not expected.']);
    end
end


%OLD METHOD:
% 
% %Calculate potential solar radiation:
% if regexpbl(sMeta.dt,'month')
%     tempRaad = nan([numel(angleDec),size(sHydro.latGridC)]);
%     for ii = 1 : numel(angleDec)
%         %Integrate solar radiation over each day, then divide by hours in day
%         hlfHrs = acosd(-tand(sHydro.latGrid)*tand(angleDec(ii)))/15;
%         hourAngle = 15*hlfHrs+lonDiff;
%         cZenith = (2*hlfHrs.*sind(sHydro.latGridC).*sind(angleDec(ii)) + 3.8197*cosd(sHydro.latGridC).*cosd(angleDec(ii)).*(sind(hourAngle)-sind(-hourAngle)))/24; %Unitless
%         tempRaad(ii,:,:) = (I0./radVec(ii))*cZenith; %Potential Solar irradiation on flat surface
%     end
%     sAtm.rsdt = squeeze(nanmean(tempRaad,1)); %Use average radiation (this only relevent if monthly time-step used.
% elseif regexpbl(sMeta.dt,{'day','daily'})
%     hlfHrs = acosd(-tand(sHydro.latGrid)*tand(angleDec))/15;
%     hourAngle = 15*hlfHrs+lonDiff;
%     cZenith = (2*hlfHrs.*sind(sHydro.latGridC).*sind(angleDec(ii)) + 3.8197*cosd(sHydro.latGridC).*cosd(angleDec(ii)).*(sind(hourAngle)-sind(-hourAngle)))/24; %Unitless
%     sAtm.rsdt = (I0./radVec(ii))*cZenith; %Potential Solar irradiation on flat surface
% elseif regexpbl(sMeta.dt,{'hour'})
%     %Ref: DeWalle and Rango (2008)
%     %cos(Zenith) ~ sin(latitute)*sin(solar declination) + cos(latitude)*cos(solar declination)*cos("hour angle")
%     cZenith = sind(sHydro.latGridC).*sind(angleDec) + cosd(sHydro.latGridC).*cosd(angleDec).*cosd(15*(sMeta.currTime(4) - 12)+lonDiff);
%     sAtm.rsdt = (I0./radVec).*cZenith; %Potential Solar irradiation on flat surface
% else
%     error('toa_rad_DeWalle:dt',['The specificied time step is ' char(39) sMeta.dt char(39) ', which has not been coded for.']);
% end
% 
% %Scale by fitting parameter:
% sAtm.rsdt = (10^a)*sAtm.rsdt;
