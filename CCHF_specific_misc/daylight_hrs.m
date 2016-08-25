function hrs = daylight_hrs(sHydro,dayJ)
%Hours of daylight based on latitude and day of year

%Formulation comes from
%DeWalle, D. R., & Rango, A. (2008). Principles of snow hydrology. 
%Cambridge University Press. (Appendix B)

%Solar Zenith Angle = 90 at Sunrise and sunset

% %cosine of day of year, including Julian Day correction:
% cTheta = cosd(0.98565*(dayJ+10)+1.914*sind(0.98565*(dayJ-2)));
% %Solar Declination Angle:
% angleDec = -asind(0.39779*cTheta);

%Declination angle:
%See Holbert_ASU-solarCalcs.pdf in SETI folder
angleDec = 23.45*sind((dayJ+284)*360/365);

if ~isfield(sHydro,'latGridC')
    indFlat = find(sHydro.aspectFdr == -1); %-1 is value assigned to flat surfaces.
    %Correct latitude grid to account for slope and aspect:
    sHydro.latGridC = asind(sind(-sHydro.slopeFdrA).*cosd(sHydro.aspectFdr).*cosd(sHydro.latGrid) + cosd(-sHydro.slopeFdrA).*sind(sHydro.latGrid));
    indLatP = find(isnan(sHydro.latGridC));
    %Set cells that are nan and cells that are flat to the original latitude
    sHydro.latGridC(indLatP) = sHydro.latGrid(indLatP);
    sHydro.latGridC(indFlat) = sHydro.latGrid(indFlat);
end


% %Difference in longitude between sloped surface and equivalent horizontal
% %surface:
% lonDiff = atand((sind(sHydro.aspectFdr).*sind(sHydro.slopeFdrA))./(cosd(sHydro.slopeFdrA).*cosd(sHydro.latGrid) - cosd(sHydro.aspectFdr).*sind(sHydro.slopeFdrA).*sind(sHydro.latGrid)));
% lonDiff(isnan(lonDiff)) = 0;

hrs = nan([numel(dayJ), size(sHydro.latGridC)],'single');
for ii = 1 : numel(dayJ)
    hrsFlat = real(2*acosd(-tand(sHydro.latGrid)*tand(angleDec(ii)))/15);
    hrsSlope = real(2*acosd(-tand(sHydro.latGridC)*tand(angleDec(ii)))/15);
    hrsSlope(hrsSlope > hrsFlat) = hrsFlat(hrsSlope > hrsFlat);
    hrs(ii,:,:) = hrsSlope;
end

%Correct for unacceptable values:
hrs(hrs > 24) = 24;
hrs(hrs <  0) =  0;

hrs = squeeze(hrs);

%Old formulation, which didn't work well:
% %http://mathforum.org/library/drmath/view/56478.html
% %Formula from: Ecological Modeling, volume 80 (1995) pp. 87-95, "A Model 
% %Comparison for Daylength as a Function of Latitude and Day of the 
% %Year."
% 
% P = asin(0.39795*cos(0.2163108 + 2*atan(.9671396*tan(0.00860*(dayJ-186)))));
% 
% D = 24 - (24/pi)*acos((sind(0.8333) + sind(lat)*sin(P))./(cosd(lat)*cos(P)));