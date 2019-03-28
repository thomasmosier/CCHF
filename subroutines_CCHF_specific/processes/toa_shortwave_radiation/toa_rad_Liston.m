% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function varargout = toa_rad_Liston(date, sHydro, varargin)

global sAtm

%REFERENCES:
%Liston's MicroMet routine ('solar_rad' subroutine):
%Liston, G. E., & Elder, K. (2006). A meteorological distribution system 
    %for high-resolution terrestrial modeling (MicroMet). Journal of 
    %Hydrometeorology, 7(2), 217–234.
%Garnier, B. J., & Ohmura, A. (1968). A method of calculating the direct 
    %shortwave radiation income of slopes. Journal of Applied Meteorology, 
    %7(5), 796–800.
%Coops, N. C., Waring, R. H., & Moncrieff, J. B. (2000). Estimating mean 
    %monthly incident solar radiation on horizontal and inclined slopes from 
    %mean monthly temperatures extremes. International Journal of 
    %Biometeorology, 44(4), 204–211.


if isempty(varargin(:))
	varargout{1} = cell(0,6);
 
    return
else
    sMeta = varargin{1};
end
    



%Solar Constant (W/m^2):
I0 =  find_att(sMeta.global, 'solar_constant'); 


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
    error('toa_rad_DeWalle:monthTs',[char(39) 'toa_rad_Coops' char(39) ...
        ' can only be run for daily time-steps.']);
else
    error('toa_rad_DeWalle:tUnit',['Current time units are ' sMeta.dt ' but days or months are required.  Program for this.'])
end

latCan = find_att(varargin{1}.global, 'lat_cancer'); 

%Solar declination (degrees)
    %23.45 = tropic of cancer
solDec = latCan*sind((dayJ+284)*360/365.25);
%From Liston (does not agree):
    %latCan*sin(2*pi*(dayJ-173)/365);
%From internet (does agree):
    %asind(sind(latCan)*sind(360/365*(dayJ-81)))

%Consider aspect of each grid cell:

%Create southern exposure azimuth grid:
slopeAzS = nan(size(sHydro.latGrid),'single');
indN = find(sHydro.latGrid > 0);
if ~isempty(indN)
    indW = find(sHydro.aspectFdr(indN) >= 180);
    slopeAzS(indN(indW)) = sHydro.aspectFdr(indN(indW)) - 180;
    indE = setdiff(indN,indN(indW));
    slopeAzS(indN(indE)) = sHydro.aspectFdr(indN(indE)) + 180;
end

indS = setdiff((1:numel(sHydro.latGrid)),indN);
slopeAzS(indS) = sHydro.aspectFdr(indS);


%Initialize output grid:
sAtm.rsdt = nan([numel(dayJ),size(sHydro.latGrid)],'single');

for jj = 1 : numel(dayJ)
    %Loop over hour (this is all from Liston):
    hr = (1:24);
    
    toaRadDirTemp = nan([numel(hr), size(sHydro.latGrid)],'single');
    toaRadDifTemp = nan([numel(hr), size(sHydro.latGrid)],'single');
    for ii = 1 : 24 
        %Hour angle (degrees)
        hrAng = hr(ii)*15-180;

        %Solar zenith angle (from Garnier and Ohmura, 1968):
        cosZ = sind(solDec(jj))*sind(sHydro.latGrid) + cosd(solDec(jj))*cosd(sHydro.latGrid)*cosd(hrAng);
        cosZ(cosZ < 0) = 0;

        sinZ = sqrt(1-cosZ.*cosZ);

        %Calculate azimuth angle of sun (degrees):
        sunAz = cosd(solDec(jj))*sind(hrAng)./sinZ;
        %Ensure within limits
            sunAz(sunAz > 1) = 1;
            sunAz(sunAz < -1) = -1;
        sunAz = asind(sunAz);

        %Flip solar azimuth if in Southern hemisphere:
        indSouth = find(sHydro.latGrid < 0);
        if ~isempty(indSouth)
            sunAz(indSouth) = -sunAz(indSouth);
        end

        %Check position of sun relative to slope:
        if hrAng < 0
            indX = find(hrAng < sunAz);
            sunAz(indX) = -180 - sunAz(indX);
        elseif hrAng > 0
            indX = find(hrAng > sunAz);
            sunAz(indX) = 180 - sunAz(indX);
        end
        
        %Calculate the angle of incidence (degrees):
        cosI = cosd(sHydro.slopeFdrA).*cosZ + sind(sHydro.slopeFdrA).*sinZ.*cosd(sunAz - slopeAzS);
        cosI(cosI < 0) = 0;
        cosI(cosZ <= 0) = 0;
        
        toaRadDirTemp(ii,:,:) = I0*cosI.*(0.6+0.2*cosZ);
        toaRadDifTemp(ii,:,:) = I0*cosZ.*(0.3+0.1*cosZ);
    end
    
    sAtm.rsdt(jj,:,:) = squeeze(nanmean(toaRadDirTemp,1));
    sAtm.rsdf(jj,:,:) = squeeze(nanmean(toaRadDifTemp,1));
end

sAtm.datersdt = datersdt;












% if ~isfield(sHydro,'latGridC')
%     indFlat = find(sHydro.aspectFdr == -1); %-1 is value assigned to flat surfaces.
%     %Correct latitude grid to account for slope and aspect:
%     sHydro.latGridC = asind(sind(sHydro.slopeFdrA).*cosd(sHydro.aspectFdr).*cosd(sHydro.latGrid) + cosd(sHydro.slopeFdrA).*sind(sHydro.latGrid));
%     indLatP = find(isnan(sHydro.latGridC));
%     %Set cells that are nan and cells that are flat to the original latitude
%     sHydro.latGridC(indLatP) = sHydro.latGrid(indLatP);
%     sHydro.latGridC(indFlat) = sHydro.latGrid(indFlat);
% end
% 
% %Declination angle:
% %See Holbert_ASU-solarCalcs.pdf in SETI folder
% solDec = 23.45*sind((dayJ+284)*360/365);
% 
% %Eccentricity of Earth's orbit:
% %radVec = (1 - e^2)/(1+e*cosd(theta))
% radVec = 1 + 0.033*cosd(dayJ*360/365);
% 
% %Difference in longitude between sloped surface and equivalent horizontal
% %surface:
% lonDiff = atand((sind(sHydro.aspectFdr).*sind(sHydro.slopeFdrA))./(cosd(sHydro.slopeFdrA).*cosd(sHydro.latGrid) - cosd(sHydro.aspectFdr).*sind(sHydro.slopeFdrA).*sind(sHydro.latGrid)));
% lonDiff(isnan(lonDiff)) = 0;
% 
% 
% 
% % D = asind(0.39785*sind(278.9709 + 0.9856*JDay) + 1.9163*sind(356.6153+ 0.9856*jDay));
% 
% %cosine of day of year, including Julian Day correction:
% cTheta = cosd(0.98565*(dayJ+10)+1.914*sind(0.98565*(dayJ-2)));
% 
% %Solar Declination Angle:
% solDec = -asind(0.39779*cTheta);
% 
% %Seconds of Sun:
% dtSun = 3600*daylight_hrs(sHydro,dayJ);
% 
% 
% %Calculate daily radiation:
% tempRaad = nan([numel(solDec),size(sHydro.latGridC)]);
% for ii = 1 : numel(solDec)
%     hourAngle1 = -7.5*hrs + lonDiff;
%     hourAngle2 = 7.5*hrs + lonDiff;
% 
%     %Integrated cos(Z)
%     Z = acosd(sind(sHydro.latGridC).*sind(solDec(ii)) ...
%         + cosd(sHydro.latGridC).*cosd(solDec(ii)).*(sind(hourAngle2)-sind(hourAngle1))); 
% 
%     %Optical Air Mass:
%     M = 1/(cosd(Z)) + 1*10^(-7); 
% %     %z is the zenith angle and must be modified by the angle of incidence:
%     %Only modify for sub-daily calculations:
% %     cos(i) = cos(sHydro.slopeFdr)*cos(Z)+sin(sHydro.slopeFdr)*sin(Z)*cos(solarAzimuth + sHydro.aspectFdr);
% 
%     %I is direct beam solar radiation:
%     I = cosd(sHydro.latGrid)*(I0*dtSun*(sAtm.tasrange^M));
% 
%     %Diffuse radiation:
%     dF = ((((I0*dtSun*cos(Z))^2)*(sAtm.rstran^M))^0.5)*(1-I0*dtSun*(sAtm.rstran^M)*cos(Z))^0.5;
%     dSlope = dF*cos(sHydro.slopeFdr/2)^2;  
%     %Integrate above in 20 min intervals:
% 
%     sAtm.rsdt = dSlope + I;
% end
% for ii = 1 : 
%     sCryo.heatSW = sCryo.heatSW + D + I;
% 
% end
% 
% sAtm.datersdt = datersdt;
% sAtm.rsdt = tempRaad;