function angleDec = solar_dec_angle(dayJ)
%Declination angle:

%Method from NOAA (Most Accurate here; accounts for eliptical orbits):
%Website: https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
yrFrac = ((2*pi)/365)*(dayJ-1);
angleDec = (180/pi)*(0.006918 ...
    - 0.399912*cos(  yrFrac) + 0.070257*sin(  yrFrac) ...
    - 0.006578*cos(2*yrFrac) + 0.000907*sin(2*yrFrac) ...
    - 0.002697*cos(3*yrFrac) + 0.001480*sin(3*yrFrac));


% %Another Method that accounts for elliptical orbit:
% %See pdf by REUK (or website: http://www.reuk.co.uk/wordpress/solar/solar-declination/)
% a = 360/365.24;
% b = (360*0.0167)/pi;
% angleDec = asind(sind(-23.44)*cosd(a*(dayJ+10) + b*sind(a*(dayJ - 2))));

% %Method From Liston's Micromet:
% tropCanc = 23.5;
% solst = 173;
% daysYr = 365.25;
% 
% angleDec = tropCanc*cos(2*pi*(dayJ-solst)/daysYr);

%See Holbert_ASU-solarCalcs.pdf distributed with CCHF
% angleDec = 23.45*sind((dayJ+284)*360/365);

% angleDec = asind(0.39795*cosd(0.98563*(dayJ - 173)));

%TABLE FROM DeWalle and Rango (pg. 395)
%Date,              Julian Day,     Solar declination,  (radius vector)^2
%6/22 (solstice),   173,            +23.5,              1.03297
%6/1 and 7/12,      152 and 193,    +21.967,            1.03090
%5/18 and 7/27,     138 and 208,    +19.333,            1.02728
%5/3 and 8/10,      123 and 222,    +15.583,            1.02190
%4/19 and 8/25,     109 and 237,    +10.917,            1.01528
%4/4 and 9/9,       94 and 252,     +5.633,             1.00739
%3/21 and 9/23,     80 and 266,     +0,                 0.99960
%3/7 and 10/8,      66 and 281,     -5.633,             0.99154
%2/20 and 10/22,    51 and 295,     -10.917,            0.98404
%2/7 and 11/5,      38 and 309,     -15.583,            0.97790
%1/24 and 11/19,    24 and 323,     -19.333,            0.97285
%1/10 and 12/3,     10 and 337,     -21.967,            0.96938
%12/22 (solstice),  356,            -23.5,              0.96759