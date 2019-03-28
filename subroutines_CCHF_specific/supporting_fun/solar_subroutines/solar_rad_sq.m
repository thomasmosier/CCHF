function e2 = solar_rad_sq(dayJ)

%Radius vector (variation) of Earth's orbit
%radius = a*(1 - e^2)/(1+e*cosd(theta)), where a is semi-major axis, e
%is eccentricity, and theta is angle from perihelion
%Theta = (Julian date)*360/365.25
%eccentricity varies between 0 and 0.06. Currently e ~0.0167

%For testing, compare function to values from Table copied below:
% dayJ = [173; 193; 208; 222; 237; 252; 266; 281; 295; 309; 323; 337; 356];
% e2Table = [1.03297; 1.03090; 1.02728; 1.02190; 1.01528; 1.00739; 0.99960; 0.99154; 0.98404; 0.97790; 0.97285; 0.96938; 0.96759];

e2 = 1 + 0.03297*cos(2*pi*(dayJ-174)/365.25);
%Note: Solstice offset of 174 empirically obtained through minimizing
%the mean absolute percentage difference with values from table below.
% mean(abs(100*(e2 - e2Table)./e2Table))

%[e2, e2Table, 100*(e2 - e2Table)./e2Table]
%mean(abs(100*(e2 - e2Table)./e2Table))

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

%http://curious.astro.cornell.edu/about-us/41-our-solar-system/the-earth/orbit/80-how-can-i-find-the-distance-to-the-sun-on-any-given-day-advanced