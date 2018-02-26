function trans = transmissivity_DeWalle(tasRng, a, b, c)

%Eq. 6.7 on pg. 151 of DeWalle and Rango 2008

trans = b*(1-exp(-0.01*a*tasRng.^c));

trans(trans > 1) = 1;
trans(trans < 0) = 0; %This line shouldn't be necessary