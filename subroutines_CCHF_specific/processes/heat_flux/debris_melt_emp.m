function relMelt = debris_melt_emp(thick)

%Define empirically derived function to modify heat based on debris cover. 
%See Jupyter notebook ('Glacier Melt Under Debris.ipynb') for explanation
%of data and derivation
%Empirical input data used for derivation is from Kraaijenbrink et al.
%(2017), Fig. S5.

relMelt = ones(size(thick), 'single');

%Two regimes: linear scaling and 1/x^2

indLin = find(thick > 0 & thick <= 0.82);
indQuad = find(thick > 0.82);

if ~isempty(indLin)
    relMelt(indLin) = 1 + 0.391341*thick(indLin);
end

if ~isempty(indQuad)
    relMelt(indQuad) = 1.3209 + 63.7515046863/((thick(indQuad)+6.46889438).^2) - 1.1765100993;
end

% %Transform from percentage to fraction and only keep real part (natural log
% %may introduce imaginary component)
% relMelt = real(relMelt)/100;