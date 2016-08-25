function [gWE, gThick] = glacier_WE_init(glacierBl, slope, sMeta)
%output is ice thickness (in terms of water equivalent), assuming glacier-
%climate equilibrium.

%Equilibrium glacier thickness (Eq. 8) of
%Shea, J. M., Immerzeel, W. W., Wagnon, P., Vincent, C., & Bajracharya, S. 
%(2015). Modelling glacier change in the Everest region, Nepal Himalaya. 
%The Cryosphere, 9(3), 1105–1128. http://doi.org/10.5194/tc-9-1105-2015

tauNaught = find_att(sMeta.global,'equil_shear'); 
rhoI = find_att(sMeta.global,'density_ice'); 
rhoW = find_att(sMeta.global,'density_water');
g = find_att(sMeta.global,'grav_accel');

%Set minimal angle to 1.5 degrees (sind(1.5) = 0.0262) (From Shea et al.)

slope(slope > -0.0262 & slope < 0) = -0.0262; 
slope(slope >= 0 & slope <= 0.0262) = 0.0262;

gThick = abs(tauNaught./(-rhoI*g*slope));
gThick(gThick < 0) = 0;
gThick(isnan(glacierBl)) = nan;
gThick(glacierBl == 0) = 0;

if issparse(glacierBl)
    gThick = sparse(double(gThick));
else
    gThick = single(gThick);
end

gWE = (rhoI/rhoW)*gThick;