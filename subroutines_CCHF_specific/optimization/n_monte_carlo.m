function nRCarlo = n_monte_carlo(nPrm)

%Determine number of Monte Carlo runs (based on number of parameters)
if nPrm <= 4
    nRCarlo = 500;
elseif nPrm <= 7
    nRCarlo = 1000;
elseif nPrm > 7 && nPrm < 10
    nRCarlo = 2000;
else
    nRCarlo = 3000;
end