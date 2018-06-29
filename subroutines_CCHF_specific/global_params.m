function output = global_params
%This function is a repository of all constants used throughout the script
%(Not currently complete. Continue to gather constants from functions here.)

%retrieve constants using function 'value = find_att(listAtts, strAtt)'
output = ...
    {...
    'solar_constant', 1360; ... %(W/m^2) DeWalle and Rango, Principles of Snow Hydrology
    'gas_constant_univ', 8.3144; ...%(J K-1 mol-1) universal gas constant
    'albedo_ground', 0.15; ... %(unitless)
    'albedo_ice', 0.35; ... %From Hock's energy balance notes (unitless) %Liston uses ice albedo = 0.4
    'albedo_debris', 0.15; ... %(unitless)
    'albedo_water', 0.05; ... %(unitless)
    'albedo_snow_fresh', 0.9; ... %From Hock's energy balance notes (unitless)
    'albedo_snow_old', 0.55; ... %From Hock's energy balance notes (unitless)
    'heat_cap_water', 4180; ... %(J/kg/K)
    'heat_cap_ice', 2097; ... %(J/kg/K) %From Hock's energy balance notes (unitless)
    'heat_cap_snow', 2009; ... %(J/kg/K) %From Hock's energy balance notes (unitless)
    'heat_cap_debris', 750;...%(J/kg/K) %From Rounce and McKinney (2014)
    'latent_water', 334000; ... %(J/kg) 
    'density_water', 1000; ... %(kg/m^3) 
    'density_ice', 917; ... %(kg/m^3) DeWalle and Rango, Principles of Snow Hydrology
    'density_ice_Bolch', 850; ... (kg/m^3) long-term density of combined snow and ice used by Tobias Bolch for geodetic measurement uncertainty
    'density_debris', 2700; ... (kg/m^3) %From Rounce and McKinney (2014)
    'debris_threshold', 0.04;... (m) Approx from Fig S5 in Kraaijenbrink, P. D. A., Bierkens, M. F. P., Lutz, A. F., & Immerzeel, W. W. (2017). Impact of a global temperature rise of 1.5 degrees Celsius on Asia?s glaciers. Nature, 549(7671), 257?260. https://doi.org/10.1038/nature23878
    'thermal_conduct_ice', 2.1; ... %(W/m/K; at 0 Deg C) Andy Aschwanden Glacier Thermodynamics notes
    'thermal_conduct_debris', 0.96; ... (W/m/K) %From Rounce and McKinney (2014)
    'thermal_diffus_ice', 1.09*10^(-6); ... %(m^2/s; at 0 Deg C) Andy Aschwanden Glacier Thermodynamics notes
%     'density_snow_new', 400; ... %(kg/m^3)
%     'density_snow_wet', 700; ... %(kg/m^3)
%     'thermal_conduct_snow_dry', 0.26; ... %(W/m/K; density = 400) Loosely based on Sturm, M., Holmgren, J., K�nig, M., & Morris, K. (1997). The thermal conductivity of seasonal snow. Journal of Glaciology, 43(143), 26-41.
%     'thermal_conduct_snow_wet', 1.6; ... %(W/m/K; density = 700) Loosely based on Sturm, M., Holmgren, J., K�nig, M., & Morris, K. (1997). The thermal conductivity of seasonal snow. Journal of Glaciology, 43(143), 26-41.
    'lat_cancer', 23.45; ... %(decimal degrees) Latitude of the tropic of cancer
    'emm_clear', 0.7; ... %(clear sky emmissivity)
    'emm_cloud', 1; ... %(cloudy sky emmissivity)
    'equil_shear', 8*10^4; ... %(N/m^2; equilibrium sheer stress)
    'grav_accel', 9.81; ... %(m/s^2; gravitational acceleration)
    'avalanche_angle', 50; ... %(degrees; angle at which snow modelled to avalanche to downhill grid cell; value based on Shea et al., Cryosphere, 2015)
    'radius_Earth', 6371000; ... %(meters)
    'snow_temp_min', -40; ... %(deg C; minimum allowed snow temperature)
    'SCA_SWE_opt_depth', 0.05; ... %(meters; depth of SWE for grid cell to be classified as snow covered)
    'debris_opt_depth', 0.05; ... %(meters); depth of debris required for glacier to be considered optically debris-covered
    };