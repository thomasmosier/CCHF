function [varLd, varDisp] = clm_var_load(modules)


%Remove tasmin and tasmax if 'simple' version of cryosphere model being run.    
if regexpbl(find_att(modules,'heat'), {'simple','SDI'})
    varLd = {'pr','tas'};
    varDisp = {'precipitation', 'near-surface air temperature'};
elseif regexpbl(find_att(modules,'heat'), {'Pelli', 'Hock', 'ETI', 'LST', 'SETI'})
    varLd = {'pr','tas', 'tasmax', 'tasmin'};
    varDisp = {'precipitation', ...
        'mean near-surface air temperature', ...
        'maximum near-surface air temperature', ...
        'minimum near-surface air temperature'};
elseif regexpbl(find_att(modules,'heat'), {'Pritchard'})
    varLd = {'pr','tas', 'tasmax', 'tasmin', 'hurs', 'sfcwind'};
    varDisp = {'precipitation', ...
        'mean near-surface air temperature', ...
        'maximum near-surface air temperature', ...
        'minimum near-surface air temperature', ...
        'relative humidity'...
        'wind speed', ...
        };
else
    error('clm_var_ld:unknownRepresentation',...
        ['The heat process representation ' find_att(modules,'heat') ...
        ' has not been programmed for.']);
end

%Some variable abbreviation that may be of interest:
%(see 'CMIP5 Standard Output' for more field abbreviations)
%'pr' = total precipitation
%'tas' = near-surface air temperature
%'tasmax' = daily minimum near-surface air temperature
%'tasmin' = daily maximum near-surface air temperature
%'hurs' = near-surface relative humidity
%'huss' = near-surface specific humidity
%'uas' = eastward near-surface wind
%'vas' = northward near-surface wind
%'sfcwind' = near-sufrace wind speed