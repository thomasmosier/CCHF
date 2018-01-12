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