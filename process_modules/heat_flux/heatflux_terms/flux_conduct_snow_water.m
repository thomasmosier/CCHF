function argout = flux_conduct_snow_water(tsis, varargin)

%Positive values indicate heat flowing from internal snowpack to snow-air interface
%Negative values indicate heat flowing from snow-air interface to internal snowpack


global sCryo

if isempty(varargin(:))
	argout = cell(0,6);
    argout = cat(1, argout, {'cond_pwr', -2, 2, 1, 'heat_conduction_snow','cryo'});  
    argout = cat(1,argout, ['tsn',cell(1,5)]);
    return
else
    c = find_att(varargin{1}.coef,'cond_pwr');
    sMeta = varargin{1};
    
    if numel(varargin(:)) > 1
       mode = varargin{2};
    else
        mode = '';
    end
end

argout = nan(size(sCryo.snw));

%Set maximum value for conduction:
condtMax = 500; %Units = Watts per meter squared

%Find global attributes:
rhoDry = find_att(sMeta.global, 'density_snow_dry');
rhoWet = find_att(sMeta.global, 'density_snow_wet');
rhoWater = find_att(sMeta.global, 'density_water');
kDry = find_att(sMeta.global, 'thermal_conduct_snow_dry');
kWet = find_att(sMeta.global, 'thermal_conduct_snow_wet'); 

%Find indices of dry and wet snow:
indWet = find(sCryo.snlq > 0.01 * sCryo.snw); %this value somewhat arbitrary
indDry = setdiff((1:numel(sCryo.snw)), indWet);
        

if ~regexpbl(mode,'deriv') %Normal mode:
    if isfield(sCryo, 'tsn')
        argout(indWet) = -10^(c)*(2*kWet*rhoWet/rhoWater) ...
            *((tsis(indWet) - sCryo.tsn(indWet))./ sCryo.snw(indWet));
        argout(indDry) = -10^(c)*(2*kDry*rhoDry/rhoWater) ...
            *((tsis(indDry) - sCryo.tsn(indDry))./ sCryo.snw(indDry));
    else
        error('heat_conduct_snow:missingTmp','Surface or internal snow temperature not defined.');
    end
else %Derivative of above:
    argout(indWet) = -10^(c)*(2*kWet*rhoWet/rhoWater) ./ sCryo.snw(indWet);
    argout(indDry) = -10^(c)*(2*kDry*rhoDry/rhoWater) ./ sCryo.snw(indDry);
end

%Ensure no odd values (due to there not being any snow:
argout(isnan(argout)) = 0;
argout(sCryo.snw == 0) = 0;
argout(argout > condtMax) = condtMax;
argout(argout < -condtMax) = -condtMax;

%If temperatures equal, then no conduction
argout(sCryo.tsn == tsis) = 0;