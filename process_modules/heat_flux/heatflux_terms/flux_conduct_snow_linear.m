function argout = flux_conduct_snow_linear(tsis, varargin)

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
    
    if numel(varargin(:)) > 1
       mode = varargin{2};
    else
        mode = '';
    end
end

%Set maximum value for conduction:
condtMax = 500; %Units = Watts per meter squared


if ~isfield(sCryo,'rhosn')
    sCryo.rhosn = 400*ones(size(sCryo.snw),'single');
end

if ~regexpbl(mode,'deriv') %Normal mode:
    if isfield(sCryo, 'tsn')
        argout = -10^(c).*(tsis - sCryo.tsn) ./ sCryo.snw;
    else
        error('heat_conduct_snow:missingTmp','Surface or internal snow temperature not defined.');
    end
else %Derivative of above:
    argout = -10^(c) ./ sCryo.snw;
end

%Ensure no odd values (due to there not being any snow:
argout(isnan(argout)) = 0;
argout(sCryo.snw == 0) = 0;
argout(argout > condtMax) = condtMax;
argout(argout < -condtMax) = -condtMax;

%If temperatures equal, then no conduction
argout(sCryo.tsn == tsis) = 0;