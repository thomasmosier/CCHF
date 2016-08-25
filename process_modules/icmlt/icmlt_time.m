function varargout = icmlt_time(varargin)
%%MELT ICE:
    
global sCryo


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr' , -2, 2, -0.2, 'icmlt_remainder','cryo'});  %Units of mm/hr/Deg C
    
    if isfield(sCryo,'icdbr')
        varargout{1} = cat(1,varargout{1}, {'debris_scalar' , 0, 2, 1, 'icmlt_remainder','cryo'});  %Unitless
    end
    
    return
else
    dii = find_att(varargin{1}.coef,'ice_melt_pwr');  %Units of mm/hr/Deg C
    
%     if isfield(sCryo,'tsn')
%         tsnPwr = find_att(varargin{1}.coef,'tsn_pwr');
%     end
    
    if isfield(sCryo,'icdbr')
        debScl = find_att(varargin{1}.coef,'debris_scalar');  %Unitless
    end
    
    sMeta = varargin{1};
end

if ~isfield(sCryo, tic)
    sCryo.tic = nan(size(sCryo.snw),'single');
end

sCryo.lhicme  = zeros(size(sCryo.snw),'single');
sCryo.icdwe   = zeros(size(sCryo.snw),'single');
sCryo.iclr    = zeros(size(sCryo.snw),'single'); %Ice liquid release


%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'icx')
    indIceMlt = find( sCryo.icx ~= 0 & sCryo.hfneti > 0 & sCryo.snw < eps & sCryo.lhsnme < eps);
else
    indIceMlt = [];
end

%Calculate changes using residual energy at locations with ice:
if ~isempty(indIceMlt)
    indNMlt = setdiff((1:numel(sCryo.lhsnme)), indIceMlt);
    
    densW = find_att(sMeta.global,'density_water'); 
    cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3
    
    %Calculate how much time was required to melt snow:
    timeSn = cLate*sCryo.lhsnme./sCryo.hfnet;
    timeSn(indNMlt) = 0;
    
    tStep = time2sec(1,sMeta.dt,sMeta.dateCurr);
    
    timeIc = tStep - timeSn;
    timeIc(timeIc < 0 ) = 0;
    timeIc(timeIc > tStep ) = tStep;
    
    sCryo.lhicme(indIceMlt) = 10^(dii)*sCryo.hfneti(indIceMlt).*timeIc(indIceMlt)/cLate;
    sCryo.lhicme(sCryo.lhicme < 0) = 0;
    sCryo.lhpme(indIceMlt) = 0;
        
    %Scale melt for debris covered glacier cells:
    if isfield(sCryo,'icdbr')
        indDebrisMlt = intersect(indIceMlt, find(~isnan(sCryo.icdbr))); 
    	sCryo.lhicme(indDebrisMlt) = debScl*sCryo.lhicme(indDebrisMlt);
    end
    
    sCryo.lhicme( isnan(sCryo.lhicme)) = 0;
    %*'iclr' accounts for fact that glacier coverage of each cell may only
    %be partial
    sCryo.iclr(indIceMlt) = sCryo.icx(indIceMlt).*sCryo.lhicme(indIceMlt);
    %*'icdwe' is depth change per unit area of ice (indpenedent of 
    %fractional coverage)
    sCryo.icdwe(indIceMlt) = -sCryo.lhicme(indIceMlt);
end