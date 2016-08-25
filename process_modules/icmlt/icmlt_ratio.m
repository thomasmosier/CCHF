function varargout = icmlt_ratio(varargin)
%%MELT ICE:
    
global sCryo


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr' , -2, 2, -0.2, 'icmlt_ratio','cryo'});  %Units of mm/hr/Deg C
    
    if isfield(sCryo,'icdbr')
        varargout{1} = cat(1,varargout{1}, {'debris_scalar' , 0, 2, 1, 'icmlt_ratio','cryo'});  %Unitless
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
sCryo.icdwe       = zeros(size(sCryo.snw),'single');
sCryo.iclr        = zeros(size(sCryo.snw),'single'); %Ice liquid release


%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'icx')
    indIceMlt = find( sCryo.icx ~= 0 & sCryo.lhpme > 0);
else
    indIceMlt = [];
end

%Calculate changes using residual energy at locations with ice:
if ~isempty(indIceMlt)
    %Calculate ice melt, which involves different degree-index factor
    ratioMelt = sCryo.hfneti(indIceMlt)./sCryo.hfnet(indIceMlt);
    
    %Ensure ratio of melt is positive:
    %Set to zero if negative because ice heat flux is zero:
    ratioMelt(sCryo.hfneti(indIceMlt) < 0) = 0;
    %If negative because snow surface heat flux is zero, switch sign:
    ratioMelt(ratioMelt < 0) = -ratioMelt(ratioMelt < 0);
    %Set maximum scaling factor:
    maxRat = 4;
    ratioMelt(ratioMelt > maxRat) = maxRat;
    
    %Calculate energy available for ice melt:
    sCryo.lhicme(indIceMlt) = 10^(dii)*sCryo.lhpme(indIceMlt).*ratioMelt;
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


% %Set ice Temp (locations where no snow):
% indBIce = find(sCryo.icx ~= 0 & sCryo.snw == 0);
%     indBIce3 = ind2_3d(size(sAtm.tas),indBIce, sAtm.indCurr);
% %     sCryo.tic(indBIce) = 0;
% %     sCryo.tic(indBIce) = min(0, squeeze(sAtm.tas(indBIce3)));
% %Temprature of bare ice is min(0, air temp):
% if isfield(sAtm,'tasmin')
%     sCryo.tic(indBIce) = min(0, squeeze(sAtm.tasmin(indBIce3)));
% else
%     sCryo.tic(indBIce) = min(0, squeeze(sAtm.tas(indBIce3)));
% end