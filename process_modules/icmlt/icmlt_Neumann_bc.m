function varargout = icmlt_Neumann_bc(varargin)
%%MELT ICE:
    
global sCryo


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'tmp_depth', -2, 2, 0, 'icmlt_Neumann_bc','cryo'});
    
    return
else
	tmpDepth = find_att(varargin{1}.coef,'tmp_depth');
    sMeta = varargin{1};
end


kGlac = find_att(sMeta.global,'thermal_conduct_ice'); 
densW = find_att(sMeta.global,'density_water'); 
cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3


%Calculate conduction from bottom boundary condition of glacier upwards:
%conduct = -k * temp grad / distance
sCryo.hficbc = -kGlac*13/tmpDepth;

sCryo.lhicme  = zeros(size(sCryo.snw),'single');
sCryo.icdwe   = zeros(size(sCryo.snw),'single');
sCryo.iclr    = zeros(size(sCryo.snw),'single'); %Ice liquid release


%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'icx')
    eps = 0.001;
    indIceMlt = find( sCryo.icx ~= 0 & sCryo.hfneti > 0 & sCryo.snw < eps & sCryo.lhsnme < eps);
else
    indIceMlt = [];
end

%Calculate changes using residual energy at locations with ice:
if ~isempty(indIceMlt)
    %Calculate ice melt, which involves different degree-index factor

    sCryo.lhicme(indIceMlt) = time2sec(1,sMeta.dt,sMeta.dateCurr)*sCryo.hfneti(indIceMlt)/cLate;
        sCryo.lhpme(indIceMlt) = 0;
        
    sCryo.lhicme( sCryo.lhicme < 0) = 0;
    sCryo.iclr(indIceMlt) = sCryo.icx(indIceMlt).*sCryo.lhicme(indIceMlt);
    sCryo.icdwe(indIceMlt) = -sCryo.lhicme(indIceMlt);
end



% %Set ice Temp (locations where no snow):
% sCryo.tic(sCryo.icx ~= 0 & sCryo.snw == 0) = 0;