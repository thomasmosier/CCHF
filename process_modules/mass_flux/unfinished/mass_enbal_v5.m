function varargout = mass_enbal_v5(varargin)

    
global sCryo sAtm


if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'gen_melt_pwr' , -2, 2, 0, 'melt_enbal','cryo'});
    varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr' , -2, 2, -0.2, 'melt_enbal','cryo'});  %Units of mm/hr/Deg C
%     varargout{1} = cat(1,varargout{1}, {'tsn_pwr', -2,   2,    0, 'melt_enbal','cryo'});
    
    if isfield(sCryo,'icdbr')
        varargout{1} = cat(1,varargout{1}, {'debris_scalar' , 0, 2, 1, 'melt_cc','cryo'});  %Unitless
    end
    
    return
else
    dgi = find_att(varargin{1}.coef,'gen_melt_pwr');  %Units of mm/hr/Deg C
    dii = find_att(varargin{1}.coef,'ice_melt_pwr');  %Units of mm/hr/Deg C
    
%     if isfield(sCryo,'tsn')
%         tsnPwr = find_att(varargin{1}.coef,'tsn_pwr');
%     end
    
    if isfield(sCryo,'icdbr')
        debScl = find_att(varargin{1}.coef,'debris_scalar');  %Unitless
    end
    
    sMeta = varargin{1};
end


%This conversion complements a similar one in 'ETI_Pellicciotti' and
%'simple_degree'. Makes units consistent with SETI
densW = find_att(sMeta.global,'density_water'); 
cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3
cSens = densW*find_att(sMeta.global,'heat_cap_snow');
minTemp = find_att(sMeta.global,'snow_temp_min'); 

%Initialize release field:
sCryo.lhsnme = zeros(size(sCryo.snw),'single');
sCryo.lhicme  = zeros(size(sCryo.snw),'single');
sCryo.icdwe       = zeros(size(sCryo.snw),'single');
sCryo.iclr        = zeros(size(sCryo.snw),'single'); %Ice liquid release
if ~isfield(sCryo, 'lwsnl')
    sCryo.lwsnl = zeros(size(sCryo.snw),'single');
end



% if isequal(sMeta.dateCurr, [2003,5]) || isequal(sMeta.dateCurr, [2003,6]) 
%    keyboard 
% end

%Calculate melt potential using simple degree indec formulation (units of m): 
sCryo.lhpme = 10^(dgi)*time2sec(1,sMeta.dt,sMeta.dateCurr)*sCryo.hfnet/cLate;


%Update snow temperature based on conduction:
% dTemp = 10^(tsnPwr)*sCryo.hfnetSC./(cSens*sCryo.snw);
% dTemp(isnan(dTemp)) = 0;
% sCryo.tsn = sCryo.tsn + dTemp;



   
%Melt potential all goes into melting snow but is limited by amount of snow
%available
indMelt = find( sCryo.lhpme > 0);
if ~isempty(indMelt)
    sCryo.lhsnme(indMelt) = sCryo.lhpme(indMelt);
    indMaxSnow = find(sCryo.lhsnme(indMelt) > sCryo.snw(indMelt));
    sCryo.lhsnme(indMelt(indMaxSnow)) = sCryo.snw(indMelt(indMaxSnow));
    sCryo.lhpme(indMelt) = sCryo.lhpme(indMelt) - sCryo.lhsnme(indMelt);
    
    %Add melted snow to 'release' field and remove from 'solid':
    sCryo.lwsnl(indMelt) = sCryo.lwsnl(indMelt) + sCryo.lhsnme(indMelt);
    sCryo.snw(indMelt) = sCryo.snw(indMelt) - sCryo.lhsnme(indMelt);
end


%Modify snowpack temperature:
if isfield(sCryo,'hfsnc')
    %Negative means that heat is conducted from inner area to surface
    sCryo.tsn = sCryo.tsn - sCryo.hfsnc./(cSens*sCryo.snw);
else
   error('mass_enbal:noConduct','The enbal conservation of energy and mass formulation requires a heat conduction term.'); 
end


%Refreeze liquid in snow if internal temperature < 0
%Don't freeze rain, because that water is already in snow liquid content
%field
indLiqFrz = find(sCryo.tsn < 0 & sCryo.lwsnl > 0 & sCryo.snw > 0); %Indices where there is positive cold content and there is liquid water in snow
if ~isempty(indLiqFrz)
    %Indices where there's more liquid water than energy to freeze
    indFrzMax = find(sCryo.lwsnl(indLiqFrz) > -(cSens/cLate)*sCryo.snw(indLiqFrz).*sCryo.tsn(indLiqFrz) );
    if ~isempty(indFrzMax)
        frzMax = -(cSens/cLate)*sCryo.snw(indLiqFrz(indFrzMax)).*sCryo.tsn(indLiqFrz(indFrzMax));
        sCryo.snw(indLiqFrz(indFrzMax)) = sCryo.snw(indLiqFrz(indFrzMax)) + frzMax;
        sCryo.lwsnl(indLiqFrz(indFrzMax)) = sCryo.lwsnl(indLiqFrz(indFrzMax)) - frzMax;
        sCryo.tsn(indLiqFrz(indFrzMax)) = 0;
    end
    
    %Indices where there's sufficient energy to freeze all liquid:
    if isempty(indFrzMax)
        indFrzNMax = indLiqFrz;
    else
        indFrzNMax = setdiff(indLiqFrz, indLiqFrz(indFrzMax));
    end
    if ~isempty(indFrzNMax)
        sCryo.snw(indFrzNMax) = sCryo.snw(indFrzNMax) + sCryo.lwsnl(indFrzNMax);
        sCryo.tsn(indFrzNMax) = sCryo.tsn(indFrzNMax) + (cLate/cSens)*sCryo.lwsnl(indFrzNMax)./sCryo.snw(indFrzNMax);
        sCryo.lwsnl(indFrzNMax) = 0;
    end
end





%%MELT ICE:
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
    maxRat = 4;
    ratioMelt(ratioMelt > maxRat) = maxRat;
    ratioMelt(ratioMelt < 0) = 0; %Can't have negative ice energy
    sCryo.lhicme(indIceMlt) = 10^(dii)*sCryo.lhpme(indIceMlt).*ratioMelt;
        sCryo.lhpme(indIceMlt) = 0;
        
    %Scale melt for debris covered glacier cells:
    if isfield(sCryo,'icdbr')
        indDebrisMlt = intersect(indIceMlt, find(~isnan(sCryo.icdbr))); 
    	sCryo.lhicme(indDebrisMlt) = debScl*sCryo.lhicme(indDebrisMlt);
    end
    sCryo.lhicme( isnan(sCryo.lhicme)) = 0;
    sCryo.iclr(indIceMlt) = sCryo.icx(indIceMlt).*sCryo.lhicme(indIceMlt);
    sCryo.icdwe(indIceMlt) = -sCryo.lhicme(indIceMlt);
end


%ADJUST SNOW INTERNAL TEMP AND OTHER PARAMETERS:
%Set netagive solid snow values to 0:
sCryo.snw(sCryo.snw < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.lwsnl(sCryo.lwsnl < 0 ) = 0;

sCryo.tsn(sCryo.tsn > 0) = 0;
sCryo.tsn(sCryo.tsn < minTemp) = minTemp; 
sCryo.tsn(sCryo.lwsnl > 0.005*sCryo.snw) = 0;
sCryo.tsn(sCryo.snw == 0) = nan; 


%Set ice Temp (locations where no snow):
indBIce = find(sCryo.icx ~= 0 & sCryo.snw == 0);
    indBIce3 = ind2_3d(size(sAtm.tas),indBIce, sAtm.indCurr);
%     sCryo.tsn(indBIce) = 0;
%     sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tas(indBIce3)));
%Temprature of bare ice is min(0, air temp):
if isfield(sAtm,'tasmin')
    sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tasmin(indBIce3)));
else
    sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tas(indBIce3)));
end