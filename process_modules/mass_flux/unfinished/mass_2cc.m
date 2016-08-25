function varargout = mass_2cc(varargin)

    
global sCryo sAtm

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'cc_pwr', -3,   2,    -1, 'melt_cc','cryo'}); %Unitless; impacts accumulation of cold-content (negative energy)
    varargout{1} = cat(1,varargout{1}, {'ice_melt_pwr' , -2, 2, 0.11, 'melt_cc','cryo'});  %Units of mm/hr/Deg C
    
    if isfield(sCryo,'tsn') || isfield(sCryo,'tsnmin')
        varargout{1} = cat(1,varargout{1}, {'tsn_pwr', -3,   2,    -1, 'melt_cc','cryo'});
    end
    
    if isfield(sCryo,'debris')
        varargout{1} = cat(1,varargout{1}, {'debris_scalar' , 0, 2, 1, 'melt_cc','cryo'});  %Unitless
    end
    
    return
else
    ccPwr = find_att(varargin{1}.coef,'cc_pwr'); 
    dii = find_att(varargin{1}.coef,'ice_melt_pwr');  %Units of mm/hr/Deg C
    
    if isfield(sCryo,'tsn')
        tsnPwr = find_att(varargin{1}.coef,'tsn_pwr');
    end
    
    if isfield(sCryo,'debris')
        debScl = find_att(varargin{1}.coef,'debris_scalar');  %Unitless
    end
    
    sMeta = varargin{1};
end

%This conversion complements a similar one in 'ETI_Pellicciotti' and
%'simple_degree'. Makes units consistent with SETI
cLate = 3.34*10^8; %Density water * Latent heat fusion ice; %Units = J/m^3

%Initialize release field:
sCryo.latentSnow = zeros(size(sCryo.solid),'single');
sCryo.latentIce  = zeros(size(sCryo.solid),'single');
sCryo.dIce       = zeros(size(sCryo.solid),'single');
sCryo.ilr        = zeros(size(sCryo.solid),'single'); %Ice liquid release
if ~isfield(sCryo, 'snlq')
    sCryo.snlq = zeros(size(sCryo.solid),'single');
end
if ~isfield(sCryo, 'sncc')
    sCryo.sncc = zeros(size(sCryo.solid),'single');
end
if ~isfield(sCryo, 'tic')
    sCryo.tic = zeros(size(sCryo.solid),'single');
end


% if isequal(sMeta.currTime, [2003,5]) || isequal(sMeta.currTime, [2003,6]) 
%    keyboard 
% end

%Calculate melt potential using simple degree indec formulation (units of m): 
sCryo.lhme = time2sec(1,sMeta.dt,sMeta.currTime)*sCryo.heat/cLate;


%Update cold-conent:
%Add thermal inertia to the cold-content parameter:
indNegEn = find(sCryo.lhme < 0);
if ~isempty(indNegEn)
   sCryo.sncc(indNegEn) = sCryo.sncc(indNegEn) - 10^(ccPwr)*sCryo.lhme(indNegEn);
end

%Remove cold-content and melt energy:
indPosEn = find(sCryo.lhme > 0 & sCryo.sncc > 0);
if ~isempty(indPosEn)
    tempDcc = sCryo.sncc(indPosEn) - sCryo.lhme(indPosEn);
    indMaxCC = find(tempDcc < 0);
    if ~isempty(indMaxCC)
        sCryo.lhme(indPosEn(indMaxCC)) = -tempDcc(indMaxCC);
        sCryo.sncc(indPosEn(indMaxCC)) = 0;
    end
    indNMaxCC = setdiff((1:numel(indPosEn)),indMaxCC);
    if ~isempty(indNMaxCC)
        sCryo.sncc(indPosEn(indNMaxCC)) = tempDcc(indNMaxCC);
        sCryo.lhme(indPosEn(indNMaxCC)) = 0;
    end
end


%Reset cold content if no snow:
sCryo.sncc(sCryo.sncc ~= 0 & sCryo.solid == 0) = 0;



%
%Refreeze liquid in snow; if snow contains liquid and there is
%rain, release rain, otherwise hold in snowpack

indLiqFrz = find(sCryo.sncc > 0 & sCryo.snlq > 0 & sCryo.solid > 0); %Indices where there is positive cold content and there is liquid water in snow
if ~isempty(indLiqFrz)
    %Indices where there's more liquid water than energy to freeze
    indFrzMax = find(sCryo.snlq(indLiqFrz) > sCryo.sncc(indLiqFrz) );
    if ~isempty(indFrzMax)
        sCryo.solid(indLiqFrz(indFrzMax)) = sCryo.solid(indLiqFrz(indFrzMax)) + sCryo.sncc(indLiqFrz(indFrzMax));
        sCryo.snlq(indLiqFrz(indFrzMax)) = sCryo.snlq(indLiqFrz(indFrzMax)) - sCryo.sncc(indLiqFrz(indFrzMax));
        sCryo.sncc(indLiqFrz(indFrzMax)) = 0;
    end
    
    %Indices where there's sufficient energy to freeze all liquid:
    if isempty(indFrzMax)
        indFrzNMax = indLiqFrz;
    else
        indFrzNMax = setdiff(indLiqFrz, indLiqFrz(indFrzMax));
    end
    if ~isempty(indFrzNMax)
        sCryo.solid(indFrzNMax) = sCryo.solid(indFrzNMax) + sCryo.snlq(indFrzNMax);
        sCryo.sncc(indFrzNMax) = sCryo.sncc(indFrzNMax) - sCryo.snlq(indFrzNMax);
        sCryo.snlq(indFrzNMax) = 0;
    end
end


%Melt potential all goes into melting snow but is limited by amount of snow
%available
indMelt = find( sCryo.lhme > 0 & sCryo.sncc <= 0);
if ~isempty(indMelt)
    sCryo.latentSnow(indMelt) = sCryo.lhme(indMelt);
    indMaxSnow = find(sCryo.latentSnow(indMelt) > sCryo.solid(indMelt));
    sCryo.latentSnow(indMelt(indMaxSnow)) = sCryo.solid(indMelt(indMaxSnow));
    sCryo.lhme(indMelt) = sCryo.lhme(indMelt) - sCryo.latentSnow(indMelt);
    
    %Add melted snow to 'release' field and remove from 'solid':
    sCryo.snlq(indMelt) = sCryo.snlq(indMelt) + sCryo.latentSnow(indMelt);
    sCryo.solid(indMelt) = sCryo.solid(indMelt) - sCryo.latentSnow(indMelt);
end

%%MELT ICE:
%if snowpack has ice field, allow additional melt at those locations:
%Only finds indices where ice exists and where melt potential greater than
%snowpack:
if isfield(sCryo,'ice')
    indIceMlt = find( sCryo.ice ~= 0 & ~isnan(sCryo.ice) & sCryo.lhme > 0);
else
    indIceMlt = [];
end

%Calculate changes using residual energy at locations with ice:
if ~isempty(indIceMlt)
    %Calculate ice melt, which involves different degree-index factor
    sCryo.latentIce(indIceMlt) = 10^(dii)*sCryo.lhme(indIceMlt).*sCryo.heatI(indIceMlt)./sCryo.heat(indIceMlt);
        sCryo.lhme(indIceMlt) = 0;
        
    %Scale melt for debris covered glacier cells:
    if isfield(sCryo,'debris')
        indDebrisMlt = intersect(indIceMlt, find(~isnan(sCryo.debris))); 
    	sCryo.latentIce(indDebrisMlt) = debScl*sCryo.latentIce(indDebrisMlt);
    end
    sCryo.latentIce( isnan(sCryo.latentIce)) = 0;
    sCryo.ilr(indIceMlt) = sCryo.latentIce(indIceMlt);
    sCryo.dIce(indIceMlt) = -sCryo.latentIce(indIceMlt);
end


%Set negative cold content values to 0:
sCryo.sncc(sCryo.sncc < 0) = 0;
%Set netagive solid snow values to 0:
sCryo.solid(sCryo.solid < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.snlq(sCryo.snlq < 0 ) = 0;


%ADJUST SNOW TEMPERATURE BASED ON COLD-CONTENT
if isfield(sCryo,'tsn')
    cSensW = 4.18*10^6; %= c_water*density_water: J / m^3 / Deg C
    cSensI = 2.027*10^6; %= c_ice*density_water: J / m^3 / Deg C (specific heat capacity of ice/snow ranges from about 2030-2090 J/kg/K)
    
    indSnow = find(sCryo.solid ~= 0);
    %From def. of cold-content, T = (cc*rho_w*L)/(rho_s*c*SWE)
    sCryo.tsn(indSnow) = -10^(tsnPwr)*(cLate*sCryo.sncc(indSnow)./(cSensI*sCryo.solid(indSnow) + cSensW*sCryo.snlq(indSnow)));
        %In the above, there shouldn't actually be any liquid water at
        %locations where cc is positive
        
    %Find indices of bare ice:
    indBIce = setdiff(find(sCryo.ice ~= 0), indSnow);
        indBIce3 = ind2_3d(size(sAtm.tas),indBIce, sAtm.indCurr);
    %Temprature of bare ice is min(0, air temp):
    if isfield(sAtm,'tasmin')
        sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tasmin(indBIce3)));
    else
        sCryo.tsn(indBIce) = min(0, squeeze(sAtm.tas(indBIce3)));
    end

%     %Temperature of debris covered ice is mean air temperature:
%     if isfield(sCryo,'debris')
%         indDebris = intersect(indBIce, find(~isnan(sCryo.debris)));
%             indDebris3 = ind2_3d(size(sAtm.tas),indDebris, sAtm.indCurr);
% 		if isfield(sAtm,'tasmin')	
% 			sCryo.tsn(indDebris) = squeeze(sAtm.tasmin(indDebris3));
% 		else
% 			sCryo.tsn(indDebris) = squeeze(sAtm.tas(indDebris3));
% 		end
%     end
    
    %%MAKE CORRECTIONS TO VARIOUS PROPERTIES:
    %Set temperature to 0 where there is liquid water:
    sCryo.tsn(sCryo.snlq > 0) = 0;
    %Set temp to zero if solid snow exists and temp positive:
    sCryo.tsn(sCryo.tsn > 0 & sCryo.solid > 0) = 0;
    %Set snow temperature to nan where there is no snow or ice:
    sCryo.tsn(sCryo.solid == 0 & sCryo.ice <= 0) = nan;

    % %Set negative cold content to 0:
    % sCryo.sncc(sCryo.sncc < 0) = 0;

    if any(any(sCryo.tsn < sCryo.tsnmin & sCryo.solid > 0.01))
        if ~regexpbl(sMeta.mode, 'calib')
            warning('snowpack_enbal:negSnowTmp',['Snowpack temperature is less than the minimum, which is ' ...
                num2str(sCryo.tsnmin) ' for at least one cell. It is being set to the minimum.']);
        end
        sCryo.tsn(  sCryo.tsn < sCryo.tsnmin) = sCryo.tsnmin; 
    end 
    
end

