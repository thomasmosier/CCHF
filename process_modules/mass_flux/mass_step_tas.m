function varargout = mass_step_tas(varargin)

    
global sCryo sAtm

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'tmp_melt_offset', 0,   4,    1, 'heat_melt_threshold','cryo'}); %Units of depth melt
        
    return
else
    tmpOffset = find_att(varargin{1}.coef,'tmp_melt_offset'); 
    sMeta = varargin{1};
end

%This conversion complements a similar one in 'ETI_Pellicciotti' and
%'simple_degree'. Makes units consistent with SETI
densW = find_att(sMeta.global,'density_water'); 
cLate = densW*find_att(sMeta.global,'latent_water');  %Density water * Latent heat fusion ice; %Units = J/m^3

%Initialize release field:
sCryo.lhsnme = zeros(size(sCryo.snw),'single');
% sCryo.lhicme  = zeros(size(sCryo.snw),'single');
if ~isfield(sCryo, 'lwsnl')
    sCryo.lwsnl = zeros(size(sCryo.snw),'single');
end


%%ORIGINAL METHOD:
% %Freeze rain if temperature below threshold:
% indRainFrz = find(sAtm.rain > 0 & squeeze(sAtm.tas(sAtm.indCurr,:,:)) <= tMelt);
% if ~isempty(indRainFrz)
%     sCryo.snw(indRainFrz) = sCryo.snw(indRainFrz) + sAtm.rain(indRainFrz);
%     sAtm.rain(indRainFrz) = 0;
% end
%
% %Set melt potential to zero at cells where temperature less than threshold:
% sCryo.lhpme(squeeze(sAtm.tas(sAtm.indCurr,:,:)) <= tMelt ) = 0;

%%NEW TEST METHOD:
%Calculate melt potential using simple degree indec formulation (units of m): 
sCryo.lhpme = time2sec(1,sMeta.dt,sMeta.dateCurr)*sCryo.hfnet/cLate;
% %Set melt potential to zero at cells where temperature less than threshold:
% sCryo.lhpme(squeeze(sAtm.tas(sAtm.indCurr,:,:)) <= tMelt ) = 0;
sCryo.lhpme(squeeze(sAtm.tas(sAtm.indCurr,:,:)) <= tmpOffset ) = 0;
sCryo.lhpme(sCryo.lhpme < 0 ) = 0;


% %If liquid water exists in snowpack and sCryo.lhpme == 0, add to solid:
% indReFrz = find(sCryo.lwsnl > 0 & sCryo.lhpme == 0);
% if ~isempty(indReFrz)
%     sCryo.snw(indReFrz) = sCryo.snw(indReFrz) + sCryo.lwsnl(indReFrz);
%     sCryo.lwsnl(indReFrz) = 0;
% end



%Melt potential all goes into melting snow but is limited by amount of snow
%available
indMelt = find( sCryo.lhpme > 0);
if ~isempty(indMelt)
    sCryo.lhsnme(indMelt) = sCryo.lhpme(indMelt);
    indMaxSnow = find(sCryo.lhsnme(indMelt) > sCryo.snw(indMelt));
    sCryo.lhsnme(indMelt(indMaxSnow)) = sCryo.snw(indMelt(indMaxSnow));
    
    sCryo.lhpme(indMelt) = sCryo.lhpme(indMelt) - sCryo.lhsnme(indMelt);
    
    %Add melted snow to 'liquid' field and remove from 'solid':
    sCryo.lwsnl(indMelt) = sCryo.lwsnl(indMelt) + sCryo.lhsnme(indMelt);
    sCryo.snw(indMelt) = sCryo.snw(indMelt) - sCryo.lhsnme(indMelt);
end


%Set netagive solid snow values to 0:
sCryo.snw(sCryo.snw < 0 ) = 0;
%Set negative snow liquid values to 0:
sCryo.lwsnl(sCryo.lwsnl < 0 ) = 0;


% %%MELT ICE:
% %if snowpack has ice field, allow additional melt at those locations:
% %Only finds indices where ice exists and where melt potential greater than
% %snowpack:
% if isfield(sCryo,'icx')
%     indIceMlt = find( sCryo.icx ~= 0 & sCryo.lhpme > 0);
% else
%     indIceMlt = [];
% end
% 
% %Calculate changes using residual energy at locations with ice:
% if ~isempty(indIceMlt)
%     %Calculate ice melt, which involves different degree-index factor
%     sCryo.lhicme(indIceMlt) = 10^(dii)*sCryo.lhpme(indIceMlt).*sCryo.hfneti(indIceMlt)./sCryo.hfnet(indIceMlt);
%         sCryo.lhicme( isnan(sCryo.lhicme)) = 0;
%     sCryo.iclr(indIceMlt) = sCryo.icx(indIceMlt).*sCryo.lhicme(indIceMlt);
%     sCryo.icdwe(indIceMlt) = -sCryo.lhicme(indIceMlt);
% end