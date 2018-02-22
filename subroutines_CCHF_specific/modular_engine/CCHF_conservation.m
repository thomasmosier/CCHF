function CCHF_conservation(sMeta)
%Function to check that mass and energy are conserved during each timestep
global sCryo sAtm sLand

keyboard


%%CONSERVATION OF MASS:
%Precip = rain + snow
if ~all(isequal(squeeze(sAtm.pr(sAtm.indpr,:,:)), sAtm.rain + sAtm.prsn))
    strDate = date_2_string(sMeta.dateCurr);
    error('cchfConservation:precipitationNotConserved', ['Precipitation does not equal rainfall plus snowfall for ' strDate]);
end
      
%Snowpack changes
if ~all(isequal(sCryo.sndwe, sAtm.prsn - sCryo.snsb - sCryo.snlr + sCryo.sndlw))
    
end

%Soil moisture / Runoff
if isfield(sLand, 'sm')
    if ~all(isequal(sLand.smin - sLand.smout, sLand.rnrf + sCryo.snlr + sCryo.iclr - sLand.et - sAtm.mrro))
        strDate = date_2_string(sMeta.dateCurr);
        error('cchfConservation:soilMoistureNotConserved', ['Change in soil moisture does not match inflow and outflow sources for ' strDate]);
    end
end

%Annual discharge?


% %Find heat fields:
% fldHeat = 'fldsheat';
% if ~isfield(sCryo, fldHeat)
%     fldsCryo = fieldnames(sCryo);
%     
%     for ii = numel(fldsCryo) : -1 : 1
%         if ~strcmpi(fldsCryo{ii}(1:2), 'hf') || regexpbl(fldsCryo{ii}, 'net') || strcmpi(fldsCryo{ii}(end), 'i')
%             fldsCryo(ii) = [];
%         end
%     end
% 
%     sCryo.(fldHeat) = fldsCryo;
% end
% 
% %What is 'lhpme'?
% fldLatentSnow = 'lhsnme';
% fldLatentIce = 'lhicme';
% fldSensibleSnow = 1;
% %Find latent heat energy fields:
% fldLatent = 'fldslatent';
% if ~isfield(sCryo, fldLatent)
%     fldsCryo = fieldnames(sCryo);
%     
%     for ii = numel(fldsCryo) : -1 : 1
%         if ~strcmpi(fldsCryo{ii}(1:2), 'lh') || strcmpi(fldsCryo{ii}(end), 'i')
%             fldsCryo(ii) = [];
%         end
%     end
% 
%     sCryo.(fldLatent) = fldsCryo;
% end

%Find sensible heat energy fields





%Precip = rain + snow

%Snowpack changes

%runoff

%Discharge
