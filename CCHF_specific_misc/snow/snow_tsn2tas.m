function snow_tsn_reset
%If no snow at cell, set snow temperature to air temperature 
%(may matter in energy balance):

global sCryo


if ~isfield(sCryo,'tsn')
    sCryo.tsn = zeros(size(sCryo.solid),'single');
end
% if ~isfield(sCryo,'dTmp')
%     sCryo.dTmp = zeros(size(sCryo.solid),'single');
% end

% ind2NoSnow = find(sCryo.solid == 0 & sCryo.liquid == 0);
% [snowRow, snowCol] = ind2sub(size(sCryo.solid), ind2NoSnow);
% szTmp = size(sAtm.tas);
% ind3NoSnow = sAtm.indCurr + (snowRow-1)*szTmp(1) + (snowCol-1)*szTmp(2)*szTmp(1);
% sCryo.tsn(ind2NoSnow) = sAtm.tas(ind3NoSnow);

% sCryo.dTmp(ind2NoSnow) = 0;

sCryo.tsn(sCryo.solid == 0 & sCryo.liquid == 0) = nan;

%If snow temp = nan and there is snow, set to 0
sCryo.tsn(isnan(sCryo.tsn) & (sCryo.solid ~= 0 | sCryo.liquid ~= 0)) = 0;