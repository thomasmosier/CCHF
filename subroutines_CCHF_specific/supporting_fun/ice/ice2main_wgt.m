function [wgtMnF, varargout] = ice2main_wgt(igrdLon, igrdLat, mainLon, mainLat, igrd2main, indIgrdF, indIgrdT, varargin)

%Optional varargin can be used to provide fractional weighting along a
%given path
if isempty(varargin(:))
    igrdWgt = [];
else
    igrdWgt = varargin{1};
    
    nIgrd = numel(igrdLat)*numel(igrdLon);
    if isequal(size(igrdWgt), [numel(igrdLat), numel(igrdLon)])
        igrdWgt = spdiags(igrdWgt(:),0,nIgrd,nIgrd);
    elseif ~isequal(size(igrdWgt), [nIgrd, nIgrd])
        error('ice2mainWgt:wrongWgtSize','The optional input weighting argument has the incorrect size.')
    end
    
    mxWgt = max(igrdWgt(:));
    if mxWgt > 1
       error('ice2mainWgt:igrdWgtMax', ['The maximum igrd weight is ' ...
           num2str(mxWgt) '. The maximum value should be 1.']); 
    end
end

%translation only matters when igrd indices cross main hydrologic grid 
indCross = find(igrd2main(indIgrdF) ~= igrd2main(indIgrdT));
indIgrdF = indIgrdF(indCross);
indIgrdT = indIgrdT(indCross);

%Find unique pairs of main grid indices (multiple ice grid cells may flow to same set of main grid cells):
%NOTE: the same main grid cell may flow to multiple main grid cells
%due to differences in flow direction at fine grid scale.
[temp, ~, ~] = unique([full(igrd2main(indIgrdF)), full(igrd2main(indIgrdT))],'rows');
indMnF = temp(:,1);
indMnT = temp(:,2);


%%Calculate area-weighting factor for each main grid cell where avalanching crosses border. 
%These factors will be used during each time step to redistribute
%snow

%Need areas for the cells that snow is avalanching from
indMnArea  = unique(indMnF);
indIceArea = unique(indIgrdF);

%Initialize area arrays (same size as whole domain)
areaMain = sparse(numel(mainLat), numel(mainLon));
areaIce  = sparse(numel(igrdLat), numel(igrdLon));

szIgrd = size(areaIce);

%Calculate areas
%It could be the case this this function is not correct
areaMain(indMnArea) = area_grid_pt( indMnArea, mainLon, mainLat);
areaIce(indIceArea) = area_grid_pt(indIceArea, igrdLon, igrdLat);

%%Calculate fractional contributions for cells that snow is avalanching from:
%Note: fractional contributions come from sum of areas in ice grid
%cells that avalance
wgtMnF = zeros(numel(indMnF), 1, 'single');

indIgrd2MainAvF = igrd2main(indIgrdF);
indIgrd2MainAvT = igrd2main(indIgrdT);

%Loop over all main grid cells that are being looped over
for ii = 1 : numel(indMnF)
    %Find ice grid indices contributing to current grid cell and
    %calculate fraction
    if ~isempty(igrdWgt)
        indCurr = indIgrd2MainAvF == indMnF(ii) & indIgrd2MainAvT == indMnT(ii);
        
        wgtMnF(ii) = sum(full(igrdWgt(sub2ind(szIgrd,indIgrdF(indCurr),indIgrdT(indCurr)))*areaIce(indIgrdF(indCurr))))...
            /areaMain(indMnF(ii));
    else
        wgtMnF(ii) = sum(full(areaIce(indIgrdF(indIgrd2MainAvF == indMnF(ii) & indIgrd2MainAvT == indMnT(ii)))))...
            /areaMain(indMnF(ii));
    end
end

if nargout == 3
    varargout{1} = indMnF;
    varargout{2} = indMnT;
end