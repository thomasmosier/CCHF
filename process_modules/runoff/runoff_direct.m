function varargout = runoff_direct(sHydro, varargin)

global sCryo sLand

if isempty(varargin(:))
	varargout{1} = cell(1,6);
    varargout{1}(1,:) = {'runoff_factor', 0, 1, 0.928, 'direct_runoff', 'routing'}; %percent runoff (unitless)
        
    return
else
    sMeta = varargin{1};
    rf = find_att(sMeta.coef,'runoff_factor');   
end


%Reset runoff field:
sLand.mrro = zeros(size(sLand.mrro));

indX = [];
if isfield(sCryo, 'ice')
   indX = find(sCryo.icx > 0);
end

%Account for PET
if isfield(sLand,'pet')
    %Set ET so that it at most equals available liquid:
    indPET = find(ismember(sLand.datepet(:,2:end),sMeta.dateCurr(2:end), 'rows') == 1);
    if isempty(indPET)
        error('groundwater_bucket:noPETcurr','No PET grid was calculated for the current time-step');
    end
    
    sLand.et = squeeze(sLand.pet(indPET,:,:));
    
    %Set ET to 0 when snowpack or glaciers exist:
    sLand.et(sCryo.snw > 0) = 0;
    if ~isempty(indX)
        sLand.et(indX) = sCryo.icx(indX).*sLand.et(indX);
    end
else
    error('direct_runoff:noPET','There should be a potential evapotranspiration field.');
end

%Calculate total input water as snowmelt plus rain (will account for ET and direct runoff later):
%**Rain added to sCryo.release in snow melt functions
sLand.mrro = rf*sHydro.area.*(sLand.rnrf + sCryo.snlr + sCryo.iclr - sLand.et);
    sLand.mrro(sLand.mrro < 0 ) = 0;
