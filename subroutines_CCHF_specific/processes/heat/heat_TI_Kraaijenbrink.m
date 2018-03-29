function varargout = heat_TI_Kraaijenbrink(varargin)
%Reference:
%Kraaijenbrink, P. D. A., Bierkens, M. F. P., Lutz, A. F., &
%Immerzeel, W. W. (2017). Impact of a global temperature rise of 1.5
%degrees Celsius on Asia?s glaciers. Nature, 549(7671), 257?260.
%https://doi.org/10.1038/nature23878


%Notes:
%Different degree index factors for snow and ice
%Ice ponds melt at 10 times the rate as thick debris
%Uses empirical curve to estimate melt under debris as function of debris
%thickness



global sCryo sAtm


if isempty(varargin(:))
    varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'watt_per_deg_snow', 10, 30, 20, 'heat_TI_Kraaijenbrink','cryo'});
    varargout{1} = cat(1,varargout{1}, {'watt_per_deg_ice',  20, 50, 40, 'heat_TI_Kraaijenbrink','cryo'}); %Units of depth melt

    %The degree index factor here is approx 3.9 times degree-day factors
    %with units of mm / C / day. Shea et al., 2015 find snow DD factor of
    %about 5 mm / C / day, which is approximately 20, and clean ice DD
    %factor of 9.7, which is 40 (range in their paper is 1 SD)
    return
else
    wattperdegS = find_att(varargin{1}.coef,'watt_per_deg_snow');
    wattperdegI = find_att(varargin{1}.coef, 'watt_per_deg_ice');
end


% fConvert = 3.34*10^8/3600; %Latent heat of fusion divided by seconds in hour. Comes from (L_f*density_water)/(1000mm/m * 3600 s/hr)
% %Units of J-hr-m^{-2}-s^{-1}

szGrid = size(sAtm.rain);

indPos = squeeze(sAtm.tas(sAtm.indtas,:,:)) > 0;

%Calculate equivalent of 'heatflux' using simple degree index model
sCryo.hfnet = zeros(szGrid, 'single');
sCryo.hfnet(indPos) = wattperdegS; %units to Watts per m^2

%Heat flux for ice
sCryo.hfneti = zeros(szGrid, 'single');
sCryo.hfneti(indPos) = wattperdegI;

% %For testing the empirical melt function:
% test = (0:0.2:50);
% meltTest = debris_melt_empirical(test, 'cm');
% figure; scatter(test, meltTest);
%
% test = (0:0.005:0.5);
% meltTest = debris_melt_empirical(test, 'm');
% figure; scatter(test, meltTest);

%Modify heat flux based on debris covered area
if isfield(sCryo, 'icdbr') %If debris cover information available
    if ~isfield(sCryo, 'icdbrmelt')
        sCryo.icdbrmelt = debris_melt_empirical(sCryo.icdbr, 'm');
    end

    sCryo.hfneti = sCryo.icdbrmelt.*sCryo.hfneti;
end

if isfield(sCryo, 'icpndx')
    %Melt is enhanced by factor of 10 at locations that are entirely ponds.
    %Use fractional relation to account for partially ponded grid cells
    %(when fraction is 0, there is no impact; when fraction is 1, the
    %factor is 10)
    sCryo.hfneti = (1+9*sCryo.icpndx).*sCryo.hfneti;
end
