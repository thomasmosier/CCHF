function varargout = atm_transmit_Coops(sHydro,varargin)

%Coops, N. C., Waring, R. H., & Moncrieff, J. B. (2000). Estimating mean 
%monthly incident solar radiation on horizontal and inclined slopes from 
%mean monthly temperatures extremes. International Journal of 
%Biometeorology, 44(4), 204–211.

 
global sAtm

%VERSION WITHOUT FITTING PARAMETERS:
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'clear_trans_pwr',   -2, 2, 0, 'atm_tansmit_Coops','cryo'});
    varargout{1} = cat(1,varargout{1}, {'atm_scl',   0, 2, 1, 'atm_tansmit_Coops','cryo'});
    return
else
    demFact = find_att(varargin{1}.coef,'clear_trans_pwr'); 
    scl = find_att(varargin{1}.coef,'atm_scl'); 
end


%Calculate current temperature range and
if isfield(sAtm,'tasmin') && isfield(sAtm,'tasmax')
    sAtm.tasrng = squeeze(sAtm.tasmax(sAtm.indCurr,:,:)) - squeeze(sAtm.tasmin(sAtm.indCurr,:,:)); 
else
    error('atm_transmit_Coops:tasMaxMin','tasmax and tasmin are not present, but are needed for present function.');
end

tClear = (0.65+0.001*demFact*sHydro.dem);
tClear(tClear > 1) = 1;

sAtm.rstran = tClear.*(1-exp(-scl*(0.031+0.2*exp(-.185*sAtm.tasrng)).*sAtm.tasrng.^1.5));



