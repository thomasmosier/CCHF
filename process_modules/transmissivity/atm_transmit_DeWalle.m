function varargout = atm_transmit_DeWalle(varargin)

%Eq. 6.7 on pg. 151 of DeWalle and Rango 2008
 
global sAtm

%VERSION WITHOUT FITTING PARAMETERS:
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1},{'atm_lin_scalar', 0, 30, 16.2, 'atm_transmit_DeWalle','cryo'});
    varargout{1} = cat(1,varargout{1},{'clear_trans', 0.5, 1, 0.65, 'atm_transmit_DeWalle','cryo'});
 
    return
else
    a = find_att(varargin{1}.coef,'atm_lin_scalar'); 
    b = find_att(varargin{1}.coef,'clear_trans'); 
end



%Calculate current temperature range and
if isfield(sAtm,'tasmin') && isfield(sAtm,'tasmax')
    sAtm.tasrng = squeeze(sAtm.tasmax(sAtm.indCurr,:,:)) - squeeze(sAtm.tasmin(sAtm.indCurr,:,:)); 
else
    error('atm_transmit_DeWalle:tasMaxMin','tasmax and tasmin are not present, but are needed for present function.');
end
    
    
sAtm.rstran =  b*(1-exp(-0.01*a*sAtm.tasrng.^2.4));

sAtm.rstran(sAtm.rstran > 1) = 1;
sAtm.rstran(sAtm.rstran < 0) = 0; %This line shouldn't be necessary



