function argout = flux_long_Hock(varargin)
%R_long from atm is parameterized using Deacon (1970)

global sCryo sAtm

%VERSION WITH ONLY ONE FITTING PARAMETER
if isempty(varargin(:))
	argout = cell(0,6);
    argout = cat(1,argout, {'lw_pwr_tsn', -2, 2, 0, 'flux_long_Hock', 'cryo'});
    argout = cat(1,argout, ['tsn',cell(1,5)]);
    return
else
    scaleOut = find_att(varargin{1}.coef,'lw_pwr'); 
end

eSky = 0.7; %Emmissivity of sky
% stefan = 5.67*10^(-8);

%Primary reference is:
%Hock, R. (2005). Glacier melt: a review of processes and their modelling. 
%Progress in Physical Geography, 29(3), 362–391.

argout = (10^scaleOut)*5.67*10^(-8)*(eSky*(squeeze(sAtm.tas(sAtm.indCurr,:,:)) + 273.15).^4 ...
    - (sCryo.tsn + 273.15).^4);