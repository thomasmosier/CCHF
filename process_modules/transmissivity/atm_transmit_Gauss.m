function varargout = atm_transmit_Gauss(varargin)

global sAtm
%Model is 'Gaussian' parameterization of cloud factor from:
%Pellicciotti, F., Raschle, T., Huerlimann, T., Carenzo, M., & Burlando, P. 
%(2011). Transmission of solar radiation through clouds on melting 
%glaciers: a comparison of parameterizations and their impact on melt 
%modelling. Journal of Glaciology, 57(202), 367-381.
%optimal values in reference are:
    %a = 0.9561, b = 11.51, and c = 10.62

    
% if isempty(varargin(:))
% 	argout = cell(3,6);
%     argout(1,:) = {'atm_lin_scalar', 0.5, 1.5, 0.9561, 'atm_transmit_Gauss','cryo'}; 
%     argout(2,:) = {'atm_temp_corr', 0, 20, 11.51, 'atm_transmit_Gauss','cryo'}; 
%     argout(3,:) = {'atm_temp_scale', 5, 20, 10.62, 'atm_transmit_Gauss','cryo'}; 
%     return
% else
%     a = find_att(varargin{1}.coef,'atm_lin_scalar'); 
%     b = find_att(varargin{1}.coef,'atm_temp_corr'); 
%     c = find_att(varargin{1}.coef,'atm_temp_scale'); 
% end

% argout = a*exp(-((dTemp-b)/c).^2);

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'atm_temp_scalar', 3, 30, 11.5, 'atm_transmit_Gauss','cryo'}); 
    return
else
    a = find_att(varargin{1}.coef,'atm_temp_scalar'); 
end


varargout = 0.9561*exp(-((sAtm.tasRange-a)/10.62).^2);