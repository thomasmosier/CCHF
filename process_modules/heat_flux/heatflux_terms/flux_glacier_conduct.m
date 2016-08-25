function varargout = flux_glacier_conduct(varargin)

global sCryo

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1, varargout{1}, {'glacier_temp', -25, 0, -15, 'heat_conduction_const','cryo'});  
    varargout{1} = cat(1, varargout{1}, {'glacier_iso_depth', 0, 30, 10, 'heat_conduction_const','cryo'});  
    varargout{1} = cat(1, varargout{1}, ['tsn',cell(1,5)]);
    return
else
    tGlac = find_att(varargin{1}.coef,'glacier_temp');
    z = find_att(varargin{1}.coef,'glacier_iso_depth');
end

k = find_att(varargin{1}.global,'thermal_conduct_ice');

sCryo.hficc = zeros(size(sCryo.snw));


indIce = find(sCryo.icx > 0); 
sCryo.hficc(indIce) = -k *(sCryo.tsn(indIce) - tGlac)/z;
