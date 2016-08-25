function varargout = flux_conduct_c(varargin)

global sCryo

if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1, varargout{1}, {'ground_watts', -50, 0, -20, 'heat_conduction_const','cryo'});  
    return
else
    c = find_att(varargin{1}.coef,'ground_watts');
end

if ~isfield(sCryo, 'hfgc')
    sCryo.hfgc = c*ones(size(sCryo.solid));
end