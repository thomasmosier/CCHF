function varargout = scx_max(varargin)
%Assumes snow is uniformly distributed over grid cell


global sCryo

if isempty(varargin(:))
    varargout{1} = cell(0,6);
    
    return
end

thresh = 0.05; %SWE threshold is 5 cm of water equivalent for pixel to be classified as snow covered

%reset each time step
sCryo.scx = zeros(size(sCryo.snw),'single');

sCryo.scx(sCryo.snw >= thresh) = 1;