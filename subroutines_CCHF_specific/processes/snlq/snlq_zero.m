function varargout = snlq_zero(varargin)

%holdCap = decimal between 0 and 1


global sCryo

%WITHOUT HOLDING CAPACITY
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    return
end



%Initialize melt release array:
sCryo.snlr = zeros(size(sCryo.snw),'single');

%%Release liquid in excess of snow holding capacity:
%Find holding capacity (percentage of solid snow):
sCryo.snlh = zeros(size(sCryo.snw));
%Find indices where it's exceeded:
indRelease = find(sCryo.lwsnl > 0);
if ~isempty(indRelease)
    %Amount of release equals exceedance of liquid water holding capacity:
    sCryo.snlr(indRelease) = sCryo.lwsnl(indRelease);
    sCryo.sndwe(indRelease) = sCryo.sndwe(indRelease) - sCryo.snlr(indRelease);
    %Remove drained water from snowpack liquid water content:
    sCryo.lwsnl = zeros(size(sCryo.snw));
end
