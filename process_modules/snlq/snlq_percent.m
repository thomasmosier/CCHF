function varargout = snlq_percent(varargin)

%holdCap = decimal between 0 and 1


global sCryo

%WITHOUT HOLDING CAPACITY
if isempty(varargin(:))
	varargout{1} = cell(0,6);
    varargout{1} = cat(1,varargout{1}, {'sn_hold', 0,   .1,    0.077, 'snlq_percent','cryo'}); %Units of depth melt
    return
else
    holdCap = find_att(varargin{1}.coef,'sn_hold','no_warning'); 
%     if isempty(holdCap)
%        holdCap = 0; 
%     end
%     a = 0.03;  
end


%Initialize melt release array:
sCryo.snlr = zeros(size(sCryo.snw),'single');

%%Release liquid in excess of snow holding capacity:
%Find holding capacity (percentage of solid snow):
sCryo.snlh = holdCap*sCryo.snw;
%Find indices where it's exceeded:
indRelease = find(sCryo.lwsnl > sCryo.snlh);
if ~isempty(indRelease)
    %Amount of release equals exceedance of liquid water holding capacity:
    sCryo.snlr(indRelease) = sCryo.lwsnl(indRelease) - sCryo.snlh(indRelease);
    sCryo.sndwe(indRelease) = sCryo.sndwe(indRelease) - sCryo.snlr(indRelease);
    %Remove drained water from snowpack liquid water content:
    sCryo.lwsnl(indRelease) = sCryo.snlh(indRelease);
end
