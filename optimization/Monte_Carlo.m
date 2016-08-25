function coef = Monte_Carlo(prmBnds, varargin)

if ~isempty(varargin) && isnumeric(varargin{1})
   n = varargin{1};
else
   n = 1;
end

coef = (repmat(prmBnds(:,2) - prmBnds(:,1),[1,n]).*rand([numel(prmBnds(:,1)),n]) + repmat(prmBnds(:,1),[1,n]))';