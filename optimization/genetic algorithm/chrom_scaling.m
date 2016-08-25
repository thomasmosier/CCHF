function paramLeng = chrom_scaling(prmSpecs, varargin)

if ~isempty(varargin)
   if numel(varargin{1}) == 1
       nPts = varargin{1}*ones(numel(prmSpecs(:,1)),1);
   elseif isequal(numel(varargin{1}), numel(prmSpecs(:,1)))
       nPts = varargin{1};
   else
       error('chrom_scaling:nDelta','The delta scaling parameter has an unexpected number of entries.')
   end
else
    nPts = 100*ones(numel(prmSpecs(:,1)),1);
end

% if iscell(prmSpecs)
%     prmSpecs = cell2mat(prmSpecs);
% end

paramLeng = ceil(log(nPts)/log(2));
nPts = 2.^paramLeng(:,1) - 1;
    paramLeng( nPts == 1 ) = 0;