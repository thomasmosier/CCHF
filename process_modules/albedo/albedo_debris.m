function varargout = albedo_debris(varargin)

global sCryo


if isempty(varargin(:))
    varargout{1} = cell(0,6);
%     varargout{1} = cat(1,varargout{1}, {'albedo_fresh', 0.5, 1, 0.71, 'albedo_Pellicciotti','cryo'});

    return
else
    if isfield(sCryo,'icdbr')
        aDebris = find_att(sMeta.global,'albedo_debris');

        sCryo.ialb(~isnan(sCryo.icdbr)) = aDebris;
    end
end