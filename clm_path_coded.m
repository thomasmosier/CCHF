function sPath = clm_path_coded(sPath, sMeta)

    
if isfield(sPath, 'dem') %Use dem path (if available)
    [startPath, ~, ~] = fileparts(sPath.dem);
else
 	nmsTemp = fieldnames(sPath);
    if ~isempty(nmsTemp) %Pick first available field
        [startPath, ~, ~] = fileparts(sPath.(nmsTemp{1}));
    else
        startPath = pwd;
    end
end


for ii = 1 : numel(sMeta.varLd)
    switch sMeta.varLd{ii}
        case 'pr'
            sPath.(sMeta.varLd{ii}) = fullfile(startPath, 'CFSR_daily', 'pre');
        case 'tas'
            sPath.(sMeta.varLd{ii}) = fullfile(startPath, 'CFSR_daily', 'tasmean');
        case 'tasmax'
            sPath.(sMeta.varLd{ii}) = fullfile(startPath, 'CFSR_daily', 'tasmax');
        case 'tasmin'
            sPath.(sMeta.varLd{ii}) = fullfile(startPath, 'CFSR_daily', 'tasmin');
        otherwise
            error('clim_path_coded:unknownVar',['The variable ' sMeta.varLd{ii} ' has not been coded for.']);
    end
end
