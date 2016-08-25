function sPath = clm_path_ui(sPath, sMeta)

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
    uiwait(msgbox(sprintf(['Select the folder containing ' sMeta.varLdDisp{ii} ...
        ' time-series data for ' sMeta.region '.\n']), ...
        '(Click OK to Proceed)','modal'));
    sPath.(sMeta.varLd{ii}) = uigetdir(startPath,['Select the folder containing ' sMeta.varLdDisp{ii} ...
        ' time-series data for ' sMeta.region]);
    disp([char(39) sPath.(sMeta.varLd{ii}) char(39) ' has been chosen as the ' sMeta.varLdDisp{ii} ' time-series data.']);
    
    %Update search path:
    [startPath, ~, ~] = fileparts(sPath.(sMeta.varLd{ii}));
%     indStart = regexpi(startPath,filesep);
%     startPath = startPath(1:indStart(end)-1);
end
