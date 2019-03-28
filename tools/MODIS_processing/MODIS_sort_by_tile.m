%MODIS move script

foldFrom = '/Volumes/geospatial_research/geodata/MODIS/MOD10A2_v6_High_Mountain_Asia';

%Find all files:
fileMODIS = dir(fullfile(foldFrom, '*.hdf'));
fileMODIS = extractfield(fileMODIS, 'name');

tileUniq = cell(numel(fileMODIS), 1);
%Find all tiles:
disp(['Assessing tile names present in ' num2str(numel(fileMODIS)) ' files.']);
for ii = 1 : numel(fileMODIS)
    if mod(ii,1000) == 0
       disp(['On file ' num2str(ii)]);
    end
    
    indPer = regexpi(fileMODIS{ii}, '\.');
    tileUniq{ii} = fileMODIS{ii}(indPer(2)+1:indPer(3)-1);
end
tileUniq = unique(tileUniq);
disp([num2str(numel(tileUniq(:))) ' unique tiles have been identified.']);

%%CREATE UNIQUE FOLDERS FOR EACH TILE:
foldTile = cell(numel(tileUniq(:)), 1);
for ii = 1 : numel(tileUniq(:))
    foldTile{ii} = fullfile(foldFrom, tileUniq{ii});

    if ~exist(foldTile{ii}, 'dir')
        mkdir(foldTile{ii});
    end
end




%%MOVE FILES:
disp(['Sorting ' num2str(numel(fileMODIS)) ' files.']);
for ii = 1 : numel(fileMODIS)
    if mod(ii,1000) == 0
       disp(['On file ' num2str(ii)]);
    end
    
    %Identify which folder to move to:
    indPer = regexpi(fileMODIS{ii}, '\.');
    indTile = strcmpi(tileUniq, fileMODIS{ii}(indPer(2)+1:indPer(3)-1));
    
    %Move file:
    movefile(fullfile(foldFrom, fileMODIS{ii}), foldTile{indTile});
end