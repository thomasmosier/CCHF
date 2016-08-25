%MODIS move script

foldFrom = '/Users/thomas/Desktop/research/geodata/MODIS/MODIS-8day-500m';

foldMove = 'repeat';

foldOut = fullfile(foldFrom, foldMove);

if ~exist(foldOut, 'dir')
    mkdir(foldOut);
end

fileMODIS = dir(fullfile(foldFrom, '*.hdf'));
fileMODIS = extractfield(fileMODIS, 'name');



disp(['There are ' num2str(numel(fileMODIS)) ' files to assess.']);

for ii = 1 : numel(fileMODIS)
    if mod(ii,1000) == 0
       disp(['On file ' num2str(ii)]);
    end
    
    if regexpbl(fileMODIS{ii}, '(')
        movefile(fullfile(foldFrom, fileMODIS{ii}), foldOut);
    end
end