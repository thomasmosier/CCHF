%MODIS move script

foldFrom = '/Users/thomas/Desktop/research/geodata/MODIS/MODIS-8day-500m/New Folder With Items 2';
fileMvMODIS = dir(fullfile(foldFrom, '*.hdf'));
fileMvMODIS = extractfield(fileMvMODIS, 'name');

foldTo = '/Users/thomas/Desktop/research/geodata/MODIS/MODIS-8day-500m';

disp(['There are ' num2str(numel(fileMvMODIS)) ' files to move.']);

for ii = 1 : numel(fileMvMODIS)
    if mod(ii,1000) == 0
       disp(['On iteration ' num2str(ii)]);
    end
    movefile(fullfile(foldFrom, fileMvMODIS{ii}), foldTo);
end