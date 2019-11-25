%
%test script for new DHM format streamgage data at
%read_gagedata level.  Output should be a valid CCHF-formatted file.
%
%list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/Donyian/' ...
%    'observations/Donyian_streamflow_measurement_info/*.xlsx']);
%list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/Langtang/' ...
%    'observations/langtang_streamgage/*DHM*.csv']);
%list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/' ...
%    'AM_Vakhsh_at_Komsomolabad/observations/' ...
%    'AM_Vakhsh_at_Komsomolabad.day_discharge.dat']);
list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/' ...
    'SY_Naryn_at_NarynTown/observations/' ...
    'SY_Naryn_at_NarynTown.day_discharge.dat']);

inFileName = fullfile(list.folder, list.name);
fprintf(inFileName);
% Vakhsh
%lon = linspace(69.5, 70.5, 10);
%lat = linspace(38.5, 39.5, 10);

% Naryn
lon = linspace(75.5, 76.5, 10);
lat = linspace(41.0, 42.0, 10);

sObs = read_gagedata(inFileName, lon, lat)

fprintf("test");

 
