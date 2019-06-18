%
%test script for new DHM format streamgage data at
%read_gagedata level.  Output should be a valid CCHF-formatted file.
%
%list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/Donyian/' ...
%    'observations/Donyian_streamflow_measurement_info/*.xlsx']);
list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/Langtang/' ...
    'observations/langtang_streamgage/*DHM*.csv']);

inFileName = fullfile(list.folder, list.name);
fprintf(inFileName);
lon = linspace(85., 86., 10)
lat = linspace(28., 29., 10);

sObs = read_gagedata(inFileName, lon, lat)

fprintf("test");

 
