%
%test script for new DHM format streamgage data
%
%list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/Donyian/' ...
%    'observations/Donyian_streamflow_measurement_info/*.xlsx']);
list = dir(['/pl/active/PMESDR/CCHF_data/assessment_sites/Langtang/' ...
    'observations/langtang_streamgage/*DHM*.csv']);

inFileName = fullfile(list.folder, list.name);
unitsIn = 'm3/s';
fprintf(inFileName);


[data, date] = read_ts_table_data(inFileName, unitsIn);
fprintf("test");

 
