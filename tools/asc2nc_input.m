%
%MUST UPDATE CODE RELATED TO 'timeVec' BASED ON FORMAT OF ESRI ASCIIs
%

extInput = 'asc';
varOut = 'pr'; %Variable name of output. Must select to be consistent with 
        %CCHF expectations (typically follows CMIP5 naming convention).
prec = 3; %Precision of output files: 
                %-1 = unsigned integer, 
                %0 = signed integer, 
                %1:5 = single, 
                %>5 = double  
        
                
dirInput = uigetdir(pwd,['Select the folder containing '...
        'the homogeneous ascii grid files.']);
    
filesInput = find_files(dirInput,extInput);

dirOutput = fullfile(dirInput,'nc');
mkdir(dirOutput)

for ii = 1 : numel(filesInput)
    pathCurr = fullfile(dirInput, filesInput{ii});
    [~, fileCurr, ~] = fileparts(pathCurr);

    [dataCurr, hdrESRI, hdrMeta] = read_ESRI(pathCurr);
    [lat, lon] = ESRI_hdr2geo(hdrESRI, hdrMeta);
    
    pathOut = fullfile(dirOutput, [fileCurr, '.nc']);
    
    %Populate 'timeVec'
        %Format = [year, month, day] (It can also include hour as well if 
        %that's the resolution of your data...)
    %For example: This code will work if the file naming format is 'root_year_month_day.ext'
    indUnd = regexpi(fileCurr,'_');
    indPer = regexpi(fileCurr,'\.');
    timeVec = [str2double(fileCurr(indUnd(end-2)+1:indUnd(end-1)-1)), ...
        str2double(fileCurr(indUnd(end-1)+1:indUnd(end)-1)), ...
        str2double(fileCurr(indUnd(end)+1:indPer(end)-1))];
    
    print_grid_NC(pathOut, dataCurr, varOut, lon, lat, timeVec, prec)
end