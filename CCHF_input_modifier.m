%Input modifier. 
%This script loads a previously generated CCHF input file (e.g.
%CCHF_saved_input_paths.mat) and modifies it for a different set of
%conditions. 
%It is principally designed for a different set of climate
%conditions.


%There are two cell arrays in the input file: 'pathGage' and 'sPath'
%These cell arrays have the same length as the number of regions. The
%contents are defined as:
    %pathGage{i} = path to CCHF formatted observation file (only used in
        %calibration, validation, or default run types; not used in 
        %simulation run type).
    %sPath{i} = structure array with fields:
        %'dem' = path to DEM file (defined externally, e.g. ArcMap)
        %'fdr' = path to FDR file (defined externally, e.g. ArcMap)
        %climate var ('pr', 'tas', 'tasmin', 'tasmax', 'bcdep', etc.) = 
            %the set of variables depends on model formulation; this points
            %to folder containing time-series files for each variable
        %[clim var 'File'] = the set of climate variable fields here 
            %matches the variables being used. It is a cell array of 
            %specific cliamte files where the number of entries equals 
            %points to the climate file to be loaded in each time step 
            %(this is more efficient for parallel instances); DELETE FIELD
        %'coef' = path to file containing coefficients to use in CCHF
            %model; DELETE FIELD (there is path input in CCHF main script 
            %that points to this)
        %'ice' = path to shapefile defining glacier outlines
        %'icwe' = path to file defining glacier initial water equivalent
        %'icdbr' = path to file defining glacier debris cover thickness
        %'icpndx' = path to file defining glacier pond fraction
        %'regclip' = obsolete; DELETE FIELD
        %'outputMain' = Overarching output folder for current model run; DELETE FIELD
        %'output' = output folder for current region within model run; DELETE FIELD
        
pathLd = '/Users/thomas/Downloads/Hunza_test/model_runs/CCHF_saved_input_paths_4_Hunza-bcdep.mat';
pathSv = '/Users/thomas/Downloads/Hunza_test/model_runs/CCHF_saved_input_paths_4_Hunza-bcdep_test.mat';

%NEW CLIMATE VARIABLE FIELDS
%These are cell arrays of same length as already present in 'sPath'
%These should all be folders (i.e. should not point to specific climate
%file within folder)
pathPr     = {'test'};
pathTas    = {'test'};
pathTasmax = {'test'};
pathTasmin = {'test'};
pathBcdep     = {'test'};

%%Processing
load(pathLd, 'pathGage', 'sPath');

%Check correct length of inputs:
if ~isequal(numel(sPath(:)), numel(pathPr(:)))
   error('cchfInputMod:diffLength','There is a length mismatch between input structure and user defied inputs.');
end

clmVar = {'pr', 'tas', 'tasmax', 'tasmin', 'bcdep'};

%Modify 'sPath':
for ii = 1 : numel(sPath(:))
    %Loop over possible climate variables
    for jj = 1 : numel(clmVar(:))
        %Change climate variable path
        if isfield(sPath{ii}, clmVar{jj})
            switch clmVar{jj}
                case 'pr'
                    sPath{ii}.(clmVar{jj}) = pathPr{ii};
                case 'tas'
                    sPath{ii}.(clmVar{jj}) = pathTas{ii};
                case 'tasmax'
                    sPath{ii}.(clmVar{jj}) = pathTasmax{ii};
                case 'tasmin'
                    sPath{ii}.(clmVar{jj}) = pathTasmin{ii};
                case 'bcdep'
                    sPath{ii}.(clmVar{jj}) = pathBcdep{ii};
                otherwise
                    warning('cchfInputMod:unknownClmVar', ...
                        ['Climate variable ' clmVar{jj} ' has not been programmed for.'])
            end
        end
        %Remove file pointers (these will be redefined in CCHF)
        varClmFile = [clmVar{jj} 'File'];
        if isfield(sPath{ii}, varClmFile)
            sPath{ii} = rmfield_x(sPath{ii}, varClmFile);
        end
    end
    
    %Remove additional fields (these will be redefined in CCHF)
    sPath{ii} = rmfield_x(sPath{ii}, {'regclip', 'outputMain', 'output', 'coef'});
end

%Save input file:
save(pathSv, 'pathGage', 'sPath');
