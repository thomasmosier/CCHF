nmRegion = 'Karora';
latOut = [];
lonOut = [];
yrsOut = [];

%Note: The DEM grid is used to set the grid for all other parameters
%(through resampling if they differ)


%%Mandatory inputs (paths must have file):
%DEM @ desired model res (file path):
pathInDem = '';
%Precipitation (path to individual file in folder; other files in folder will be searched, too):
pathInPre = '';
%Mean temperature (path to individual file in folder; other files in folder will be searched, too):
pathInTmp = '';

%%Optional Inputs:
%Glaciers:

%FDR  (file path):
pathInFdr = '';
%Debris  (file path):
pathInDebris = '';
%Glacier Lakes  (file path):
pathInGlacLakes = '';
%Glacier Thickness  (file path):
pathInGlacWe = '';

%MODIS SCA:

%Maximum Temperature (folder):
pathInTmx = '';
%Minimum Temperature (folder):
pathInTmn = '';
%Even higher resolution DEM (can be used for glacier processes):
pathInGlacDem = '';



%%Options used in processing:
fr = 0;
crop = 'in';
precPre = '1';
varTime = 'time';
varDate = 'date';


%%Implementation:

%Create output folder:
[foldDem, ~, ~] = fileparts(pathInDem);
foldOut = fullfile(foldDem, nmRegion);
if ~exist(foldOut, 'dir')
   mkdir(foldOut);
end


%Process DEM:
varDem = 'z';
sDem = read_geodata_v2(pathInDem, varDem, lonOut, latOut, nan(1,2), fr, crop);
pathOutDem = fullfile(foldOut, [nmRegion '_dem.asc']);
write_geodata_v2(pathOutDem, sDem, 0, 'asc');

%Process Flow direction grid:
if exist(pathInFdr, 'file')
    varFdr = 'fdr';
    sFdr = read_geodata_v2(pathInFdr, varFdr, lonOut, latOut, nan(1,2), fr, crop);
    pathOutFdr = fullfile(foldOut, [nmRegion '_fdr.asc']);
    write_geodata_v2(pathOutFdr, sFdr, 0, 'asc');
end

%Process debris cover:
if exist(pathInDebris, 'file')
    varDebris = 'debris';
    sDebris = read_geodata_v2(pathInDebris, varDebris, lonOut, latOut, nan(1,2), fr, crop);
    pathOutDebris = fullfile(foldOut, [nmRegion '_debris_thickness.asc']);
    if isequal(round2(sDebris.(varLon), 5), round2(sDem.(varLon), 5)) && isequal(round2(sDebris.(varLat), 5), round2(sDem.(varLat), 5)) 
        debrisOut = sDebris.(varDebris)(indCurr, :, :);
    else
        sDebris.(varDebris) = geodata_area_wgt(sDebris.(varLon), sDebris.(varLat), sDebris.(varDebris)(indCurr, :, :), sDem.(varLon), sDem.(varLat));
        sDebris.(varLon) = sDem.(varLon);
        sDebris.(varLat) = sDem.(varLat);
    end
    write_geodata_v2(pathOutDebris, sDebris, 2, 'asc');
end

%Process glacier lakes:
if exist(pathInGlacLakes, 'file')
    varGlacLakes = 'water';
    sGlacLakes = read_geodata_v2(pathInGlacLakes, varGlacLakes, lonOut, latOut, nan(1,2), fr, crop);
    pathOutGlacLakes = fullfile(foldOut, [nmRegion '_glacier_lakes.asc']);
    if isequal(round2(sGlacLakes.(varLon), 5), round2(sDem.(varLon), 5)) && isequal(round2(sGlacLakes.(varLat), 5), round2(sDem.(varLat), 5)) 
        debrisOut = sGlacLakes.(varGlacLakes)(indCurr, :, :);
    else
        sGlacLakes.(varGlacLakes) = geodata_area_wgt(sGlacLakes.(varLon), sGlacLakes.(varLat), sGlacLakes.(varGlacLakes)(indCurr, :, :), sDem.(varLon), sDem.(varLat));
        sGlacLakes.(varLon) = sDem.(varLon);
        sGlacLakes.(varLat) = sDem.(varLat);
    end
    write_geodata_v2(pathOutGlacLakes, sGlacLakes, 2, 'asc');
end

%Process glacier thickness (water equivalent):
if exist(pathInGlacWE, 'file')
    varGlacWE = 'we';
    sGlacWE = read_geodata_v2(pathInGlacWE, varGlacWE, lonOut, latOut, nan(1,2), fr, crop);
    pathOutGlacWE = fullfile(foldOut, [nmRegion '_glacier_thickness_we.asc']);
    if isequal(round2(sGlacWE.(varLon), 5), round2(sDem.(varLon), 5)) && isequal(round2(sGlacWE.(varLat), 5), round2(sDem.(varLat), 5)) 
        debrisOut = sGlacWE.(varGlacWE)(indCurr, :, :);
    else
        sGlacWE.(varGlacWE) = geodata_area_wgt(sGlacWE.(varLon), sGlacWE.(varLat), sGlacWE.(varGlacWE)(indCurr, :, :), sDem.(varLon), sDem.(varLat));
        sGlacWE.(varLon) = sDem.(varLon);
        sGlacWE.(varLat) = sDem.(varLat);
    end
    write_geodata_v2(pathOutGlacWE, sGlacWE, 2, 'asc');
end


%Process Precipitation:
varPre = 'pre';
sPre = read_geodata_v2(pathInPre, varPre, lonOut, latOut, yrsOut, fr, crop);

dateMnths = unique(sPre.(varDate)(:,2:3));
foldOutPre = fullfile(foldOut, 'climate', 'pre');
if ~exist(foldOutPre, 'dir')
   mkdir(foldPre) 
end
for ii = 1 : numel(dateMnths)
    [indCurr, ~] = ismember(sPre.(varDate)(:,2:3), dateMnths(ii,:), 'rows');
    nDy = sum(indCurr);
        
    pathOutPre = fullfile(foldPreOut, [varPre '_' nmRegion '_' num2str(dateMnths(ii,1)) '_' num2str(dateMnths(ii,2)) '.nc']);
    %If file already exists, delete:
    if exist(pathOutPre, 'file')
        delete(pathOutPre);
    end
    %Create vector of output times:
    timeVec = sPre.date(indCurr,:);
    
    if isequal(round2(sPre.(varLon), 5), round2(sDem.(varLon), 5)) && isequal(round2(sPre.(varLat), 5), round2(sDem.(varLat), 5)) 
        preOut = sPre.(varPre)(indCurr, :, :);
    else
        preOut = geodata_area_wgt(sPre.(varLon), sPre.(varLat), sPre.(varPre)(indCurr, :, :), sDem.(varLon), sDem.(varLat));
    end
    %Print data to file:
    print_grid_NC_v2(pathOutPre, preOut, sDem.(varLon), sDem.(varLat), ...
        timeVec, timeVec, unitsOut, precPre);
end


%Process Mean Temperature:
varTmp = 'tas';
sTmp = read_geodata_v2(pathInTmp, varTmp, lonOut, latOut, yrsOut, fr, crop);

dateMnths = unique(sTmp.(varDate)(:,2:3));
foldOutTmp = fullfile(foldOut, 'climate', 'pre');
if ~exist(foldOutTmp, 'dir')
   mkdir(foldTmp) 
end
for ii = 1 : numel(dateMnths)
    [indCurr, ~] = ismember(sTmp.(varDate)(:,2:3), dateMnths(ii,:), 'rows');
    nDy = sum(indCurr);
        
    pathOutTmp = fullfile(foldTmpOut, [varTmp '_' nmRegion '_' num2str(dateMnths(ii,1)) '_' num2str(dateMnths(ii,2)) '.nc']);
    %If file already exists, delete:
    if exist(pathOutTmp, 'file')
        delete(pathOutTmp);
    end
    %Create vector of output times:
    timeVec = sTmp.date(indCurr,:);
    
    if isequal(round2(sTmp.(varLon), 5), round2(sDem.(varLon), 5)) && isequal(round2(sTmp.(varLat), 5), round2(sDem.(varLat), 5)) 
        preOut = sTmp.(varTmp)(indCurr, :, :);
    else
        preOut = geodata_area_wgt(sTmp.(varLon), sTmp.(varLat), sTmp.(varTmp)(indCurr, :, :), sDem.(varLon), sDem.(varLat));
    end
    %Print data to file:
    print_grid_NC_v2(pathOutTmp, preOut, sDem.(varLon), sDem.(varLat), ...
        timeVec, timeVec, unitsOut, precTmp);
end


%Process Minimum Temperature:
if exist(pathInTmn, 'file')
    varTmn = 'tasmin';
    sTmn = read_geodata_v2(pathInTmn, varTmn, lonOut, latOut, yrsOut, fr, crop);

    dateMnths = unique(sTmn.(varDate)(:,2:3));
    foldOutTmn = fullfile(foldOut, 'climate', 'pre');
    if ~exist(foldOutTmn, 'dir')
       mkdir(foldTmn) 
    end
    for ii = 1 : numel(dateMnths)
        [indCurr, ~] = ismember(sTmn.(varDate)(:,2:3), dateMnths(ii,:), 'rows');
        nDy = sum(indCurr);

        pathOutTmn = fullfile(foldTmnOut, [varTmn '_' nmRegion '_' num2str(dateMnths(ii,1)) '_' num2str(dateMnths(ii,2)) '.nc']);
        %If file already exists, delete:
        if exist(pathOutTmn, 'file')
            delete(pathOutTmn);
        end
        %Create vector of output times:
        timeVec = sTmn.date(indCurr,:);

        if isequal(round2(sTmn.(varLon), 5), round2(sDem.(varLon), 5)) && isequal(round2(sTmn.(varLat), 5), round2(sDem.(varLat), 5)) 
            preOut = sTmn.(varTmn)(indCurr, :, :);
        else
            preOut = geodata_area_wgt(sTmn.(varLon), sTmn.(varLat), sTmn.(varTmn)(indCurr, :, :), sDem.(varLon), sDem.(varLat));
        end
        %Print data to file:
        print_grid_NC_v2(pathOutTmn, preOut, sDem.(varLon), sDem.(varLat), ...
            timeVec, timeVec, unitsOut, precTmn);
    end
end

%Process Maximum Temperature:
if exist(pathInTmn, 'file')
    varTmx = 'tasmax';
    sTmx = read_geodata_v2(pathInTmx, varTmx, lonOut, latOut, yrsOut, fr, crop);

    dateMnths = unique(sTmx.(varDate)(:,2:3));
    foldOutTmx = fullfile(foldOut, 'climate', 'pre');
    if ~exist(foldOutTmx, 'dir')
       mkdir(foldTmx) 
    end
    for ii = 1 : numel(dateMnths)
        [indCurr, ~] = ismember(sTmx.(varDate)(:,2:3), dateMnths(ii,:), 'rows');
        nDy = sum(indCurr);

        pathOutTmx = fullfile(foldTmxOut, [varTmx '_' nmRegion '_' num2str(dateMnths(ii,1)) '_' num2str(dateMnths(ii,2)) '.nc']);
        %If file already exists, delete:
        if exist(pathOutTmx, 'file')
            delete(pathOutTmx);
        end
        %Create vector of output times:
        timeVec = sTmx.date(indCurr,:);

        if isequal(round2(sTmx.(varLon), 5), round2(sDem.(varLon), 5)) && isequal(round2(sTmx.(varLat), 5), round2(sDem.(varLat), 5)) 
            preOut = sTmx.(varTmx)(indCurr, :, :);
        else
            preOut = geodata_area_wgt(sTmx.(varLon), sTmx.(varLat), sTmx.(varTmx)(indCurr, :, :), sDem.(varLon), sDem.(varLat));
        end
        %Print data to file:
        print_grid_NC_v2(pathOutTmx, preOut, sDem.(varLon), sDem.(varLat), ...
            timeVec, timeVec, unitsOut, precTmx);
    end
end

