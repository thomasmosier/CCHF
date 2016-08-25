clear all


pathDEM = 'E:\ARCGIS\Falls_Creek\manual\dive_dem_wc.asc';
pathFDR = 'E:\ARCGIS\Falls_Creek\manual\dive_fdr_wc.asc';
    ptInterest = [-122.34986,  44.37944];
    nm = 'Diversion';
    sMeta.region = 'FC';

% pathDEM = 'E:\ARCGIS\Cook_Inlet\Wolverine\edited\wolverine_dem_wc.asc';
% pathFDR = 'E:\ARCGIS\Cook_Inlet\Wolverine\edited\wolverine_fdr_manual.asc';
%     ptInterest = [-148.896670,  60.370556];
%     nm = '15236900';
%     sMeta.region = 'Wolverine';

% pathDEM = 'E:\ARCGIS\Cook_Inlet\Gulkana\edited\gulkana_dem_wc.asc';
% pathFDR = 'E:\ARCGIS\Cook_Inlet\Gulkana\edited\gulkana_fdr_wc.asc';
%     ptInterest = [-145.467500,  63.240833];
%     nm = '15478040';
%     sMeta.region = 'Gulkana';


%Load DEM and calculate watershed geometry fields:
sFDR = read_geodata(pathFDR, sMeta, 'none');
    sHydro.lat = sFDR.lat;
    sHydro.lon = sFDR.lon;
    sHydro.fdrESRI = sFDR.data;
sDem = read_geodata(pathDEM, sMeta, 'none');
    sHydro.dem = sDem.data;
clear('sDem', 'sFDR');
%Calculate fdr, fac, flow calculation order, slope, aspect, distance between neighboring centroids
sHydro = watershed(sHydro);
% 
% [R, S] = dem_flow(sHydro.dem);
% T = flow_matrix(sHydro.dem, R);

[area, fdrSteps, totFdr] = area_upslope(sHydro);

[fold, file, ~] = fileparts(pathDEM);
pathFAC = fullfile(fold, [file '_fac.asc']);
hdr = ESRI_hdr(sHydro.lon, sHydro.lat, 'corner');
write_ESRI_v4(area, hdr, pathFAC, 0);

[~, gageRow] = min(abs(ptInterest(2)-sHydro.lat));
[~, gageCol] = min(abs(ptInterest(1)-sHydro.lon));

indI = round(sub2ind(size(sHydro.fdrESRI), gageRow, gageCol));

shed = reshape(totFdr(:,indI),size(sHydro.dem));
shed(shed == 0) = nan;
pathShed = fullfile(fold, [file '_shed' nm '.asc']);
write_ESRI_v4(shed, hdr, pathShed, 0);

