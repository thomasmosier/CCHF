%Add constant factors to climate data.

scen = 'rcp85';



foldInput = uigetdir(pwd, 'Select the folder containing gridded NetCDF data to plot.');


foldOut = fullfile(foldInput,scen);
    mkdir(foldOut);

%Find files to edit
filesEdit = dir([foldInput, filesep, '*', '.nc']);
    filesEdit = extractfield(filesEdit,'name');

%Find Variable:
finfo = ncinfo(fullfile(foldInput, filesEdit{1}));
allVar = extractfield(finfo.Variables,'Name');
    allVar(strcmpi(allVar,'lon')) = [];
    allVar(strcmpi(allVar,'longitude')) = [];
    allVar(strcmpi(allVar,'lat')) = [];
    allVar(strcmpi(allVar,'latitude')) = [];
    allVar(strcmpi(allVar,'time')) = [];
    allVar(strcmpi(allVar,'time_bnds')) = [];
if numel(allVar) == 1
    varLd = allVar{1};
else
    error('plot_gridded:autoVarFail','Cannot auto-detect variable to load.');
end
    

%climate perturbation scenarios:
if regexpbl(scen,'RCP45') %Rcp4.5: +2.5c and +12% precip
    dTas = 2.5; %Deg C
    dPre = 1.12; %20% increase
elseif regexpbl(scen,'RCP85') %Rcp8.5: +4.3c and +21% precip
    dTas = 4.3; %Deg C
    dPre = 1.21; %20% increase
else
    error('clim_perturb:scen',['Scenario ' scen ' has not been programmed for.']);
end

%Select perturbation based on variable:
if regexpbl(varLd,'pr')
    scale = 'mult';
    pert = dPre;
    prec = 0;
elseif regexpbl(varLd,{'tas','tmn','tmin','tmx','tmax'})
    scale = 'add';
    pert = dTas;
    prec = 1;
else
    error('clim_pertubrb:perturbVar',['A perturbation has not been defined for ' varLd '.']);
end

%Open waitbar:
hWait = waitbar(0,'Processing starting.');  

%Loop over all files:
for ii = 1 : numel(filesEdit)
        waitbar(ii/numel(filesEdit), hWait, ...
            ['Processing file ' num2str(ii) ' of ' num2str(numel(filesEdit)) '.']);
    
    
    pathCurr = fullfile(foldInput,filesEdit{ii});
    
    %READ FIELDS:
    varTemp = ncinfo(pathCurr); 
    varLp = squeeze(struct2cell(varTemp.Variables));
        varLp = varLp(1,:)';
    
    sCurr = struct;
    for jj = 1 : numel(varLp)
        tempAtt = ncinfo(pathCurr, varLp{jj});
        if ~isempty(tempAtt.Attributes)
            sCurr.att.(varLp{jj}) = squeeze(struct2cell(tempAtt.Attributes))';
        else
            sCurr.att.(varLp{jj}) = cell(0,0);
        end
%             sCurr.attData = squeeze(struct2cell(sCurr.att.(varLp{jj}).Attributes))';
%         
        sCurr.var.(varLp{jj}) = ncread(pathCurr, varLp{jj});
    end
    
    if regexpbl(scale,'add')
        sCurr.var.(varLd) = sCurr.var.(varLd) + pert;
    elseif regexpbl(scale,'mul')
        sCurr.var.(varLd) = sCurr.var.(varLd)*pert;
    else
        error('clim_perturb:scaleNotKnown',['The scale form ' scale ' has not been programmed for.']);
    end
    
    %WRITE PERTURBED NetCDF FILE:
    write_NC_specify_att(fullfile(foldOut,filesEdit{ii}), sCurr, prec);
end

close(hWait);




