% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function [varargout] = CCHF_modules(sHydro, sMeta)
global sLand sAtm sCryo


%List of models that do not use top-of-atmosphere radiation, transmissivity, etc.
lsSimpleMod = {'simple','hock', 'SDI', 'kraaijenbrink'};
blSimpleMod = regexpbl(find_att(sMeta.module, 'heat', 'no_warning'), lsSimpleMod);


%If parameter mode, initialize output array
if regexpbl(sMeta.mode,'parameter')
   coef = cell(0,6); 
end

%Set current model run dates (used for solar raad and PET)
datesUse = sMeta.dateRun;
if iscell(datesUse)
    if isfield(sMeta, 'siteCurr')
        datesUse = datesUse{sMeta.siteCurr};
    else
        error('CchfModules:noSiteCurr', ['The date field is a cellarray. '...
            'This requires the presence of a siteCurr field to determine which index to use.']);
    end
end


%Initialize ice thickness on first iteration:
if ~regexpbl(sMeta.mode,'parameter') && sMeta.indCurr == 1
    if ~strcmpi(sMeta.iceGrid, 'none') 
        iceWEMod = find_att(sMeta.module, 'glacier0', 'no_warning');
        if regexpbl(iceWEMod, 'Shea')
            %Estimate thickness (force balance based on equilibrium
            %assumption from Shea et al. 2015)
            glacier0_Shea(sMeta);
        elseif regexpbl(iceWEMod, 'Chen')
            %Estimate thickness (empirical surface area scaling based on Chen and Ohmura, 1990)
            glacier0_Chen(sHydro, sMeta);
        elseif regexpbl(iceWEMod, 'external')
            %Set ice thickness using data from file
            glacier0_external(sHydro, sMeta);
        else
            error('cchfModules:iceWE', ['The ice water equivalent initialization ' iceWEMod ' is not recognized.']);
        end 
    end
end


%Track glacier mass balance
%Initialize grid to track glacier mass balance:
if sMeta.indCurr == 1
    sCryo.icmb = zeros(size(sHydro.dem),'single');

    %Display glacier mass balance date on first iteration
    if  ~regexpbl(sMeta.mode, 'calib')
        disp(['Calculating glacier dynamics each year of model run on ' ...
            num2str(sMeta.dateGlac(1)) '/' num2str(sMeta.dateGlac(2)) '.']);
    end
end


%Reset mass balance on same day of the year at model start date
%Note: mass balance is updated daily at the bottom of this script
%Note: glacier mb calculation dates are based on model start date
if isequal(sMeta.dateCurr(2:end), sMeta.dateStart{sMeta.siteCurr}(2:end))
    sCryo.icmb = zeros(size(sHydro.dem),'single');
    %Record glacier mass balance start date
    sCryo.icmb_dateStart = sMeta.dateCurr;
    sCryo.icmb_dateEnd = nan(size(sMeta.dateCurr));
end
%See recording of glacier mass balance at the end of this subroutine


%PARTITION PRECIPITATION INTO RAIN AND SNOWFALL:
partMod = find_att(sMeta.module, 'partition', 'no_warning');
if regexpbl(partMod, 'ramp') 
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef,partition_ramp());
    else
        partition_ramp(sMeta); %Partitions precipitation to snowfall and rain.
    end
else
    error('cchfModules:partitionUnknown', ['The precipitation partitioning representation ' partMod ' is not recognized.']);
end


%SNOWPACK ACCUMULATION (ADD SNOWFALL AND RAIN ON SNOW)
%Adds snowfall to snowpack and rain to liquid snow content
if ~regexpbl(sMeta.mode,'parameter')
    snow_accum;
end


%SNOW ALBEDO:
if ~blSimpleMod %Albedo not used in simple degree-index model
    snalbMod = find_att(sMeta.module, 'snalbedo', 'no_warning');
    if regexpbl(snalbMod, 'Pelli')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, snalbedo_Pellicciotti());
        else
            snalbedo_Pellicciotti(sMeta);
        end
    elseif regexpbl(snalbMod, 'Brock')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, snalbedo_Brock());
        else
            snalbedo_Brock(sMeta);
        end
    else
        error('cchfModules:snalbedoUnknown',['Ice albedo representation ' snalbMod ' not recognized.']);
    end
end


%ICE ALBEDO
if ~strcmpi(sMeta.iceGrid, 'none') && ~blSimpleMod
    icalbMod = find_att(sMeta.module, 'icalbedo', 'no_warning');
    if regexpbl(icalbMod, 'constant')
        %If debris grid included, set albedo to 0.15 at debris covered ice locations:
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, icalbedo_constant());
        else
            icalbedo_constant(sMeta);
        end
    elseif ~regexpbl(icalbMod, 'none')
        error('cchfModules:icalbedoUnknown', ['The ice albedo representation ' icalbMod ' is not recognized.']);
    end
end


%BLACK CARBON IMPACT ON ALBEDO:
%NOTE: THIS MUST COME AFTER SNOW AND ICE ALBEDO CALCULATIONS AND BEFORE
%HEAT CALCULATION
if ~blSimpleMod
    bcalbMod = find_att(sMeta.module, 'bcalbedo', 'no_warning');
    if isempty(bcalbMod)
        bcalbMod = 'none';
        sMeta.module(end+1,:) = {'bcalbedo', 'none'}; 
        if sMeta.indCurr == 1
            warning('CCHFModules:noBcModule',['No black carbon process '...
                'representation was included. Therefore it is being set to ' char(39) 'none' char(39) '.']);
        end
    end
    
    if regexpbl(bcalbMod, 'Ming')
        %If debris grid included, set albedo to 0.15 at debris covered ice locations:
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, bcalbedo_Ming());
        else
            bcalbedo_Ming(sMeta);
        end
    elseif ~regexpbl(bcalbMod, 'none')
        error('cchfModules:bcalbedoUnknown', ...
            ['The black carbon albedo representation ' bcalbMod ...
            ' is not recognized.']);
    end
end


%TOP OF ATMOSPHERE SOLAR RADIATION (not used in simple degree index):
if ~blSimpleMod
    toaMod = find_att(sMeta.module,'toa', 'no_warning');

    if regexpbl(toaMod, 'DeWalle')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, toa_rad_DeWalle([], sHydro));
        else
            if ~isfield(sAtm,'rsdt')
                toa_rad_DeWalle(datesUse, sHydro, sMeta);
            end
        end
    elseif regexpbl(toaMod,'Liston')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, toa_rad_Liston([], sHydro));
        else
            if ~isfield(sAtm,'rsdt')
                toa_rad_Liston(datesUse, sHydro, sMeta);
            end
        end
        %
        %EDIT HERE
        %ADD A NEW RADIATION MODEL!
        %
%     elseif regexpbl(toaMod,'Pritchard')
%         if regexpbl(sMeta.mode,'parameter')
%             coef = cat(1,coef, toa_rad_Pritchard([], sHydro));
%         else
%             if ~isfield(sLand,'rsdt')
%                 toa_rad_Pritchard(sMeta.dateRun, sHydro, sMeta);
%             end
%         end
    else
        error('cchfModules:toaUnknown', ['Top-of-Atmosphere radiation representation ' toaMod ' not recognized.']);
    end
end



%TRANSMISSIVITY (not used in simple degree index):
if ~blSimpleMod
    transMod = find_att(sMeta.module,'atmtrans', 'no_warning');
    if regexpbl(transMod, 'DeWalle')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, atmtrans_DeWalle());
        else
            atmtrans_DeWalle(sMeta);
        end
    elseif regexpbl(transMod, 'Coops')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, atmtrans_Coops(sHydro));
        else
            atmtrans_Coops(sHydro,sMeta);
        end
    elseif regexpbl(transMod, 'dem_decay')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, atmtrans_dem_decay(sHydro));
        else
            atmtrans_dem_decay(sHydro,sMeta);
        end
        elseif regexpbl(transMod, 'Gauss')
            if regexpbl(sMeta.mode,'parameter')
                coef = cat(1,coef, atmtrans_Gauss(sHydro));
            else
                atmtrans_Gauss(sHydro,sMeta);
            end
    else
        error('cchfModules:transmissivityUnknown',['Atmospheric transmissivity representation ' transMod ' not recognized.']);
    end
end

    

%CRYOSPHERE HEAT:
heatMod = find_att(sMeta.module, 'heat', 'no_warning');
if strcmpi(heatMod, 'SDI') || strcmpi(heatMod, 'STI') || strcmpi(heatMod, 'simple')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_simple_degree());
    else
        heat_simple_degree(sMeta);
    end
elseif regexpbl(heatMod, 'Pelli')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_ETI_Pellicciotti());
    else
        heat_ETI_Pellicciotti(sMeta);
    end
elseif regexpbl(heatMod, 'Hock')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_ETI_Hock());
    else
        heat_ETI_Hock(sMeta);
    end
elseif regexpbl(heatMod, {'LSD', 'LST'})
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_LST_Mosier());
    else
        heat_LST_Mosier(sMeta)
    end
elseif regexpbl(heatMod, {'TI','Kraaijenbrink'},'and')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_TI_Kraaijenbrink());
    else
        heat_TI_Kraaijenbrink(sMeta)
    end
elseif regexpbl(heatMod, {'ETI','debris','Reid'},'and')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_ETI_debris_Reid());
    else
        heat_ETI_debris_Reid(sMeta)
    end
elseif regexpbl(heatMod, {'ETI','debris','emp'},'and')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_ETI_debris_emp());
    else
        heat_ETI_debris_emp(sMeta)
    end
elseif regexpbl(heatMod,'Mosier')
%         if regexpbl(sMeta.mode,'parameter')
%             coef = cat(1,coef, ETI_Mosier());
%         else
%             %Calculate melt as being proportional to modelled shortwave radiation and temperature:
%             heat_ETI_Mosier(sMeta);
%         end
%     elseif regexpbl(sMeta.module,'energy-SETIbasic')
%         if regexpbl(sMeta.mode,'parameter')
%             coef = cat(1,coef, heat_SETI_basic(sHydro));
%         else
%             %Heat transfer:
%             heat_SETI_basic(sHydro,sMeta);
%         end
elseif regexpbl(heatMod, 'SETI')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_SETI_Mosier());
    else
        heat_SETI_Mosier(sMeta);
    end
else
    error('cchfModules:heatUnknown',['Cryosphere heat representation ' heatMod ' not recognized.']);
end



%SNOWPACK DENSITY:
%Comes after heat flux because depends on cryosphere surface temperature
if regexpbl(find_att(sMeta.module, 'heat', 'no_warning'), 'SETI') || regexpbl(find_att(sMeta.module, 'sndrain', 'no_warning'), 'Liston')%Albedo not used in simple degree-index model
    dnsMod = find_att(sMeta.module, 'density', 'no_warning');
    
    if regexpbl(dnsMod, 'Liston')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, density_Liston());
        else
            density_Liston(sMeta);
        end
    elseif regexpbl(dnsMod, 'Cosima')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, density_Cosima());
        else
            density_Cosima(sMeta);
        end
    elseif ~regexpbl(icalbMod, 'none')
        error('cchfModules:densityUnknown',['Snow density representation ' dnsMod ' not recognized.']);
    end
end



%SNOW ACCUMULATION AND MELT:
massMod = find_att(sMeta.module, 'snmass', 'no_warning');
if regexpbl(massMod, 'step') %Calculate snow and glacier melt using heat threshold
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, snmass_step());
    else
        snmass_step(sMeta);
    end
elseif regexpbl(massMod, 'cc')
    if regexpbl(sMeta.mode,'parameter') %Calculate snow and melt using CC
        coef = cat(1,coef, snmass_cc_v2());
    else
        snmass_cc_v2(sMeta);
    end
elseif regexpbl(massMod, {'enbal', 'energy'}) %Conduction raises and lowers snow temperature
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, snmass_enbal());
    else
        snmass_enbal(sMeta);
    end
elseif regexpbl(massMod, 'Liston') %Conduction raises and lowers snow temperature (no fitting parameters)
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, snmass_Liston());
    else
        snmass_Liston(sMeta);
    end
else
    error('cchfModules:snmassUnknown', ['Snow mass flux representation ' massMod ' not recognized.']);
end



%SNOW AVALANCHING (using glacier grid):
avMod = find_att(sMeta.module, 'avalanche', 'no_warning');
if isempty(avMod)
    avMod = 'none';
    sMeta.module(end+1,:) = {'avalanche', 'none'}; 
    if sMeta.indCurr == 1
        warning('CCHFModules:noAvalancheModule',['No avalanching process '...
            'representation was included. Therefore it is being set to ' char(39) 'none' char(39) '.']);
    end
end
if regexpbl(avMod, 'angle') 
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, avalanche_angle(sHydro));
    else
        avalanche_angle(sHydro, sMeta);
    end
    
elseif ~regexpbl(avMod, 'none')
    error('cchfModules:avalancheUnknown', ['The avalanching representation ' avMod ' not recognized.']);
end



%SNOW SUBLIMATION:
sublimateMod = find_att(sMeta.module,'sublimate', 'no_warning');
if isempty(sublimateMod)
    sublimateMod = 'none';
    sMeta.module(end+1,:) = {'sublimate', 'none'}; 
    if sMeta.indCurr == 1
        warning('CCHFModules:noSublimateModule',['No sublimation process '...
            'representation was included. Therefore it is being set to ' char(39) 'none' char(39) '.']);
    end
end
if regexpbl(sublimateMod, 'Lutz')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, sublimate_Lutz(sHydro));
    else
        sublimate_Lutz(sHydro, sMeta);
    end
elseif ~regexpbl(sublimateMod, 'none')
    error('cchfModules:sublimateUnknown', ['Sublimation representation ' sublimateMod ' not recognized.']);
end



%SNOW LIQUID WATER HOLDING CAPACITY:
liqMod = find_att(sMeta.module,'sndrain', 'no_warning');
if regexpbl(liqMod, 'percent') %Release melt water greater than percent holding capacity:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, sndrain_percent());
    else
        sndrain_percent(sMeta);
    end
elseif regexpbl(liqMod, 'density') %Release melt water greater than density threshold:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, sndrain_density());
    else
        sndrain_density(sMeta);
    end
else
    error('cchfModules:sndrainUnknown', ['The snow liquid holding representation ' liqMod ' not recognized.']);
end



%ICE MELT AND LIQUID RELEASE:
icMeltMod = find_att(sMeta.module,'icmlt', 'no_warning');
if regexpbl(icMeltMod,'ratio') %Remainder of energy from snow melt goes into melting ice
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_ratio());
    else
        icmlt_ratio(sMeta);
    end
elseif regexpbl(icMeltMod,'Liston') %Ice melt only occurs during time steps when no snow present
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_Liston());
    else
        icmlt_Liston(sMeta);
    end
elseif regexpbl(icMeltMod,'Neumann') %Ice melt only occurs during time steps when no snow present
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_Neumann_bc());
    else
        icmlt_Neumann_bc(sMeta);
    end
elseif regexpbl(icMeltMod,'time') %Ice melt only occurs during time steps when no snow present
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_time());
    else
        icmlt_time(sMeta);
    end
elseif regexpbl(icMeltMod, 'gradient') %Conduction raises and lowers snow temperature (no fitting parameters)
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_ablation_grad(sHydro));
    else
        icmlt_ablation_grad(sHydro, sMeta);
    end
else
    error('cchfModules:icmltUnknown', ['The ice melt and release representation ' icMeltMod ' not recognized.']);
end



%SNOW COVERED AREA:
scaMod = find_att(sMeta.module,'sca', 'no_warning');
if ~isempty(scaMod) && regexpbl(scaMod, 'max') %Assume snow uniformly distributed over grid cell:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, scx_max());
    else
        scx_max(sMeta);
    end
else
    error('cchfModules:scaUnknown', ['The snow covered area representation ' scaMod ' not recognized.']);
end   



%FIRN COMPACTION (SNOW BECOMES ICE):
if ~strcmpi(sMeta.iceGrid, 'none') 
    firnMod = find_att(sMeta.module,'firn', 'no_warning');
    if regexpbl(firnMod, 'threshold')
        if regexpbl(sMeta.mode,'parameter') 
            coef = cat(1,coef, firn_threshold());
        else
            firn_threshold(sMeta);
        end
    elseif ~regexpbl(firnMod, {'none', 'static'})
        error('cchfModules:flowUnknown', ['The flow routing representation ' firnMod ' not recognized.']);
    end
end



%PET:
petMod = find_att(sMeta.module,'pet', 'no_warning');
if regexpbl(petMod, 'Hammon') %Use Hammon formulation
    if regexpbl(sMeta.mode, 'parameter')
        coef = cat(1,coef, PET_Hammon(datesUse, sHydro));
    else
        if ~isfield(sLand, 'pet') || ~isfield(sLand, 'datepet')
            PET_Hammon(datesUse, sHydro, sMeta);
        end
    end
elseif regexpbl(petMod, 'Hargreaves') %Use Hargreaves formulation
    if regexpbl(sMeta.mode, 'parameter')
        coef = cat(1,coef, PET_Hargreaves(sMeta.dateCurr, sHydro));
    else
        PET_Hargreaves(sMeta.dateCurr, sHydro, sMeta);
    end
elseif regexpbl(petMod, 'Makkink') %Use Makkink formulation
    if regexpbl(sMeta.mode, 'parameter')
        coef = cat(1,coef, PET_Makkink(sMeta.dateCurr, sHydro));
    else
        PET_Makkink(sMeta.dateCurr, sHydro, sMeta);
    end
else
    warning('cchfModules:petUnknown', ['The PET representation ' petMod ' not recognized.']);
end



%RUNOFF:
runoffMod = find_att(sMeta.module, 'runoff', 'no_warning');
if regexpbl(runoffMod, 'bucket') %Calculate groundwater holding and release (using bucket model)
    if regexpbl(sMeta.mode, 'parameter')
        coef = cat(1,coef, runoff_bucket());
    else
        runoff_bucket(sMeta);
    end
elseif regexpbl(runoffMod, 'direct') %All water (precip and snowpack release) runs off immediately
    if regexpbl(sMeta.mode, 'parameter')
        coef = cat(1,coef, runoff_direct());
    else
        runoff_direct(sMeta);
    end
else
    error('cchfModules:runoffUnknown', ['The runoff generation representation ' runoffMod ' not recognized.']);
end



%TRAVEL TIME THROUGH CELLS:
timeMod = find_att(sMeta.module,'timelag', 'no_warning');
if regexpbl(timeMod, 'Johnstone') %Calculate time lag using Johnstone
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, tLag_Johnstone(sHydro));
    else
        if ~isfield(sLand,'tLag') %only calculate once because constant
            tLag_Johnstone(sHydro,sMeta);
        end
    end
elseif regexpbl(timeMod, 'Liston')
    if regexpbl(sMeta.mode,'parameter') %Fast and slow flow function of landcover
        coef = cat(1,coef, tLag_Liston(sHydro));
    else
        tLag_Liston(sHydro, sMeta);
    end
elseif ~regexpbl(timeMod, 'none')
    error('cchfModules:timelagUnknown', ['The flow travel time representation ' timeMod ' not recognized.']);
end



%ROUTING RUNOFF AND FLOW:
flowMod = find_att(sMeta.module,'flow', 'no_warning');
%Use numerical Muskingum method (allows dispersion)
if regexpbl(flowMod, 'Muskingum')
    if regexpbl(sMeta.mode,'parameter') 
        coef = cat(1,coef, flow_Muskingum(sHydro));
    else
        flow_Muskingum(sHydro, sMeta);
    end
elseif regexpbl(flowMod, 'lumped') %No dispersion of flow
    if regexpbl(sMeta.mode,'parameter') 
        coef = cat(1,coef, flow_lumped(sHydro));
    else
        flow_lumped(sHydro, sMeta);
    end
elseif regexpbl(flowMod, 'Liston') %Fast and slow flow
    if regexpbl(sMeta.mode,'parameter') 
        coef = cat(1,coef, flow_Liston(sHydro));
    else
        flow_Liston(sHydro, sMeta);
    end
elseif regexpbl(flowMod, 'none')
    %Set flow to nan if flow module not used
    sLand.flow(:) = nan;
elseif ~regexpbl(flowMod, 'none')
    error('cchfModules:flowUnknown', ['The flow routing representation ' flowMod ' not recognized.']);
end


%GLACIER PROCESSES (INCLUDES GLACIER DYNAMICS OPTIONS): 
%Only implemented once per year. Therefore, everything here is considered
%mass balance rather than daily change in ice water equivalent ('icmb' rather than 'icdwe')

%FUNCTIONS TO ONLY BE IMPLEMENTED IF GLACIER OUTLINES PRESENT AND GLACIER
%DYNAMICS ON
if ~strcmpi(sMeta.iceGrid, 'none') && sMeta.glacierDynamics == 1   

    %On 365th day of each year (or when in parameter mode):
    if isequal(sMeta.dateCurr(2:end), sMeta.dateGlac) || regexpbl(sMeta.mode,'parameter')
        %Translates changes calculated over main grid to changes calculated
        %over ice grid:
            %(1) Translate main-grid mass balance to fine-scale grid and
            %(2) Update fractional coverage of main grid (if applicable)
        %Note: This function must be first of glacier process representations 
        %because it ensures ice calculations on main grid from current time
        %step propogate to glacier grid
        if ~regexpbl(sMeta.mode,'parameter')
            mb_main2ice_grid(sHydro, sMeta)    
        end
        
        
        %Glacier velocity:
        gvelMod = find_att(sMeta.module, 'glaciervel', 'no_warning');
        if regexpbl(gvelMod, 'sheer') 
            if regexpbl(sMeta.mode,'parameter')
                coef = cat(1,coef, glaciervel_sheer());
            else
                glaciervel_sheer(sMeta);
            end
        elseif ~regexpbl(gvelMod, 'none') && ~isempty(gvelMod)
            error('cchfModules:glaciervelUnknown', ['The glacier velocity representation ' gvelMod ' not recognized.']);
        end
            
            
        %Glacier sliding:
        glacMoveMod = find_att(sMeta.module,'glaciermove', 'no_warning');
        if ~isempty(glacMoveMod) && regexpbl(glacMoveMod, 'slide')
            %Translate velocity into redistribution of ice between cells
            %(impacts mass balance)
            if regexpbl(sMeta.mode,'parameter')
                coef = cat(1,coef, glaciermove_slide(sHydro));
            else
                glaciermove_slide(sHydro, sMeta);
            end
        elseif ~isempty(glacMoveMod) && (~regexpbl(glacMoveMod, 'static') || ~regexpbl(glacMoveMod, 'none'))
            error('cchfModules:unknwownGlacierSlide', ['The glacier movement representation ' glacMoveMod ' not recognized.']);
        end
    end
     
    
    %Enforce physical constraints
    sCryo.icx( sCryo.icx  < 0) = 0;
    sCryo.icwe(sCryo.icwe < 0) = 0;
end


%At end of each time step, record changes in mass balance at main grid (the main 
%grid mass balance is then translated to the glacier grid once per year):
sCryo.icmb = sCryo.icmb + sCryo.icdwe;
    

%Record end date of glacier glacier mass balance on the last day
%of the year
if isequal(sMeta.dateCurr(2:end), sMeta.dateGlac)
    %Record end date:
    sCryo.icmb_dateEnd = sMeta.dateCurr;
end


%%OUTPUT FOR PARAMETER MODE:
if regexpbl(sMeta.mode,'parameter')
   varargout{1} = coef; 
end