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
global sLand sAtm


%If parameter mode, create output array
if regexpbl(sMeta.mode,'parameter')
   coef = cell(0,6); 
end


%Initialize ice thickness on first iteration:
if ~regexpbl(sMeta.mode,'parameter') && sMeta.indCurr == 1
    iceWEMod = find_att(sMeta.module, 'glacier0');
    
    if regexpbl(iceWEMod, 'Shea')
        %Estimate thickness (force balance based on equilibrium
        %assumption from Shea et al. 2015)
        glacier0_Shea(sHydro, sMeta);
    elseif regexpbl(iceWEMod, 'Chen')
        %Estimate thickness (empirical surface area scaling based on Chen and Ohmura, 1990)
        glacier0_Chen(sHydro, sMeta);
    else
        error('module_implement:iceWE','No ice water equivalent process representation selected.');
    end 
end

%%IMPLEMENT FOR ALL MODULE COMBINATIONS:
%PARTITION RAIN AND SNOWFALL:
%No Choices
if regexpbl(sMeta.mode,'parameter')
    coef = cat(1,coef,part_snow());
else
    part_snow(sMeta); %Partitions precipitation to snowfall and rain.
    snow_accum; %Adds snowfall to snowpack and rain to liquid snow content
end

  



%ALBEDO:
if ~regexpbl(find_att(sMeta.module, 'heat', 'no_warning'), {'simple','hock', 'SDI'}) %Albedo not used in simple degree-index model
    albMod = find_att(sMeta.module, 'albedo');
    
    if regexpbl(albMod, 'Pelli')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, albedo_Pellicciotti());
        else
            albedo_Pellicciotti(sMeta);
        end
    elseif regexpbl(albMod, 'Brock')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, albedo_Brock());
        else
            albedo_Brock(sMeta);
        end
    else
        error('module_implement:albedo','No albedo process representation selected.');
    end


    %If debris grid included, set albedo to 0.15 at debris covered ice locations:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, albedo_debris());
    else
        albedo_debris(sMeta);
    end
end



%TOP OF ATMOSPHERE SOLAR RADIATION (not used in simple degree index):
if ~regexpbl(find_att(sMeta.module,'heat', 'no_warning'), {'simple', 'SDI'})
    toaMod = find_att(sMeta.module,'toa');
    
    if regexpbl(toaMod, 'DeWalle')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, toa_rad_DeWalle([], sHydro));
        else
            if ~isfield(sAtm,'rsdt')
                toa_rad_DeWalle(sMeta.dateRun, sHydro, sMeta);
            end
        end
    elseif regexpbl(toaMod,'Liston')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, toa_rad_Liston([], sHydro));
        else
            if ~isfield(sAtm,'rsdt')
                toa_rad_Liston(sMeta.dateRun, sHydro, sMeta);
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
        error('module_implement:toa_raad','No top of atmosphere radiation process representation was selected.');
    end
end



%TRANSMISSIVITY (not used in simple degree index):
if ~regexpbl(find_att(sMeta.module,'heat', 'no_warning'), {'simple', 'SDI'})
    transMod = find_att(sMeta.module,'atmtrans');
    
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
    elseif regexpbl(transMod, 'dem_exp_decay')
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
        error('module_implement:transmissivity','No atmospheric transmissivitity process representation was selected.');
    end
end

    

%CRYOSPHERE HEAT:
heatMod = find_att(sMeta.module, 'heat');

if regexpbl(heatMod, 'SDI')
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
elseif regexpbl(heatMod, {'simple','debris'},'and')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, heat_simple_debris());
    else
        heat_simple_debris(sMeta)
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
    error('module_implement:cryo','No cryosphere heat process reprepresentation selected.');
end



%SNOWPACK DENSITY:
%Comes after heat flux because depends on cryosphere surface temperature
if regexpbl(find_att(sMeta.module, 'heat', 'no_warning'), 'SETI') || regexpbl(find_att(sMeta.module, 'snlq', 'no_warning'), 'Liston')%Albedo not used in simple degree-index model
    dnsMod = find_att(sMeta.module, 'density');
    
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
    else
        error('module_implement:albedo','No density process representation selected.');
    end
end



%SNOW ACCUMULATION AND MELT:
massMod = find_att(sMeta.module, 'snmass');

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
    error('module_implement:cryo','No cryosphere energy and mass process representation selected.');
end


%SNOW SUBLIMATION:
if sMeta.indCurr == 1
    sublimateMod = find_att(sMeta.module,'sublimate');
else
    sublimateMod = find_att(sMeta.module,'sublimate', 'no_warning');
end

if regexpbl(sublimateMod, 'Lutz')
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, sublimate_Lutz(sHydro));
    else
        sublimate_Lutz(sHydro, sMeta);
    end
elseif ~isempty(sublimateMod)
    error('module_implement:unknownSublimate', [sublimateMod ' is an ' ...
        'unknown sublimation module that has not been programmed for.']);
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
elseif isempty(liqMod) || regexpbl(liqMod, 'zero') %Release all melt water:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, sndrain_zero());
    else
        sndrain_zero(sMeta);
    end
else
    error('module_implement:unknwownLiqFun',['The liquid holding capacity option is ' char(39) ...
        liqMod char(39) ', which is an unknown choice.']);
end



%ICE MELT AND LIQUID RELEASE:
icMelt = find_att(sMeta.module,'icmlt', 'no_warning');

if regexpbl(icMelt,'ratio') %Remainder of energy from snow melt goes into melting ice
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_ratio());
    else
        icmlt_ratio(sMeta);
    end
elseif regexpbl(icMelt,'Liston') %Ice melt only occurs during time steps when no snow present
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_Liston());
    else
        icmlt_Liston(sMeta);
    end
elseif regexpbl(icMelt,'Neumann') %Ice melt only occurs during time steps when no snow present
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_Neumann_bc());
    else
        icmlt_Neumann_bc(sMeta);
    end
elseif regexpbl(icMelt,'time') %Ice melt only occurs during time steps when no snow present
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_time());
    else
        icmlt_time(sMeta);
    end
elseif regexpbl(massMod, 'gradient') %Conduction raises and lowers snow temperature (no fitting parameters)
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, icmlt_ablation_grad(sHydro));
    else
        icmlt_ablation_grad(sHydro, sMeta);
    end
else
    error('module_implement:unknwownIceMltFun',['The ice melt option is ' ...
        char(39) icMelt char(39) ', which is an unknown choice.']);
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
    error('module_implement:unknwownSicaFun',['The snow/ice covered area option is ' char(39) ...
        scaMod char(39) ', which is an unknown choice.']);
end   



%PET:
petMod = find_att(sMeta.module,'pet');

if regexpbl(petMod, 'Hammon') %Use Hammon formulation
    if regexpbl(sMeta.mode, 'parameter')
        coef = cat(1,coef, PET_Hammon(sMeta.dateRun, sHydro));
    else
        if ~isfield(sLand, 'pet') || ~isfield(sLand, 'datepet')
            PET_Hammon(sMeta.dateRun, sHydro, sMeta);
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
    warning('module_implement:PET','No potential evapotranspiration process representation selected.');
end



%RUNOFF:
runoffMod = find_att(sMeta.module, 'runoff');

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
    error('module_implement:runoff','No runoff process representation selected.');
end



%TRAVEL TIME THROUGH CELLS:
timeMod = find_att(sMeta.module,'timelag');

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
else
    error('module_implement:tLag','No surface water time-lag process representation selected.');
end



%ROUTING RUNOFF AND FLOW:
flowMod = find_att(sMeta.module,'flow');

if ~regexpbl(sMeta.mode,'parameter') %Use numerical Muskingum method (allows dispersion)
    if regexpbl(flowMod, 'Muskingum')
        flow_Muskingum(sHydro, sMeta);
    elseif regexpbl(flowMod, 'lumped') %No dispersion of flow
        flow_lumped(sHydro, sMeta);
    elseif regexpbl(flowMod, 'Liston') %Fast and slow flow
        flow_Liston(sHydro, sMeta);
    else
        error('module_implement:routing','No flow routing process representation selected.');
    end
end


%GLACIER PROCESSES (INCLUDES GLACIER DYNAMICS OPTIONS): 
%FUNCTIONS TO ONLY BE IMPLEMENTED IF GLACIER OUTLINES PRESENT
if isfield(sHydro,'icbl') && any2d(sHydro.icbl) && sMeta.glacierDynamics == 1
    if sMeta.indCurr == 1 && ~regexpbl(sMeta.mode, 'calib')
        disp('Calculating glacier dynamics')
    end
    glacMove = find_att(sMeta.module,'glaciermove', 'no_warning');
    
    firn = find_att(sMeta.module,'firn', 'no_warning');
            
    if ~regexpbl(sMeta.mode,'parameter') && regexpbl(firn, 'simple')
        %Model firn compaction (i.e. metamorphism from snow to ice) using
        %simple threshold
        %This must come before 'ice_link'
        firn_compact_simple();
    elseif ~regexpbl(sMeta.mode,'parameter') && ~regexpbl(firn, 'static')
        error('module_implement:glacierDepth','Unknown glacier depth module selected.');
    end
    
    %Translates changes calculated over main grid to changes calculated
    %over ice grid:
    if ~regexpbl(sMeta.mode,'parameter')
        ice_link(sHydro)    %This function must be first of glacier process 
                            %representations because it zeros glacier 
                            %change for current time-step

        %Model avalanching of snow (using glacier grid):
        avalanche(sHydro, sMeta)
    end

    %Glacier ice velocity
    if ~isempty(glacMove) && regexpbl(glacMove, 'Shea')
        %Model ice velocity:
        glacier_velocity(sMeta);

        %Translate velocity into redistribution of ice between cells
        glacier_slide(sMeta);
    elseif ~isempty(glacMove) && ~regexpbl(glacMove, 'static')
        error('module_implement:unknwownGlacierSlide',['The optional glacier process module is ' ...
            glacMove ', which is an unknown choice.']);
    end
end
    
    

%%OUTPUT FOR PARAMETER MODE:
if regexpbl(sMeta.mode,'parameter')
   varargout{1} = coef; 
end