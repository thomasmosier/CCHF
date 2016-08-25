function [varargout] = module_implement(sHydro, sMeta)
global sLand sAtm


%If parameter mode, create output array
if regexpbl(sMeta.mode,'parameter')
   coef = cell(0,6); 
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
if ~regexpbl(find_att(sMeta.module, 'heat', 'no_warning'), {'simple','hock'}) %Albedo not used in simple degree-index model
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
if ~regexpbl(find_att(sMeta.module,'heat', 'no_warning'), 'simple')
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
    else
        error('module_implement:toa_raad','No top of atmosphere radiation process representation was selected.');
    end
end



%TRANSMISSIVITY (not used in simple degree index):
if ~regexpbl(find_att(sMeta.module,'heat', 'no_warning'), 'simple')
    transMod = find_att(sMeta.module,'toa');
    
    if regexpbl(transMod, 'DeWalle')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, atm_transmit_DeWalle());
        else
            atm_transmit_DeWalle(sMeta);
        end
    elseif regexpbl(transMod, 'Coops')
        if regexpbl(sMeta.mode,'parameter')
            coef = cat(1,coef, atm_transmit_Coops(sHydro));
        else
            atm_transmit_Coops(sHydro,sMeta);
        end
%         elseif regexpbl(transMod, 'Gauss')
%             if regexpbl(sMeta.mode,'parameter')
%                 coef = cat(1,coef, atm_transmit_Gauss(sHydro));
%             else
%                 atm_transmit_Gauss(sHydro,sMeta);
%             end
    else
        error('module_implement:transmissivity','No atmospheric transmissivitity process representation was selected.');
    end
end

    

%CRYOSPHERE HEAT:
heatMod = find_att(sMeta.module, 'heat');

if regexpbl(heatMod, 'simple')
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
else
    error('module_implement:cryo','No cryosphere heat process reprepresentation selected.');
end




%CRYOSPHERE ACCUMULATION AND MELT:
massMod = find_att(sMeta.module,'mass');

if regexpbl(massMod, 'step') %Calculate snow and glacier melt using heat threshold
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, mass_step_tas());
    else
        mass_step_tas(sMeta);
    end
elseif regexpbl(massMod, 'cc')
    if regexpbl(sMeta.mode,'parameter') %Calculate snow and melt using CC
        coef = cat(1,coef, mass_cc_v2());
    else
        mass_cc_v2(sMeta);
    end
elseif regexpbl(massMod, {'enbal', 'energy'}) %Conduction raises and lowers snow temperature
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, mass_enbal());
    else
        mass_enbal(sMeta);
    end
elseif regexpbl(massMod, 'Liston') %Conduction raises and lowers snow temperature (no fitting parameters)
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, mass_Liston());
    else
        mass_Liston(sMeta);
    end
else
    error('module_implement:cryo','No cryosphere energy and mass process representation selected.');
end



%SNOW LIQUID WATER HOLDING CAPACITY:
liqMod = find_att(sMeta.module,'snlq', 'no_warning');

if regexpbl(liqMod, 'percent') %Release melt water greater than percent holding capacity:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, snlq_percent());
    else
        snlq_percent(sMeta);
    end
elseif regexpbl(liqMod, 'density') %Release melt water greater than density threshold:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, snlq_density());
    else
        snlq_density(sMeta);
    end
elseif isempty(liqMod) || regexpbl(liqMod, 'zero') %Release all melt water:
    if regexpbl(sMeta.mode,'parameter')
        coef = cat(1,coef, snlq_zero());
    else
        snlq_zero(sMeta);
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
        coef = cat(1,coef, runoff_bucket(sHydro));
    else
        runoff_bucket(sHydro,sMeta);
    end
elseif regexpbl(runoffMod, 'direct') %All water (precip and snowpack release) runs off immediately
    if regexpbl(sMeta.mode, 'parameter')
        coef = cat(1,coef, runoff_direct(sHydro));
    else
        runoff_direct(sHydro,sMeta);
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
    else
        error('module_implement:routing','No flow routing process representation selected.');
    end
end
    
    

%%OUTPUT FOR PARAMETER MODE:
if regexpbl(sMeta.mode,'parameter')
   varargout{1} = coef; 
end