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

function flow_lumped(sHydro, varargin)
%sLand.flow has units of m^3 / time unit used if sLand.mrro has units of m.


global sLand

%If not variable input argument, return parameters and do not proceed.
%Clark time lag:
if isempty(varargin(:))
    return
else
    sMeta = varargin{1};
end

if sum(diag(sHydro.fdr)) > 0
   warning('flow_Clark:endoheroicFlow',['Some cells flow to themselves.  This may '...
       'cause an infinite loop or may not be handled properly.']); 
end


%Find dt (seconds):
dt = time2sec(1,sMeta.dt,sMeta.dateCurr); %Used for distributing previous time-steps streamflow (units = seconds)

if ~isfield(sLand,'totFdrR') || ~isfield(sLand,'totFdrF') || ~regexpbl(find_att(sMeta.module,'timelag'),'johnstone')
    szRun = size(sLand.mrro);
    
    blCalc = 0;
    if regexpbl(find_att(sMeta.module,'timelag'),'johnstone')
        if isfield(sMeta,'runType') && ~regexpbl(sMeta.runType, 'calib')
            if isfield(sMeta,'rtDir')
                rtDir = sMeta.rtDir;
                if iscell(rtDir)
                   rtDir = rtDir{1}; 
                end
            else
                rtDir = pwd;
            end
            rtDir = char(rtDir);
            [rtDir, ~, ~] = fileparts(rtDir);
            
            dirFlow = fullfile(rtDir, 'streamflow_model_grids');
            if ~exist(dirFlow, 'dir')
               mkdir(dirFlow); 
            end

           pathFlow = fullfile(dirFlow, ['module_flow_lumped_array_' sMeta.region{sMeta.siteCurr} '_' sMeta.dt '.mat']);

           %Variable names used for storing runoff and flow grids;
           %If these are changed, it is necessary to change the line where
           %outputs are saved.
           varR = 'runoff';
           varF = 'flow';
        else
            pathFlow = '';
        end
        
        %Boolean value to indicate that time-lag array must be calculated
        if ~isempty(pathFlow) && exist(pathFlow, 'file')
            SLd = load(pathFlow);

            if isfield(SLd, varR) && isfield(SLd, varF)
                sLand.totFdrR = SLd.(varR);
                sLand.totFdrF = SLd.(varF);
            else
                blCalc = 1;
            end
        else
            blCalc = 1;
        end
        
    else
        blCalc = 1;
    end
    
    %Calculate if not already loaded:
    if blCalc == 1
        cumTime  = zeros(szRun ,'single'); %(units = seconds)
        fdrSteps = zeros(szRun,'single');

        [sLand.totFdrR, fdrSteps, cumTime] = fdr_cum(sHydro, sLand.tlag, fdrSteps, cumTime,   dt);
        [sLand.totFdrF,        ~,       ~] = fdr_cum(sHydro, sLand.tlag, fdrSteps, cumTime, 2*dt);
        
        %Save to file
        runoff = sLand.totFdrR;
        flow   = sLand.totFdrF;
        
        if ~isempty(pathFlow)
            save(pathFlow, 'runoff', 'flow');
        end
    end
end


% if isequal(sMeta.dateCurr,[2001, 7,30])
%     keyboard
% end

%Convert NaN values to 0 in case the NaN gets routed downstream.
sLand.mrro(isnan(sLand.mrro)) = 0;

if sMeta.indCurr == 1
    %Add runoff to flow for current time step:
    sLand.flow(:) = (speye(size(sLand.totFdrR)) + sLand.totFdrR)'*double(reshape(sLand.mrro.*sHydro.area,[],1))/dt;
else
    sLand.flow(:) = (speye(size(sLand.totFdrR)) + sLand.totFdrR)'*double(reshape(sLand.mrro.*sHydro.area,[],1))/dt ...
        + sLand.totFdrF'*double(sLand.flow(:));
end
%In FDR: row is cell of origin and col is downstream cell

%Set flow to NaN at all locations where the DEM is NaN:
sLand.flow(isnan(sHydro.dem)) = nan;


%%ROUTE FLOW FROM PREVIOUS TIME-STEP
%
%ONLY ROUTE FLOW FROM CELLS HIGHER UP THAN WOULD HAVE BEEN INCLUDED IN
%ROUTING CURRENT TIMESTEP (uses fdrStep to accomplish this and does not
%reset 'ii')
% if sLand.indCurr ~= 1
% %     %Create indexes for flow from previous time-step:
% %     indFlowFrom = (sLand.indCurr*ones(nGrid,1) - 1) + (rFlow-1)*szFlow(1) + (cFlow-1)*szFlow(2)*szFlow(1);
% %     
%     %Add previous timestep flow to current timestep:
%     sLand.flow(:) = sLand.flow(:) + sLand.totFdrF'*sLand.flow(:);
% end

%Flow from multiple previous timesteps:
% if sLand.indCurr ~= 1
%     totFdr = zeros(size(sHydro.fdr));
% 
%     %First index is that of cell upstream of runoff routing from past
%     %time-step:
%     ii = min2d(fdrSteps) + 1;
%     
%     flagFl = 1; tt = 0; %Don't include flow at each cell from previous time step (but route upstream flow)
%     while flagFl == 1 && ii < nGrid %Loop over time step that the flow is flowing from
%         tt = tt + 1;
%         %Create indexes for flow from previous time-step:
%         indFlowFrom = (sLand.indCurr*ones(nGrid,1) - tt) + (rFlow-1)*szFlow(1) + (cFlow-1)*szFlow(2)*szFlow(1);
%         
%         while flagFl == 1 && ii < nGrid %Create a fdr grid for the current set of flow-from -> flow-to points
%             currFdr = full(sHydro.fdr^ii);   
% 
%             indFdrSkip = find(fdrSteps + 1 > ii);
%             indCurrFdr = setdiff(indRun, indFdrSkip);
%             %Add time lag to cumulative time array
%             tTravCurr = diag(currFdr*sLand.tlag,0);
% 
%             cumTime(indCurrFdr) = cumTime(indCurrFdr) + tTravCurr(indCurrFdr); %(units = seconds)
% 
%     %         indRem = find(cumTime <= dt | cumTime > 2*dt); %Ensure 
% 
%             currFdr(cumTime > tt*dt,:) = 0;
%             currFdr(indFdrSkip,:) = 0;
%             ii = ii + 1;
% 
%             totFdr = totFdr + currFdr;
%             fdrSteps(indCurrFdr) = fdrSteps(indCurrFdr) + 1;
% 
%             %Exit when no more water to route or tiem has been exceeded at all
%                 %points:
%             if all(all(currFdr == 0)) || all(all(cumTime > tt*dt))
%                flagFl = 0; 
%             end
%         end
% 
%         %Add previous timestep flow to current timestep:
%         sLand.flow(indFlowTo) = sLand.flow(indFlowTo) + totFdr'*sLand.flow(indFlowFrom);
%     end
% end